function [X,F,D,A,sfpts] = randomdata(numMSp,numFac,mzs,jitter1,jitter2,mods,scans, peakstd1,peakstd2,smpls,SNR)
%Generates a tensor of random values, with drift in the first and second
%retention modes. %Let's analyze it in the same function.

%example settings
%numMSp = 20;
%numFac = 3;
%mzs = 761;
%jitter1 = 1.5; %1.5
%jitter2 = 25; %25
%mods = 20;
%scans = 200;
%peakstd1 = 0.75; %0.75
%peakstd2 = 20; %20
%smpls = 8;
%SNR = 50;

A = zeros(mzs,numFac);

for aa = 1:numFac
    
    for pp = 1:numMSp
        A(:,aa) = A(:,aa) +  rand(1,1).*normpdf(1:1:761,761*rand(1,1),0.5)';
    end
 
    
    A(:,aa) = A(:,aa)./norm(A(:,aa));
end

%Generate the nominal retention times using rand
    tr12o = [(mods/3 + ((mods - mods/3) - mods/3).*rand(numFac,1)), (scans/3 + ((scans-scans/3) - scans/3).*rand(numFac,1))];

for qq = 1:smpls 

    

    tr1pm = (-jitter1 + (jitter1+jitter1).*rand(numFac,1));
    tr2pm = (-jitter2 + (jitter2+jitter2).*rand(numFac,1));
    pm = [tr1pm,tr2pm];
    
    tr12 = tr12o + pm; 
    [F(:,:,qq),D(:,:,qq),A] = gc2syn([scans,mods],tr12(:,1),repelem(peakstd1,numFac),repmat(tr12(:,2),[1,mods]),repelem(peakstd2,numFac),repelem(SNR,numFac),A);
    
    for aa = 1:numFac
        plotF(:,:,aa) = reshape(F(:,aa,qq),[scans,mods]);
        pts{aa} = plotF(:,:,aa);
    end
    
    sf = surfacfus(pts{:});
    sfpts(:,:,:,qq) = sf;
    
    Xtemp(:,:,qq) = (F(:,:,qq)*D(:,:,qq)*A' + (1/SNR)*sqrt(var(F(:,:,qq)*D(:,:,qq)*A',0,'all')).*randn(scans*mods,mzs));
    Xtemp2 = reshape(Xtemp(:,:,qq),[scans,mods,mzs]);
    X(:,:,:,qq) = permute(Xtemp2,[1,3,2]);

    X(X<0) = 0;
        
end



end

function [F,Dk,A] = gc2syn(dim,tr1s,tr1sig,tr2s,tr2sig,SNRs,mzs)
%Function to generate synthetic GCxGC-TOFMS Datasets. Call this function
%once per sample.

%X = rand(dim(1),dim(2),size(mzs,1));

%Generate the first dimension retention times


for ww = 1:size(mzs,2)
   
    ytr1(ww,:) = normpdf(1:1:dim(2),tr1s(ww),tr1sig(ww))';
    
for ll = 1:dim(2)
    
    ytr2(:,ll) = ytr1(ww,ll).*normpdf(1:1:dim(1),tr2s(ww,ll),tr2sig(ww));
    
end

F(:,ww) = ytr2(:);
F(:,ww) = F(:,ww)./norm(F(:,ww));
F(:,ww) = F(:,ww).*SNRs(ww);
Dk(ww,ww) = norm(F(:,ww));
%F(:,ww) = F(:,ww) + rand(size(F(:,ww),1),1); %Add noise.
%F(:,ww) = F(:,ww) + (max(max(max(F(:,ww))))/SNRs(ww)).*randn(size(F,1),1);
%F(F<0) = 0;
F(:,ww) = F(:,ww)./norm(F(:,ww));

A(:,ww) = mzs(:,ww)./norm(mzs(:,ww));

end

end