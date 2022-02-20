function [F,A,D,pvar,SSR] = nnparafac2x2seq(Xijkl,R)
%nnPARAFAC2x2
%
%Algorithm for linear tensor decomposition of GCxGC-TOFMS data, based on
%alternating 2tr, 1tr frontal slabs to solve for score tensor determination, which with
%the mass spec is then used to solve for the sample loadings. Agnostic to
%drift in the first and second dimensions.CIX.
%
%(c) Michael Sorochan Armstrong, James James Harynuk 2021
%
%

%Convergence properties
eps = 1e-5;
iter = 1;
maxiter = 10000;
coupit = 10;
nnls = 2; %1 for fnnls, 2 for fcnnls, 3 for cls
animate = 1;
omg = 2;
YNorm = norm(Xijkl(:))^2;

sz = size(Xijkl);

%2tr frontal plane - tensor of observations
%Xkl = reshape(Xijkl,[sz(1),sz(2),sz(3)*sz(4)]);
Xkl = reshape(Xijkl,[sz(1),sz(2),sz(3)*sz(4)]);

%1tr frontal plane - tensor of observations
%Xil = reshape(Xijkl,[sz(3),sz(2),sz(1)*sz(4)]); %reshape is broken here.
Xilp = permute(Xijkl,[3 2 1 4]);
Xil = reshape(Xilp,[sz(3),sz(2),sz(1)*sz(4)]);

%initialising values, randomised (will use ICA in the future)
% A = zeros(sz(2),R);
% for kl = 1:sz(3)*sz(4)
%     [~,~,V] = svds(Xkl(:,:,kl),R);
%     A(:,:) = A(:,:) + V;
% end
% 
% for rr = 1:R
%     A(:,rr) = A(:,rr)./norm(A(:,rr));
% end

%Akl = abs(orth(rand(sz(2),R)));
%Ail = abs(orth(rand(sz(2),R)));

Akl = rand(sz(2),R);
for ii = 1:R
    Akl(:,ii) = Akl(:,ii)./norm(Akl(:,ii));
end
%Ail = rand(sz(2),R);
Ail = Akl;

F1 = diag(ones(1,R));
Dkl = repmat(F1,[1,1,sz(3)*sz(4)]);
Bkl = rand(sz(1),R,sz(3)*sz(4));
Bskl = rand(R,R);
for rr = 1:R
    Bskl(:,rr) = Bskl(:,rr)./norm(Bskl(:,rr));
end

% for kl = 1:sz(3)*sz(4)
%     %Xhkl = Bkl(:,:,kl)*Dkl(:,:,kl)*A';
%     mkl(kl) = sum(sum(sum(Bkl(:,:,kl)*Dkl(:,:,kl)*Akl')))./sum(sum(sum(Bkl(:,:,kl))));
% end

for kl = 1:sz(3)*sz(4)
    [U,~,V] = svds(Bkl(:,:,kl)*Bskl,R);
    Pkl(:,:,kl) = U(:,:)*V(:,:)';
end

for kl = 1:sz(3)*sz(4)
   mkl(kl) = norm(Bkl(:,:,kl)*Dkl(:,:,kl)*Akl','fro')./norm(Bkl(:,:,kl),'fro'); 
end

Bil = rand(sz(3),R,sz(1)*sz(4));
Dil = repmat(F1,[1,1,sz(1)*sz(4)]);
Bsil = rand(R,R);

for il = 1:sz(1)*sz(4)
    [U,~,V] = svds(Bil(:,:,il)*Bsil,R);
    Pil(:,:,il) = U(:,:)*V(:,:)';
end


for il = 1:sz(1)*sz(4)
    mil(il) = norm(Bil(:,:,il)*Dil(:,:,il)*Ail','fro')./norm(Bil(:,:,il),'fro');
end

for rr = 1:R
    Bsil(:,rr) = Bsil(:,rr)./norm(Bsil(:,rr));
end

for il = 1:sz(1)*sz(4)
    XHil(il) = norm(Xil(:,:,il) - Bil(:,:,il)*(Dil(:,:,il))*Ail','fro')^2;
    mhil(il) = mil(il)*norm(Bil(:,:,il) - Pil(:,:,il)*Bsil,'fro')^2; 
end

%mA = 0.033; %mass spectral coupling constant

for kl = 1:sz(3)*sz(4)
    XHkl(kl) = norm(Xkl(:,:,kl) - Bkl(:,:,kl)*(Dkl(:,:,kl))*Akl','fro')^2; 
    mhkl(kl) = mkl(kl)*norm(Bkl(:,:,kl) - Pkl(:,:,kl)*Bskl,'fro')^2;
end

mA = (10^(omg))*(sum(XHkl) +sum(XHil));

mha = mA*norm(Akl - Ail,'fro')^2;

%ssr1 = sum(sum(sum(XHkl))) + sum(sum(sum(XHil)));
ssr1 = sum(XHkl + mhkl) + sum(XHil + mhil) + mha;

ssr2 = 1;

%Ckl = ones(sz(3)*sz(4),R);
%Cil = ones(sz(1)*sz(4),R);

if animate == 1
    h = figure;
    axis tight manual
    filename = strcat(date,'msdescent.gif');
else
end

while abs(ssr1 - ssr2) > abs(ssr2)*eps && abs(ssr1 - ssr2)/abs(ssr2) > eps && iter < maxiter
    
    ssr1 = ssr2;
    
    %PARAFAC2 - Xkl
    KL = sz(3)*sz(4);
    
    %Pkl Estimation
    if iter > 1
    for kl = 1:KL
        [U,~,V] = svds(Bkl(:,:,kl)*Bskl,R);
        Pkl(:,:,kl) = U(:,:)*V(:,:)';
    end
    else
    end
    
    %Bskl estimation
    for kl = 1:KL
        Bstkl(:,:,kl) = mkl(kl).*Pkl(:,:,kl)'*Bkl(:,:,kl);
    end
    
    Bskl = 1/(sum(mkl))*(sum(Bstkl,3));
    
    for rr = 1:R
        Bskl(:,rr) = Bskl(:,rr)./norm(Bskl(:,rr));
    end
    
    %Akl Estimation
    for kl = 1:sz(3)*sz(4)
        BklDkl(:,:,kl) = Bkl(:,:,kl)*Dkl(:,:,kl);
    end

    BklDklIK = nshape(BklDkl,2)';

    Xjkil = nshape(Xkl,2)';

    %for kl = 1:sz(3)*sz(4)
        %Akl = fcnnls([],[],BklDklIK'*BklDklIK,BklDklIK'*Xjkil - mA.*Ail')'; %Removed +mA.*Akl'
        %Akl = fcnnls([],[],BklDklIK'*BklDklIK + mA.*eye(R),(mA.*Ail +
        %Xjkil'*BklDklIK)')'; %%PROBLEM
    %end

        if nnls == 1
            for jj = 1:sz(2)
        %Ail = fcnnls([],[],BilDilIK'*BilDilIK + mA.*eye(R),(mA.*Akl + Xjilk'*BilDilIK)')';
                Akl(jj,:) = fnnls(BklDklIK'*BklDklIK + mA.*eye(R),(mA.*Ail(jj,:) + Xjkil(:,jj)'*BklDklIK)');
            end
            elseif nnls == 2
                Akl = fcnnls([],[],BklDklIK'*BklDklIK + mA.*eye(R),(mA.*Ail + Xjkil'*BklDklIK)')';
            elseif nnls == 3
                Akl = pinv(BklDklIK'*BklDklIK + mA.*eye(R))*(mA.*Ail + Xjkil'*BklDklIK)';
        end
            
%     if any(sum(Akl,2) == 0) %|| any(isnan(sum(Akl,2)))
%         Akl = 0.9*(rand(size(Akl,1),size(Akl,2))) + 0.1*Akl;
%     end

    for rr = 1:R
        Akl(Akl == 0) = 1e-20;
        Akl(:,rr) = Akl(:,rr)./norm(Akl(:,rr));
    end
    
    %Bkl Estimation
    for kl = 1:KL
        %for ii  = 1:sz(1)
        %    Bkl(ii,:,kl) = fnnls((Dkl(:,:,kl)*(A'*A)*Dkl(:,:,kl)*eye(R)),(Xkl(ii,:,kl)*A*Dkl(:,:,kl) + mkl(kl)*Pkl(ii,:,kl)*Bskl)');
        %end
        if nnls == 1
            for ii = 1:sz(1)
            Bkl(ii,:,kl) = fnnls((Dkl(:,:,kl)*(Akl'*Akl)*Dkl(:,:,kl) + mkl(kl)*eye(R)),(Xkl(ii,:,kl)*Akl*Dkl(:,:,kl) + mkl(kl)*Pkl(ii,:,kl)*Bskl)');
            end
        elseif nnls == 2
        Bkl(:,:,kl) = fcnnls([],[],(Dkl(:,:,kl)*(Akl'*Akl)*Dkl(:,:,kl) + mkl(kl)*eye(R)),(Xkl(:,:,kl)*Akl*Dkl(:,:,kl) + mkl(kl)*Pkl(:,:,kl)*Bskl)')';
        elseif nnls == 3
        Bkl(:,:,kl) = (pinv(Dkl(:,:,kl)*(Akl'*Akl)*Dkl(:,:,kl) + mkl(kl)*eye(R))*(Xkl(:,:,kl)*Akl*Dkl(:,:,kl) + mkl(kl)*Pkl(:,:,kl)*Bskl)')';
        end
    end
   
    for kl = 1:KL
        for rr = 1:R
            Btkl = Bkl(:,rr,kl);
            Btkl(Btkl == 0) = 1e-20;
            Bkl(:,rr,kl) = Btkl;
            Bkl(:,rr,kl) = Bkl(:,rr,kl)./norm(Bkl(:,rr,kl));
        end
    end
    
    %Dkl Estimation
    for kl = 1:KL
        %Dtkl = diag((pinv(Bkl(:,:,kl)'*Bkl(:,:,kl))*Bkl(:,:,kl)'*Xkl(:,:,kl))/Akl');
        %Dkl(:,:,kl) = diag(Dtkl);
        %Dkl(:,:,kl) = diag(diag(pinv(Bkl(:,:,kl))*Xkl(:,:,kl)*pinv(Akl'))');
        %Dkl(:,:,kl) = diag(diag((Bkl(:,:,kl)'*Xkl(:,:,kl)*Akl)*pinv((Bkl(:,:,kl)'*Bkl(:,:,kl)).*(Akl'*Akl))));
        %Dkl(:,:,kl) = diag(diag((pinv(Bkl(:,:,kl))*Xkl(:,:,kl))/Akl')');
        %%Cohen
        %Dkl(:,:,kl) = diag(diag((Akl'*Xkl(:,:,kl)'*Bkl(:,:,kl))*pinv((Bkl(:,:,kl)'*Bkl(:,:,kl))*(Akl'*Akl)))); %Calculated by me
        Dkl(:,:,kl) = diag(diag(pinv(Akl'*Akl)*(Akl'*Xkl(:,:,kl)'*Bkl(:,:,kl))*pinv(Bkl(:,:,kl)'*Bkl(:,:,kl))));
    end
    
    %mkl Estimation
    if iter == 1
        for kl = 1:sz(3)*sz(4)
            S = svds(Xkl(:,:,kl),2);% - mean(Xkl(:,:,kl)));
            SNRkl = S(1)/S(2);
            %SNRkl = 10;
            mkl(kl) = 10^(-SNRkl/10)*sum(sum(sum(sqrt((Xkl(:,:,kl) - Bkl(:,:,kl)*Dkl(:,:,kl)*Akl').^2))))/sum(sum(sum(sqrt((Bkl(:,:,kl) - Pkl(:,:,kl)*Bskl).^2))));
        end
    else
        for kl = 1:sz(3)*sz(4)
            if iter < coupit
                mkl(kl) = min(mkl(kl)*1.02,1e12);
            end
        end
    end
    
    
    %PARAFAC2 - Xil
    IL = sz(1)*sz(4);
    
    %Pil Estimation
    if iter > 1
    for il = 1:IL
        [U,~,V] = svds(Bil(:,:,il)*Bsil,R);
        Pil(:,:,il) = U(:,:)*V(:,:)';
    end
    else
    end
    
    %Bsil estimation
    for il = 1:IL
        Bstil(:,:,il) = mil(il).*Pil(:,:,il)'*Bil(:,:,il);
    end
    
    Bsil = 1/(sum(mil))*(sum(Bstil,3));
    
    for rr = 1:R
        Bsil(:,rr) = Bsil(:,rr)./norm(Bsil(:,rr));
    end
    
    %Ail estimation
    for il = 1:sz(1)*sz(4)
        BilDil(:,:,il) = Bil(:,:,il)*Dil(:,:,il);
    end

    BilDilIK = nshape(BilDil,2)';

    Xjilk = nshape(Xil,2)';

    %for ki = 1:sz(2)
        %Ail = fcnnls([],[],BilDilIK'*BilDilIK,BilDilIK'*Xjilk - mA.*Akl')'; %Removed +mA.*Ail'
        if nnls == 1
            for jj = 1:sz(2)
        %Ail = fcnnls([],[],BilDilIK'*BilDilIK + mA.*eye(R),(mA.*Akl + Xjilk'*BilDilIK)')';
                Ail(jj,:) = fnnls(BilDilIK'*BilDilIK + mA.*eye(R),(mA.*Akl(jj,:) + Xjilk(:,jj)'*BilDilIK)');
            end
            elseif nnls == 2
                Ail = fcnnls([],[],BilDilIK'*BilDilIK + mA.*eye(R),(mA.*Akl + Xjilk'*BilDilIK)')';
            elseif nnls == 3
                Ail = pinv(BilDilIK'*BilDilIK + mA.*eye(R))*(mA.*Akl + Xjilk'*BilDilIK)';
        end
            

    %end

    %if any(sum(Ail1,2) == 0) %&& iter > 1
        %Ail1 = 0.1*Ail1 + 0.9.*A1;
    %end

    for rr = 1:R
        Ail(Ail == 0) = 1e-20;
        Ail(:,rr) = Ail(:,rr)./norm(Ail(:,rr));
    end
    
%     %Bkl Estimation
%     for kl = 1:KL
%         %for ii  = 1:sz(1)
%         %    Bkl(ii,:,kl) = fnnls((Dkl(:,:,kl)*(A'*A)*Dkl(:,:,kl)*eye(R)),(Xkl(ii,:,kl)*A*Dkl(:,:,kl) + mkl(kl)*Pkl(ii,:,kl)*Bskl)');
%         %end
%         
%         Bkl(:,:,kl) = fcnnls([],[],(Dkl(:,:,kl)*(A'*A)*Dkl(:,:,kl)*eye(R)),(Xkl(:,:,kl)*A*Dkl(:,:,kl) + mkl(kl)*Pkl(:,:,kl)*Bskl)')';
%     end
%    
%     for kl = 1:KL
%         for rr = 1:R
%             Btkl = Bkl(:,rr,kl);
%             Btkl(Btkl == 0) = 1e-6;
%             Bkl(:,rr,kl) = Btkl;
%             Bkl(:,rr,kl) = Bkl(:,rr,kl)./norm(Bkl(:,rr,kl));
%         end
%     end
    
    %Bil Estimation
    for il = 1:IL
        %for kk  = 1:sz(3)
            %Bil(kk,:,il) = fnnls([],[],(Dil(:,:,il)*(A'*A)*Dil(:,:,il)*eye(R)),(Xil(kk,:,il)*A*Dil(:,:,il) + mil(il)*Pil(kk,:,il)*Bsil)');
        %end
        
        if nnls == 1
            for kk = 1:sz(3)
            Bil(kk,:,il) = fnnls((Dil(:,:,il)*(Ail'*Ail)*Dil(:,:,il) + mil(il)*eye(R)),(Xil(kk,:,il)*Ail*Dil(:,:,il) + mil(il)*Pil(kk,:,il)*Bsil)');
            end
        elseif nnls == 2
        Bil(:,:,il) = fcnnls([],[],(Dil(:,:,il)*(Ail'*Ail)*Dil(:,:,il) + mil(il)*eye(R)),(Xil(:,:,il)*Ail*Dil(:,:,il) + mil(il)*Pil(:,:,il)*Bsil)')';
        elseif nnls == 3 %problem
        Bil(:,:,il) = pinv(Dil(:,:,il)*(Ail'*Ail)*Dil(:,:,il) + mil(il)*eye(R))*(Xil(:,:,il)*Ail*Dil(:,:,il) + mil(il)*Pil(:,:,il)*Bsil)';
        end
        
    end
    
    for il = 1:IL
        for rr = 1:R
            Btil = Bil(:,rr,il);
            Btil(Btil == 0) = 1e-20;
            Bil(:,rr,il) = Btil;
            Bil(:,rr,il) = Bil(:,rr,il)./norm(Bil(:,rr,il));
        end
    end
            
    %Dil Estimation
    for il = 1:sz(1)*sz(4)
        %Dtil = diag((pinv(Bil(:,:,il)'*Bil(:,:,il))*Bil(:,:,il)'*Xil(:,:,il))/Ail');
        %Dil(:,:,il) = diag(fnnls(Bil(:,:,il)'*Bil(:,:,il),Bil(:,:,il)'*Xil(:,:,il)/Ail'));
        %Dil(:,:,il) = diag(Dtil);
        %Dil(:,:,il) = diag(diag(pinv(Bil(:,:,il))*Xil(:,:,il)*pinv(Ail'))');
        %Dil(:,:,il) = diag(diag((Bil(:,:,il)'*Xil(:,:,il)*Ail)*pinv((Bil(:,:,il)'*Bil(:,:,il)).*(Ail'*Ail))));
        %Dil(:,:,il) = diag(diag((Bil(:,:,il)\Xil(:,:,il))/(Ail'))');
        %Dil(:,:,il) = diag(diag((Ail'*Xil(:,:,il)'*Bil(:,:,il))*pinv((Bil(:,:,il)'*Bil(:,:,il))*(Ail'*Ail))));
        %Dil(:,:,il) = diag(diag(pinv(Bil(:,:,il)'*Bil(:,:,il))*(Ail'*Xil(:,:,il)'*Bil(:,:,il))*pinv(Ail'*Ail)));
        %Dil(:,:,il) = diag(diag((pinv(Bil(:,:,il))*Xil(:,:,il))/Ail')');
        Dil(:,:,il) = diag(diag(pinv(Ail'*Ail)*(Ail'*Xil(:,:,il)'*Bil(:,:,il))*pinv(Bil(:,:,il)'*Bil(:,:,il))));
    end
    
    %mil Estimation
    if iter == 1
        for il = 1:sz(1)*sz(4)
            S = svds(Xil(:,:,il),2);% - mean(Xil(:,:,il)));
            %SNRil = sum(diag(S(1,1)))/S((2),(2));
            SNRil = S(1)/S(2);
            %SNRil = 10;
            mil(il) = 10^(-SNRil/10)*sum(sum(sum(sqrt((Xil(:,:,il) - Bil(:,:,il)*Dil(:,:,il)*Ail').^2))))/sum(sum(sum(sqrt((Bil(:,:,il) - Pil(:,:,il)*Bsil).^2))));
        end
    else
        for il = 1:sz(1)*sz(4)
            if iter < coupit
                mil(il) = min(mil(il)*1.02,1e12);
            end
        end
    end
    
    %A Estimation hmm not so sure
    %for jj = 1:sz(2)
    %    A1(jj,:) = fnnls(sum(Dkl,3)*(sum(Bkl,3)'*sum(Bkl,3))*sum(Dkl,3) + sum(Dil,3)*(sum(Bil,3)'*sum(Bil,3))*sum(Dil,3),(sum(Xkl(:,jj,:),3)'*sum(Bkl,3)*sum(Dkl,3)+sum(Xil(:,jj,:),3)'*sum(Bil,3)*sum(Dil,3))')';
    %end
    
    
    %This is way too slow - we must do something to speed it up. Summing
    %every 1tr sample-wise, because there can be no drift there.
    
    %for ll = 1:sz(4)            
        
    %    XilS(:,:,ll) = sum(Xil(:,:,ll*sz(1):(ll*sz(1)+sz(1))),3);
    %    DilS(:,:,ll) = sum(Dil(:,:,ll*sz(1):(ll*sz(1)+sz(1))),3);
    %    BilS(:,:,ll) = sum(Bil(:,:,ll*sz(1):(ll*sz(1)+sz(1))),3);
    %end
    
    %Problematic area
    %XilS = permute(Xil,[3 1 2]);
    %XilS = reshape(XilS(1:sz(1)*sz(4),:,:),sz(4),sz(1),sz(3),sz(2));
    %XilS = sum(XilS,2);
    %XilS = permute(XilS,[3 4 1 2]);
    
    %DilS = permute(Dil,[3 1 2]);
    %DilS = reshape(DilS(1:sz(1)*sz(4),:,:),sz(4),sz(1),R,R);
    %DilS = sum(DilS,2);
    %DilS = permute(DilS,[3 4 1 2]);
    
    %BilS = permute(Bil,[3 1 2]);
    %BilS = reshape(BilS(1:sz(1)*sz(4),:,:),sz(4),sz(1),sz(3),R);
    %BilS = sum(BilS,2);
    %BilS = permute(BilS,[3 4 1 2]);
    
%     x1 = 1:sz(1):sz(1)*sz(4);
%     x2 = sz(1):sz(1):sz(1)*sz(4);
%     
%     %This works, summing into sample-wise matrices of 1tr x mods to reduce
%     %the number of mass spectra to calculate.
%     for ll = 1:sz(4)
%         XilS(:,:,ll) = sum(Xil(:,:,x1(ll):x2(ll)),3);
%         DilS(:,:,ll) = sum(Dil(:,:,x1(ll):x2(ll)),3);
%         BilS(:,:,ll) = sum(Bil(:,:,x1(ll):x2(ll)),3);
%     end
    
    %This is the serial loop that works
%     for kl = 1:sz(3)*sz(4)
%         for ll = 1:sz(4)
%             for jj = 1:sz(2)
%                 A1(jj,:,kl,ll) = fnnls(Dkl(:,:,kl)*(Bkl(:,:,kl)'*Bkl(:,:,kl))*Dkl(:,:,kl) + DilS(:,:,ll)*(BilS(:,:,ll)'*BilS(:,:,ll))*DilS(:,:,ll),(Xkl(:,jj,kl)'*Bkl(:,:,kl)*Dkl(:,:,kl) + XilS(:,jj,ll)'*BilS(:,:,ll)*DilS(:,:,ll))')';
%             %A1(:,:,kl,ll) = (pinv(Dkl(:,:,kl)*(Bkl(:,:,kl)'*Bkl(:,:,kl))*Dkl(:,:,kl) + DilS(:,:,ll)*(BilS(:,:,ll)'*BilS(:,:,ll))*DilS(:,:,ll))*(Xkl(:,:,kl)'*Bkl(:,:,kl)*Dkl(:,:,kl) + XilS(:,:,ll)'*BilS(:,:,ll)*DilS(:,:,ll))')';
%             end
%         end
%     end
    
    %This is the parallel loop that i am working on 
    %LL = sz(4); JJ = sz(2);
    %A1 = zeros(sz(2),R);%sz(3)*sz(4),sz(4));
    
    %This appears to work in theory, but the mass spectra are super wonky
    %so I am trying other things
    %disp('Starting parallelisation')
%     for kl = 1:KL
%         DklBkltBklDkl = Dkl(:,:,kl)*(Bkl(:,:,kl)'*Bkl(:,:,kl))*Dkl(:,:,kl);
%         DklBkl = Bkl(:,:,kl)*Dkl(:,:,kl);
%         for ll = 1:LL
%             DilBiltBilDil = DilS(:,:,ll)*(BilS(:,:,ll)'*BilS(:,:,ll))*DilS(:,:,ll); %#ok
%             DilBil = BilS(:,:,ll)*DilS(:,:,ll);
%             %for jj = 1:JJ
%             %A1(jj,:,kl,ll) = fnnls(DklBkltBklDkl + DilBiltBilDil,(Xkl(:,jj,kl)'*DklBkl+Xil(:,jj,ll)'*DilBil)')'; %#ok
%             %A1(jj,:,kl,ll) = fnnls(Dkl(:,:,kl)*(Bkl(:,:,kl)'*Bkl(:,:,kl))*Dkl(:,:,kl) + DilS(:,:,ll)*(BilS(:,:,ll)'*BilS(:,:,ll))*DilS(:,:,ll),(Xkl(:,jj,kl)'*Bkl(:,:,kl)*Dkl(:,:,kl) + XilS(:,jj,ll)'*BilS(:,:,ll)*DilS(:,:,ll))')';
%             %A1(:,:,kl,ll) = (pinv(Dkl(:,:,kl)*(Bkl(:,:,kl)'*Bkl(:,:,kl))*Dkl(:,:,kl) + DilS(:,:,ll)*(BilS(:,:,ll)'*BilS(:,:,ll))*DilS(:,:,ll))*(Xkl(:,:,kl)'*Bkl(:,:,kl)*Dkl(:,:,kl) + XilS(:,:,ll)'*BilS(:,:,ll)*DilS(:,:,ll))')';
%             %end
%             A1(:,:,kl,ll) = fcnnls([],[],DklBkltBklDkl + DilBiltBilDil,(Xkl(:,:,kl)'*DklBkl+Xil(:,:,ll)'*DilBil)')';
%         end
%         %fprintf('%d\n',kl)
%     end

%First thing, trying KRB product
%DklNS = nshape(Dkl,1);
%BklNS = nshape(Bkl,2);
%DklBkl = krb(DklNS',BklNS');

% for kl = 1:KL
%     DklBkl(:,:,kl) = Bkl(:,:,kl)*Dkl(:,:,kl);
%     DklBkltBklDkl = DklBkl(:,:,kl)'*DklBkl(:,:,kl);
% end
% 
% %DklBkltBklDkl = DklBkl'*DklBkl;
% 
% %DilNS = nshape(Dil,1);
% %BilNS = nshape(Bil,2);
% %DilBil = krb(DilNS',BilNS');
% for il = 1:IL
%     DilBil(:,:,il) = Bil(:,:,il)*Dil(:,:,il);
%     DilBiltBilDil = DilBil(:,:,il)'*DilBil(:,:,il);
% end
% 
% %DilBiltBilDil = DilBil'*DilBil;
% 
% XklNS = nshape(Xkl,2);
% DklBklNS = nshape(DklBkl,2);
% XilNS = nshape(Xil,2);
% DilBilNS = nshape(DilBil,2);
% 
% A1(:,:) = fcnnls([],[],sum(DklBkltBklDkl + DilBiltBilDil,3),(XklNS*DklBklNS' + XilNS*DilBilNS')')'; %This is weird


A1 = Ail + Akl;    
%     for rr = 1:R
%         A1(:,rr) = sum(sum(A1(:,rr,:,:),3),4);
%     end
%     
%     if any(sum(A1,2) == 0) || any(isnan(sum(A1,2)))
%         A1 = 0.1*A1 + rand(size(A1,1),size(A1,2));
%         A1(A1 == 0) = 1e-6;
%     end
    A1(A1 == 0) = 1e-6;
    A1(isnan(A1)) = 1e-6;
    
    for rr = 1:R
        A1(:,rr) = A1(:,rr)./norm(A1(:,rr));
    end
    
   A = A1;
 
    
    %A = Ail + Akl;

    
    %Create Scores Matrix, F
    
    for kl = 1:sz(3)*sz(4)
        BklDkl(:,:,kl) = Bkl(:,:,kl)*Dkl(:,:,kl);
    end
    
    BklDklR = reshape(BklDkl,[sz(1),R,sz(3),sz(4)]);
    
    for il = 1:sz(1)*sz(4)
        BilDil(:,:,il) = Bil(:,:,il)*Dil(:,:,il);
    end
    
    BilDilR = permute(BilDil,[3 2 1]);
    BilDilR = reshape(BilDilR,[sz(1),R,sz(3),sz(4)]);
    
    Bklvec = 0:sz(3):sz(3)*sz(4);
    Bilvec = 0:sz(1):sz(1)*sz(4);
    
    for ll = 1:sz(4)
        %Extract one sample - Bkl
        tr2p = BklDkl(:,:,(Bklvec(ll)+1):(Bklvec(ll+1)));
        tr1p = BilDil(:,:,(Bilvec(ll)+1):(Bilvec(ll+1)));
        
        tr1p = permute(tr1p,[1,3,2]);
        tr2p = permute(tr2p,[3,1,2]);
        
        F11(:,:,:,ll) = tr1p + tr2p;
        
    end
    
    F = permute(F11,[2,3,1,4]);
    
%     for ll = 1:sz(4)
%         for kk = 1:sz(3)
%             for rr = 1:R
%                 F(:,rr,kk,ll) = BklDklR(:,rr,kk,ll) + BilDilR(:,rr,kk,ll);
%             end
%         end
%     end

%BklT = permute(Bkl,[1 3 2]);
%BklT = reshape(BklT,[sz(1),R,sz(3),sz(4)]);

%BilT = permute(Bil,[1 3 2]);
%BilT = reshape(BilT,[sz(1),R,sz(3),sz(4)]);

%F = BklT + BilT;
    
    %Normalise scores matrices
    for ll = 1:sz(4)
        for rr = 1:R
            %F(:,rr,:,ll) = F(:,rr,:,ll)./sqrt(trace(F(:,rr,:,ll)'*F(:,rr,:,ll)));
            F(:,rr,:,ll) = F(:,rr,:,ll)./norm(squeeze(F(:,rr,:,ll)),'fro');
        end
    end
    
    %Something is terribly wrong, solving for the sample-wise loadings. I
    %know this, because I have the scores, I have the mass spectra. I need
    %the loadings.
    
    %Solve for sample-wise loadings
    %Xl = nshape(Xijkl,4);
    %FR = nshape(F,2)';
    %Refold to create an RxIKJxL tensor %Every problem I ever have is with
    %reshape.
    
    %FR_IJK_L = permute(FRA,[2,1]);
    %FRA = krb(FR,A);
    
    %FR_IJK_L = reshape(FR_IJK_L,[R,sz(1)*sz(2)*sz(3),sz(4)]);
    for ll = 1:sz(4)
        %unfold the tensor, doesn't matter at this point
        X = permute(Xijkl(:,:,:,ll),[1 3 2]);
        X = reshape(X,[sz(1)*sz(3),sz(2)]);
        
        Ft = permute(F(:,:,:,ll),[1,3,2]);
        Ft = reshape(Ft,[sz(1)*sz(3),R]);
        
        %D(:,ll) = pinv(FR_IJK_L(:,:,ll)*FR_IJK_L(:,:,ll)')*FR_IJK_L(:,:,ll)*Xl(ll,:)';
        D(:,ll) = diag(pinv(Ft'*Ft)*Ft'*X*pinv(A'));
        %Dktemp = diag((pinv(Bk(:,:,kk)'*Bk(:,:,kk))*Bk(:,:,kk)'*Xkl(:,:,kk))/A'
    end
    
    %for ll = 1:sz(4)
    %    D(:,ll) = pinv(FR_IJK_L(:,:,ll)*FR_IJK_L(:,:,ll)')*FR_IJK_L(:,:,ll)*Xl(ll,:)';
    %end
    
    %Reconstruct the matrix and measure residuals
%     for kk = 1:sz(3)
%         for ll = 1:sz(4)
%             XHijkl = F(:,:,kk,ll)*diag(D(:,ll))*A';
%         end
%     end
    
    for kl = 1:sz(3)*sz(4)
       XHkl(kl) = norm(Xkl(:,:,kl) - Bkl(:,:,kl)*(Dkl(:,:,kl))*Akl','fro')^2; 
       mhkl(kl) = mkl(kl)*norm(Bkl(:,:,kl) - Pkl(:,:,kl)*Bskl,'fro')^2;
    end
    
    for il = 1:sz(3)*sz(4)
       XHil(il) = norm(Xil(:,:,il) - Bil(:,:,il)*(Dil(:,:,il))*Ail','fro')^2;
       mhil(il) = mil(il).*norm(Bil(:,:,il) - Pil(:,:,il)*Bsil,'fro')^2; 
    end
    
    mha = mA*norm(Akl - Ail,'fro')^2;
    
    %ssr2 = norm(sum(sum(sum(XHkl)))^2 + sum(sum(sum(XHil)))^2); % + norm(mA.*(Akl - Ail),'fro')^2);%/norm(Xijkl(:));
    
    %ssr2 = sqrt(sum(sum(sum(XHkl)))^2)/norm(sum(XHkl,3));% + sqrt(sum(sum(sum(XHil)))^2)/norm(sum(XHil,3));
    
    %ssr2 = norm(XHkl(:)) + norm(XHil(:));
    
    %ssr2 = norm(XHkl(:))^2 + norm(mhkl(:))^2 + norm(XHil(:))^2 + norm(mhil(:))^2 + mha;
    
    ssr2 = (sum(XHkl) + sum(mhkl))/YNorm + (sum(XHil) + sum(mhil))/YNorm + mha;
    
    %ssr2 = sum(sum(sum(sum((Xijkl - XHijkl).^2))));
    SSR(iter) = ssr2;
    
    %Recover score tensors, wrt il and kl
    %Bkl
    %Bkl = reshape(F,[sz(1),R,sz(3)*sz(4)]);
    
%     for kl = 1:sz(3)*sz(4)
%         Bkl(:,:,kl) = Bkl(:,:,kl)./norm(Bkl(:,:,kl));
%     end
%     
    %Bil
    %Bil = permute(F,[3 2 1 4]);
    %Bil = reshape(Bil,[sz(3),R,sz(1)*sz(4)]);
    
%     for il = 1:sz(1)*sz(4)
%         Bil(:,:,il) = Bil(:,:,il)./norm(Bil(:,:,il));
%     end
    
    for kk = 1:sz(3)
        for ll = 1:sz(4)
            Xh(:,:,kk,ll) = F(:,:,kk,ll)*diag(D(:,ll))*A';
        end
    end

    %pvar = 100*(1 - sum(sum(sum(sum((Xijkl - Xh).^2))))/sum(sum(sum(sum(Xijkl.^2)))));
    %SSR(iter) = sum(sum(sum(sum((Xijkl - Xh).^2))));
    
    %ssr2 = SSR(iter);
    
    if iter == 1
        varNames = {'Iteration','Absolute Error','Relative Error','SSR','mA'};
        fprintf(1,'\n%s\t\t%s\t\t%s\t\t%s\t\t%s\n',varNames{:})
        fprintf(1,' \t\t%d\t\t%e\t\t%e\t\t%e\t\t%e\n',[iter,abs(ssr2-ssr1),abs(ssr2-ssr1)/abs(ssr2),SSR(iter),mA]);
    else
        fprintf(1,' \t\t%d\t\t%e\t\t%e\t\t%e\t\t%e\n',[iter,abs(ssr2-ssr1),abs(ssr2-ssr1)/abs(ssr2),SSR(iter),mA]);
    end
%     
     %Akl = A;
%     for rr = 1:R
%         Akl(:,rr) = Akl(:,rr)./norm(Akl(:,rr));
%     end
%     
     %Ail = A;
%     for rr = 1:R
%         Ail(:,rr) = Ail(:,rr)./norm(Ail(:,rr));
    %end
    
    %Major problem here with the coupling constant
    if iter < coupit
        if iter == 1
            mA = (norm(sum(XHkl,3),'fro') + norm(sum(XHil,3),'fro'));
        end

    %mA = min(1.05*abs(sum(sum(Ail-Akl))),1e3);
    mA = min(1.02*mA,1e20);
    %mA = 10^((-mean(SNRkl)-mean(SNRil))./10)*sqrt((sum(XHil,3) + sum(XHkl,3))^2)/norm(Ail-Akl,'fro');
    else
    end
    %mA =5*10^norm(Ail-Akl,'fro');
%     %mA = mean(mil) + mean(mkl);
     %mA = 10^((-mean(SNRkl)-mean(SNRil))./10)*sqrt((sum(XHil,3) + sum(XHkl,3))^2)/norm(Ail-Akl,'fro');
%     %mA = mA*1.05;
% %     elseif iter == 2
% %         %mA = 10000;
%     elseif iter < 10
%         mA = mA*1.1;
%     else
%     end
    %end
    %mA = min(mA*1.005,1e12);
    
    if animate == 1
        
        
        %colmat = [0 0.4470 0.7410; 0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880];
        
        colmat = linspecer(R);
        
        subplot(2,1,1);
        
        for rr = 1:R
            if rr == 1
            plot(Akl(:,rr),'color',[colmat(rr,:),0.5],'LineWidth',1.5);
            else
                hold on;
                plot(Akl(:,rr),'color',[colmat(rr,:),0.5],'LineWidth',1.5);
            end
        end
        ylabel('A_{kl}','FontSize',14);
        title(sprintf('Relative Change in Error: %.4e', abs(ssr2-ssr1)/abs(ssr2)),'FontSize',14);
        hold off;
        
        subplot(2,1,2);
        for rr = 1:R
            if rr == 1
                plot(Ail(:,rr),'color',[colmat(rr,:),0.5],'LineWidth',1.5);
            else
                hold on;
                plot(Ail(:,rr),'color',[colmat(rr,:),0.5],'LineWidth',1.5);
            end
        end
        
        title(sprintf('Mass Spectral Coupling Constant: %.4e',mA),'FontSize',14);
        ylabel('A_{il}','FontSize',14);
        xlabel('m/z (Simulated)','FontSize',14);
        hold off;
        
        set(gcf,'color','w');
        
        drawnow;
        
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        
        if iter == 1
            imwrite(imind,cm,filename,'gif','Loopcount',3);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
        
    end

    iter = iter + 1;

end

%Calculate the final model, and its SSR
for kk = 1:sz(3)
    for ll = 1:sz(4)
        Xh(:,:,kk,ll) = F(:,:,kk,ll)*diag(D(:,ll))*A';
    end
end

pvar = 100*(1 - sum(sum(sum(sum((Xijkl - Xh).^2))))/sum(sum(sum(sum(Xijkl.^2)))));

disp('!Stonks!')
end

function varargout=nshape(X,f)
%NSHAPE rearrange a multi-way array
%
% [Xf,DimXf] = nshape(X,f);
%
% Refolds an N-way array so that Xf is X with index
% f as row-index, and the remaining in succesive order. For an 
% I x J x K x L four-way array this means X1 is I x JKL, X2 is
% J x IKL, X3 is K x IJL, and X4 is L x IJK
%
%
%    K  _______             
%      /      /|           1      J     2·J    J·K
%     /______/ |         1  _____________________
%    |      |  |           |      |      |      |
%    |      | /    -->     |      |      |      |        f = (Mode) 1 (same as original array)
% I  |______|/          I  |______|______|______|
%           J
%
%                          1      I     2·I    K·I
%                        1  _____________________
%                          |      |      |      |
%                  -->     |      |      |      |        f = (Mode) 2
%                        J |______|______|______|
%
%  
%                          1      I     2·I    I·J
%                        1  _____________________
%                          |      |      |      |
%                  -->     |      |      |      |        f = (Mode) 3
%                        K |______|______|______|
%
%
% f can also indicate the order (meaning the sequence) of the modes
% [Xf,DimXf] = nshape(X,[3 2 1 4]);
% will return Xf as K x JIL
%
% If the last input is not given all rearrangements are given.
% For a fourway array this would read
% [X1,X2,X3,X4]=nshape(X);
%
% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
% $ Version 1.03 $ Date 18. July 1999 $ Not compiled $
% $ Version 1.031 $ Date 18. July 1999 $ Error in help figure and now outputs new DimX $ Not compiled $
% $ Version 2.0 $ Jan 2002 $ Not compiled $ Improved speed and added permute functionality Giorgio Tomasi
ord       = ndims(X);
DimX      = size(X);
varargout = {}; %Changed this to stop matlab from screaming at me.
if nargin < 2
   f     = 0;
   do_it = ones(1,nargout);
else
   if length(f) == 1
      do_it = [1:ord] == f;
   end
end
if length(f) == 1
   for i = 1:ord
      if do_it(i)
         varargout{end+1} = reshape(permute(X,[i 1:i-1 i+1:ord]),DimX(i),prod(DimX([1:i-1 i+1:ord])));
      end
   end
   if nargin == 2
      varargout{2} = [DimX(f) DimX([1:f-1 f+1:ord])];
   end
else
   if length(f)==ord
      DimX         = DimX(f);
      varargout{1} = reshape(permute(X,f),DimX(1),prod(DimX(2:end)));
      if nargin == 2
         varargout{2} = DimX;
      end
   else
      error(['f can either be the dimension to be put first or \n a vector containing the new order of the dimensions']);
   end
end
end






