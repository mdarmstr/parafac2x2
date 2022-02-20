function makeanimation2(Fp)

h = figure;
axis tight manual
filename = 'animation1.gif';

sz = size(Fp);

for n = 1:sz(4)
    
    for ii = 1:size(Fp,2)
        Fpn{:,:,ii} = squeeze(Fp(:,ii,:,n));
    end
    
    %temp = surfacfus(squeeze(Fp(:,1,:,n)),squeeze(Fp(:,2,:,n),squeeze(Fp(:,3,:,n))));
    temp = surfacfus(Fpn{:});
    imagesc(temp);
    set(gca,'YDir','normal');
    title('Calculated Scores, PARAFAC2\times2','FontSize',14);
    xlabel('Modulations, ^1t_R','FontSize',14);
    ylabel('Acquisitions, ^2t_R','FontSize',14);
    set(gcf,'color','w');
    
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    if n ==1
        imwrite(imind,cm,filename,'gif','Loopcount',3);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end
    
    