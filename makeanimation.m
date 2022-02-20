function makeanimation(plts)

h = figure;
axis tight manual
filename = 'animation.gif';

sz = size(plts);

for n = 1:sz(4)
    
    imagesc(plts(:,:,:,n));
    set(gca,'YDir','normal');
    title('Simulated Scores, 500 SNR','FontSize',14);
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
    
    