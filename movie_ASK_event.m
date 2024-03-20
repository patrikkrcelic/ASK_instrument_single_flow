function movie_ASK_event(datdir, event_setup)

% Function creates the movie of the ASK event.
%
% INPUTS:
%   datdir      - daata directoty
%   event_setup - event metadata
%

num = 1;

for p=event_setup.start : event_setup.step : event_setup.end
    
    sp=num2str(p,'%05.f');

    f=figure('visible','off');

    ASK = importdata(strcat(datdir,'/ASK1/ask1_',sp,'.txt'));
    ASK = image_trimming(ASK, event_setup.bck1);

    subplot(1,3,1)
    h=pcolor(ASK); 
    axis equal;
    set(h, 'EdgeColor', 'none');
    xlim([1,256])
    ylim([1,256])
    set(gca,'xtick', []);
    set(gca,'ytick', []);
    title('ASK1');

    ASK=importdata(strcat(datdir,'/ASK2/ask2_',sp,'.txt'));
    ASK= image_trimming(ASK, event_setup.bck2);

    subplot(1,3,2)
    h=pcolor(ASK); 
    axis equal;
    set(h, 'EdgeColor', 'none');
    xlim([1,256])
    ylim([1,256])
    set(gca,'xtick', []);
    set(gca,'ytick', []);
    title('ASK2');

    ASK=importdata(strcat(datdir,'/ASK3/ask3_',sp,'.txt'));
    ASK= image_trimming(ASK, event_setup.bck3);

    subplot(1,3,3)
    h=pcolor(ASK); 
    axis equal;
    set(h, 'EdgeColor', 'none');
    xlim([1,256])
    ylim([1,256])
    set(gca,'xtick', []);
    set(gca,'ytick', []);
    title('ASK2');


    M=getframe(f);
 
    if num==1
       [im,map] = rgb2ind(M.cdata,256,'nodither');
       im(1,1,1,20) = 0;
       num = num + 1;
      
    else
        
       im(:,:,1,num-1) = rgb2ind(M.cdata,map,'nodither');
       num = num + 1;
     
    end

end

imwrite(im,map, strcat(datdir,'/event.gif'),'DelayTime',0.1,'LoopCount',20); 


end