function movie_ASK2_model(datdir, event_setup, med_average, decont)

% Function creates the movie of the ASK2 model comparisson.
%
% INPUTS:
%   datdir      - daata directoty
%   event_setup - event metadata
%    med_average - Size of median filter for spatial averaging
%    decont      - Optional, decontamination ratio. Dafault is 0.03
%

num = 1;

Z=importdata(strcat(datdir, '/', num2str(event_setup.mb), '_decontamination_lut.txt'));

tstart = event_setup.warmup + event_setup.step;
tstop = event_setup.end - event_setup.step;

for p = tstart : event_setup.step : tstop
    
    sp=num2str(p,'%05.f');

    f=figure('visible','off');

    ASK1 = importdata(strcat(datdir,'/ASK1/ask1_',sp,'.txt'));
    ASK1 = image_trimming(ASK1, event_setup.bck1);

    ASK2 = importdata(strcat(datdir,'/ASK2/ask2_',sp,'.txt'));
    
    % time averaging

    if event_setup.step == 2
        sp1 = num2str(p+1,'%05.f');
        A2=importdata(strcat(datdir, '/ASK2/ask2_',sp1,'.txt'));
        ASK2 = ASK2 + A2;

        ASK2 = ASK2/2;

    else
        sp1 = num2str(p+1,'%05.f');
        sp2 = num2str(p+2,'%05.f');
        A2=importdata(strcat(datdir, '/ASK2/ask2_',sp1,'.txt'));
        ASK2 = ASK2 + A2;
    
        A2=importdata(strcat(datdir, '/ASK2/ask2_',sp2,'.txt'));
        ASK2 = ASK2 + A2;
    
        ASK2 = ASK2/3;
    end

    ASK2 = image_trimming(ASK2, event_setup.bck2);

    % decontamination

    A1 = zeros(256, 256);

    for i=1:event_setup.siz
      for j=1:event_setup.siz
          Zx=round(Z(i,j+event_setup.siz))+1;
          Zy=round(Z(i,j))+1;
          if Zx > event_setup.siz || Zx<1
              continue
          end

          if Zy > event_setup.siz || Zy<1
              continue
          end

          A1(i,j) = ASK1(Zx,Zy); 
      end
    end

    if nargin == 4
      ASK2 = ASK2 - decont*A1;
    else
      ASK2 = ASK2 - 0.03*A1;
    end

    if med_average > 1
       ASK2 = median_filter(ASK2, med_average);
    end
    
 

    subplot(1,2,1)
    h=pcolor(ASK2/max(max(ASK2))); 
    axis equal;
    set(h, 'EdgeColor', 'none');
    xlim([1,256])
    ylim([1,256])
    set(gca,'xtick', []);
    set(gca,'ytick', []);
    title('ASK2 Observed');



    datfile = strcat(datdir, '/results/modelled_emissions/ASK2/image_ask2_',sp,'.mat');
    load(datfile, 'Aproj');

    subplot(1,2,2)
    h=pcolor(Aproj/max(max(Aproj))); 
    axis equal;
    set(h, 'EdgeColor', 'none');
    xlim([1,256])
    ylim([1,256])
    set(gca,'xtick', []);
    set(gca,'ytick', []);
    title('ASK2 Modeled');



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

imwrite(im,map, strcat(datdir,'/model.gif'),'DelayTime',0.1,'LoopCount',20); 


end