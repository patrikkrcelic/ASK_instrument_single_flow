function background = event_background(datamb, cameranum)

% Function calculates the background emission in the event.
%
% INPUTS:
%    datamb    - Data directory of the event
%    cameranum - Number of the ask camera
%
% OUTPUT:
%    background - Background value of the event
%

 datadir=strcat(datamb,'/ASK', num2str(cameranum), '/');
 load(strcat(datamb,'/event_setup.mat'));

 background = 10000;


 for p=event_setup.start : event_setup.step : event_setup.end
    
    sp=num2str(p,'%05.f');

    datfile=strcat(datadir,'ask',num2str(cameranum), '_',sp,'.txt');
    A=importdata(datfile);

    s = size(A);
    AA=reshape(A,[s(1)*s(2),1]);
    Amin = prctile(AA,15);

    if Amin < background
        background = Amin;
    end
 end

 if background < 0
     background = 0;
 end

end





