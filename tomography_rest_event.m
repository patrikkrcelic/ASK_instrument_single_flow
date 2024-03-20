function [] = tomography_rest_event(datamb, cameranum, sizh, cutoff)

% Function calculates and fills in tomography volume for ASK2 or ASK3 for 
% the whole event. Tomography volumes are saved in '/results/tomography/' 
% directory. Function tomography_ASK1_event needs to be run before runing 
% this function. 
%
% INPUTS:
%    datamb    - Data directory of the event
%    cameranum - ASK camera number, either 2 or 3
%    sizh      - Number of points in volume z component
%    cutoff    - Optional, flux cutoff for ASK2 emission in range between
%                0 and 1. Default is set to 0.95


  load(strcat(datamb,'/emissions.mat'), 'prod', 'OI_em', 'OI_D_em',...
      'alt','en');
  load(strcat(datamb,'/event_setup.mat'), 'event_setup');

  datdir=strcat(datamb,'/results/energy_and_flux/');
  resdir=strcat(datamb,'/results/tomography/ASK',num2str(cameranum),'/');
 
  if cameranum == 2
      em = prod;
  elseif cameranum == 3
      em = OI_em + OI_D_em;
  else
      disp('Wrong camera number input!');
      return
  end  

  for p=event_setup.start : event_setup.step : event_setup.end
    
    sp=num2str(p,'%05.f');
    
    datfile=strcat(datdir,'/model_',sp,'.mat');
    resfile=strcat(resdir,'tomography_ask',num2str(cameranum),'_low_', sp,'.mat');
    load(datfile, 'Men', 'Mflux');
    datfile=strcat(datamb,'/results/energy_and_flux/','image_',sp,'.mat');
    load(datfile, 'Aflux');

    if nargin == 4
        tom=tomography_rest(em, alt, en, Men, Mflux, Aflux, cameranum, sizh, cutoff);
    else
        tom=tomography_rest(em, alt, en, Men, Mflux, Aflux, cameranum, sizh);
    end
   
    save(resfile, 'tom');
  end
   
end





              