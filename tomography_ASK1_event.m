function [] = tomography_ASK1_event(datamb, sizh)

% Function calculates and fills in tomography volume for ASK1 and makes 
% energy and flux maps in tomography volume reference frame for the whole 
% event. Tomography volumes are saved in '/results/tomography/' directory
% and energy and flux maps are saved in '/results/energy_and_flux/' 
% directory.
%
% INPUTS:
%    datamb    - Data directory of the event
%    sizh      - Number of points in volume z component
%

  load(strcat(datamb,'/emissions.mat'));
  load(strcat(datamb,'/densities.mat'));
  load(strcat(datamb,'/event_setup.mat'));

  datdir=strcat(datamb,'/results/energy_and_flux/');
  resdir=strcat(datamb,'/results/tomography/ASK1/');

  A1=importdata(strcat(datamb, '/', num2str(event_setup.mb), '_el_az_ask1_lut.txt'));
  El=A1(:,1:256);
  Az=A1(:,257:end);
      
  em = N2_em;

  for p=event_setup.start : event_setup.step : event_setup.end
    
    sp=num2str(p,'%05.f');
     
    datfile=strcat(datamb,'/results/energy_and_flux/','image_',sp,'.mat');
    load(datfile);

    [tom, Men, Mflux] = tomography_ASK1(Aen, Aflux, El, Az, em, alt, en, ...
        event_setup.siz, sizh, event_setup.zenith, event_setup.ds);

    resfile=strcat(datdir,'model_',sp,'.mat');
    save(resfile, 'Men', 'Mflux');
    resfile=strcat(resdir,'tomography_ask1_low_', sp,'.mat');
    save(resfile, 'tom');

  end

end









