function EN = energy_flux_event(datamb, filtnm)

% Function calculates and creates electron precipitation energy and flux 
% maps from each ASK image in the event. The maps are saved in 
% '/results/energy_and_flux/' directory of the event specified with datamb.
%
% INPUTS:
%    datamb - Data directory of the event
%    filtnm - Optional, filter number for median filter of the images

 data1=strcat(datamb,'/ASK1/');
 data3=strcat(datamb,'/ASK3/');
 resdir=strcat(datamb,'/results/energy_and_flux/');
 
 X1=importdata(strcat(datamb, '/ratio_energy.txt'));
 X2=importdata(strcat(datamb, '/energy_fluxint.txt'));

 X1=X1(2:end,:);

 load(strcat(datamb,'/emissions.mat'));
 load(strcat(datamb,'/event_setup.mat'));
 zenith=event_setup.zenith;

 EN = zeros(round((event_setup.end-event_setup.start)/event_setup.step),1);

 num = 1;

 for p=event_setup.start : event_setup.step : event_setup.end
    
    sp=num2str(p,'%05.f');
 
    datfile=strcat(data1,'ask1_',sp,'.txt');
    A1=importdata(datfile);
    datfile=strcat(data3, 'ask3_',sp,'.txt');
    A3=importdata(datfile);
    
    A1 = image_trimming(A1, event_setup.bck1);
    A3 = image_trimming(A3, event_setup.bck3);

    if nargin == 2
       [Aen, Aflux, ENN] = EnergyFlux (A1, A3, alt, en, N2_em, OI_em, zenith, X1, X2, filtnm);
    else
       [Aen, Aflux, ENN] = EnergyFlux (A1, A3, alt, en, N2_em, OI_em, zenith, X1, X2);
    end

    resfile=strcat(resdir,'image_',sp,'.mat');

    save(resfile, 'Aen', 'Aflux');

    EN(num) = ENN;
    num = num + 1;
  
 end
end

    