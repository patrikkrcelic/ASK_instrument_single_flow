function [tau1, tau2, Z, projtrace] = single_flow_setup(datdir, mb, ds,sizh)

% Function calculates variable needed for the fingle flow model 
%
% INPUTS:
%    datadir - Data directory of the event
%    mb      - Megablock of the event
%    ds      - Horizontal resolution in meters
%    sizh    - Number of points in volume z component
%    
%
% OUTPUTS:
%    tau1      - "Time of life" profile of the O+(2P) J=3/2 state
%    tau2      - "Time of life" profile of the O+(2p) J=1/2 state
%    Z         - Decontamination matrix
%    projtrace - Directory of the projection luts
%



  load(strcat(datdir,'/emissions.mat'), 'alt');
  load(strcat(datdir,'/loss_profiles.mat'), 'tau_1_2', 'tau_3_2');

  Z=importdata(strcat(datdir,'/', mb ,'_decontamination_lut.txt'));
  projtrace=strcat(datdir, '/results/projection_luts');

  a=alt(1,:);

  h=80+(0:sizh-1)*ds*2050/sizh;

  tau_1_2 = change_order(tau_1_2);
  tau_3_2 = change_order(tau_3_2);
  a = change_order(a);

  tau1 = linear_interpolation_extrapolation(tau_1_2, a, h);
  tau2 = linear_interpolation_extrapolation(tau_3_2, a, h);

end
