function [Aproj, Op2p_em, n1, n2] = ...
    single_flow_init_cond(event_setup, datdir, dt, sizh, tau1, tau2, ...
    dmin, dmax, projtrace)

% Function calculates the initial conditions for the singel folow model 
%
% INPUTS:
%    event_setup - Structure with metadata of the event
%    datadir     - Data directory of the event
%    dt          - Time step
%    sizh        - Number of points in volume z component
%    Z           - Decontamination matrix
%    tau1        - "Time of life" profile of the O+(2P) J=3/2 state
%    tau2        - "Time of life" profile of the O+(2p) J=1/2 state
%    dmin        - Lower limit of the tomographic dimension index
%    dmax        - Upper limit of the tomographic dimension index
%    projtrace   - Directory of the projection luts
%    
%
% OUTPUTS:
%    Aproj   - Modeled image of the O+(2P) state
%    Op2p_em - Thomographic 3D volume of the O+(2P) emission
%    n1 - Thomographic 3D density volume of the J=3/2 state
%    n2 - Thomographic 3D density volume of the J=1/2 state
%

  
 n1 = zeros(150,150,sizh);
 n2 = zeros(150,150,sizh);

 for pk=event_setup.start:event_setup.step:event_setup.warmup

    sp = num2str(pk,'%05.f');
  
    datfile=strcat(datdir,'/results/tomography/ASK2/tomography_ask2_low_', sp,'.mat');
    load(datfile, 'tom');
  
    [n11, n22]=tomogrphy_drift(n1,n2,[0 0],dt,event_setup.ds,sizh);
  
    [n1, n2]=Oxy_ion_evolution(n11,n22,tom, dt, tau2,tau1, dmin, dmax, sizh);
  
    Op2p_em=n1*1.07*10^(-1)+n2*5.63*10^(-2);
 
  end
  
  
  Op2p_em(isnan(Op2p_em))=0;   
  
  [A, B]=projection_fast(projtrace, 2, Op2p_em, event_setup.siz, sizh);

  Aproj=A./B;
  
  Aproj=Aproj./max(max(Aproj));

end
