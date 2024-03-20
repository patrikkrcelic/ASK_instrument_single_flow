
%
% The main code for the single flow model explained in Tuttle et al. (2020)
% and expanded upon in Krcelic et al. (2023).
%
% Before runing the code the runing of the IDL part if the code is 
% required. The IDL part of the code is located on lora in the directory:
% '/stp/raid1/workdir/pk1n20/single_flow/'. The matlab part can be run on 
% lora but it is very slow so I recommend transfering the IDL output 
% directory to your local PC and running the matlab part of the code localy.
%
% Single flow modeling and analysis is a complex process separated into
% distinct steps. Folow this code for the step by step analysis.
%


clear all;
close all;

%% Reading in the event setup 

% this is a thing to change
datadir='/home/pk1n20/Documents/Phd/matlab progrmi/Phd/test_event/20170201044056';


setup=importdata(strcat(datadir,'/event_setup.txt'));

event_setup.mb = setup.data(1);
event_setup.start = setup.data(2);
event_setup.warmup = setup.data(3);
event_setup.end = setup.data(4);
event_setup.step = setup.data(5);
event_setup.ds = setup.data(6);
event_setup.siz = setup.data(7);

event_setup.bck1 = 0;
event_setup.bck2 = 0;
event_setup.bck3 = 0;

% initial gues of the zenith
event_setup.zenith = [128, 128];

save(strcat(datadir,'/event_setup.mat'), 'event_setup');

emission_and_densities_read(datadir);

%% Zenith estimation trace
%  If there are auroral rays present use this block and skip the next. If
%  not skip this block and use the next one.

% change cam to 1 or 3, and find the correct image numbers and uncomment
% the variables below

% cam = 
% imnum1 = 
% imnum2 = 
% imnum3 = 

Xzen = zenith_estimator_trace(datdir, event_setup, imnum1, imnum2, imnum3, cam);

event_setup.zenith = Xzen;

save(strcat(datadir,'/event_setup.mat'), 'event_setup');


%% Zenith estimation model
%  Use this block if there are no auroral rays present in the event. This
%  part takes about 12 hours to complete

% Find the correct image numbers and uncomment the variables below

% imnum1 = 
% imnum2 = 
% imnum3 = 

Xzen = zenith_estimator_model(datdir, event_setup, imnum1, imnum2, imnum3);

event_setup.zenith = Xzen;

save(strcat(datadir,'/event_setup.mat'), 'event_setup');


%% Calculation of precipitating electron energies and fluxes for the event

event_setup.bck1 = event_background(datadir, 1);
event_setup.bck2 = event_background(datadir, 2);
event_setup.bck3 = event_background(datadir, 3);

save(strcat(datadir,'/event_setup.mat'), 'event_setup');

ENN = energy_flux_event(datadir);

% If you want to define a median filter use the function with the option
% filtnum:
%     ENN = energy_flux_event(datamb, filtnm);

EN = mean(ENN);

loss_profiles(datadir, EN);

%% Projection luts calculation for low resolution in z-direction

% this part takes about 5 minutes per camera

sizh = 41;

projection_trace(1, event_setup.ds, event_setup.siz, sizh, datadir, ...
    num2str(event_setup.mb), event_setup.zenith);

projection_trace(2, event_setup.ds, event_setup.siz, sizh, datadir, ...
    num2str(event_setup.mb), event_setup.zenith);

projection_trace(3, event_setup.ds, event_setup.siz, sizh, datadir, ...
    num2str(event_setup.mb), event_setup.zenith);

%% Tomography caluclation for low resolution in z-direction

sizh = 41;

% this part takes about 4 minutes for 4 second event
tomography_ASK1_event(datadir, sizh);

% this part takes about 1 minute for 4 second event
tomography_rest_event(datadir, 2, sizh);
tomography_rest_event(datadir, 3, sizh);

% If you want to define diferent ASK2 cutoff (default is 0.95) use the 
% function with the option cutoff:
%     tomography_rest_event(datadir, 2, sizh, cutoff);

%% Single flow model for low resolution in z-direction

sizh = 41;

[tau1, tau2, Z, projtrace] = single_flow_setup(datadir, ...
    num2str(event_setup.mb), event_setup.ds, sizh);

if event_setup.step == 2
    dt = 1/20 * 2;
elseif event_setup.step == 3
    dt = 1/32 * 3;
else
    disp('Wrong step size!');
    return
end

% dmin and dmax are tomography horizontal limits which are somethimes
% needed in case of the emission signal leakege on the edges. Change it at
% will, the default is set to take the whole tomography volume.
dmin = 1;
dmax = 150;

%[Aproj, Op2p_em, n1, n2] = ...
%    single_flow_init_cond(event_setup, datadir, dt, sizh, tau1, tau2, ...
%    dmin, dmax, projtrace);


%this part takes about 8 minutes
[velocities, STD] = single_flow_event(datadir, event_setup, projtrace, ...
    sizh, Z, dmin, dmax, tau1, tau2, dt, 1, 1, 2, 0);



%% Check and repeat

% Check if the results make sense! Sometimes there are some strange errors
% specific to the event, and each event needs to be inspected separatelly.
% Typical errors are due to the bad zenith estimation or bad energy and
% flux estimation so maybe you will need to play with filters to get it
% right.

% If everything is in order, repeat the last 3 steps with full resolution
% in z-direction

%% Projection luts calculation for high resolution in z-direction

% this part takes about 4.5 hours per camera

sizh = 2050;

projection_trace(1, event_setup.ds, event_setup.siz, sizh, datadir, ...
    num2str(event_setup.mb), event_setup.zenith);

projection_trace(2, event_setup.ds, event_setup.siz, sizh, datadir, ...
    num2str(event_setup.mb), event_setup.zenith);

projection_trace(3, event_setup.ds, event_setup.siz, sizh, datadir, ...
    num2str(event_setup.mb), event_setup.zenith);

%% Tomography caluclation for high resolution in z-direction

sizh = 2050;

% this part takes about 2.5 hours for 4 second event

tomography_ASK1_event(datadir, sizh);

% this part takes about 1.5 hours for 4 second event

tomography_rest_event(datadir, 2, sizh);


%tomography_rest_event(datadir, 3, sizh);

% If you want to define diferent ASK2 cutoff (default is 0.95) use the 
% function with the option cutoff:
%     tomography_rest_event(datadir, 2, sizh, cutoff);

%% Single flow model for high resolution in z-direction

sizh = 2050;

[tau1, tau2, Z, projtrace] = single_flow_setup(datadir, ...
    num2str(event_setup.mb), event_setup.ds, sizh);

if event_setup.step == 2
    dt = 1/20 * 2;
elseif event_setup.step == 3
    dt = 1/32 * 3;
else
    disp('Wrong step size!');
    return
end

% dmin and dmax are tomography horizontal limits which are somethimes
% needed in case of the emission signal leakege on the edges. Change it at
% will, the default is set to take the whole tomography volume.
dmin = 1;
dmax = 150;

%[Aproj, Op2p_em, n1, n2] = ...
%    single_flow_init_cond(event_setup, datadir, dt, sizh, tau1, tau2, ...
%    dmin, dmax, projtrace);


%this part takes about 7.5 hours
[velocities, STD] = single_flow_event(datadir, event_setup, projtrace, ...
    sizh, Z, dmin, dmax, tau1, tau2, dt, 1, 1, 2, 0);



















