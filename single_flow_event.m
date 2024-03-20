function [velocities, STD] = single_flow_event(datadir, event_setup,...
    projtrace, sizh, Z, dmin, dmax, tau1, tau2, dt, ml, temp_avr,...
    med_average, cutoff, decont) 


% Function performes the single flow modling for the entire event. 
%
% INPUTS:
%    datadir     - Data directory of the event
%    event_setup - Structure with metadata of the event
%    projtrace   - Directory of the projection luts
%    sizh        - Number of points in volume z component
%    Z           - Decontamination matrix
%    dmin        - Lower limit of the tomographic dimension index
%    dmax        - Upper limit of the tomographic dimension index
%    tau1        - "Time of life" profile of the O+(2P) J=3/2 state
%    tau2        - "Time of life" profile of the O+(2p) J=1/2 state
%    dt          - Time step
%    ml          - Multiplication factor for the velocities, in
%                  case of event with velocities larger than 4km/s
%    temp_avr    - Number of ASK images for time averaging 
%    med_average - Size of median filter for spatial averaging
%    cutoff      - Paecentege of the maximum intensity used for mean sqare
%                  error. Set to 0 to use the full image
%    decont      - Optional, decontamination ratio. Dafault is 0.03
%    
%
% OUTPUTS:
%    velovities  - Drift velocities 2D vectors, saved as Nx2 matrix, where
%                  N is the number of the timesteps
%    STD         - Standard deviations of the velocities saved in the same
%                  way as the velocities
%



% initial conditions

[~, ~, n1, n2] = ...
   single_flow_init_cond(event_setup, datadir, dt, sizh, tau1, tau2, ...
   dmin, dmax, projtrace);

N1 = n1;
N2 = n2;

tstart = event_setup.warmup + event_setup.step;
tstop = event_setup.end - event_setup.step;

ijk=1;
pcheck=10;
tmax=tstop-tstart;
velocities=zeros(round(tmax/event_setup.step),2);
STD=zeros(round(tmax/event_setup.step),2);


% The loop for images in the event
for pk=tstart:event_setup.step:tstop
  
  sp=num2str(pk,'%05.f');

  imdir=strcat(datadir,'/ASK2/');
  imdir1=strcat(datadir,'/ASK1/');
  imfile=strcat(imdir,'ask2_', sp,'.txt');
  im=importdata(imfile);

  imfile=strcat(imdir1,'ask1_', sp,'.txt');
  im2=importdata(imfile);

  im2 = image_trimming(im2, 0);

  im3=zeros(event_setup.siz, event_setup.siz);
  
  % creating a image for ASK1 noise substraction
  
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

          im3(i,j)=im2(Zx,Zy); 
      end
  end
  
  
  if temp_avr == 0
      im13=im;
  else
      if event_setup.step == 2
          sp1 =  num2str(pk+1,'%05.f');
          im11 = importdata(strcat(imdir,'ask2_',sp1,'.txt'));
          im13 = (im+im11)/2;
      elseif event_setup.step == 3
          sp1 =  num2str(pk+1,'%05.f');
          sp2 =  num2str(pk+2,'%05.f');
          im11 = importdata(strcat(imdir,'ask2_',sp1,'.txt'));
          im12 = importdata(strcat(imdir,'ask2_',sp2,'.txt'));
          im13 = (im+im11+im12)/3;
      end
  end

  im13 = image_trimming(im13, 0);
 
  if nargin == 15
      im13=im13-decont*im3;
  else
      im13=im13-0.03*im3;
  end
  
  if med_average > 1
     im13=median_filter(im13,med_average);
  end

  %im(1,:)=im(2,:);
  %im(256,:)=im(254,:);
  %im(255,:)=im(254,:);
  %im(:,1)=im(:,2);
  %im(:,256)=im(:,254);
  %im(:,255)=im(:,254);

  
  datfile=strcat(datadir,'/results/tomography/ASK2/tomography_ask2_low_', sp,'.mat');
  load(datfile, 'tom');
  

  % modell image minimisation
  if cutoff > 0
     [VVV, std, n1, n2] = single_flow_1_step(n1, n2, dt,event_setup.ds,...
         sizh, projtrace, dmin, dmax, tau1, tau2, tom, im13, ml,...
         event_setup.siz, cutoff);
  else
     [VVV, std, n1, n2] = single_flow_1_step(n1, n2, dt, event_setup.ds,...
         sizh, projtrace, dmin, dmax, tau1, tau2, tom, im13, ml, event_setup.siz);
  end

  velocities(ijk,:)=VVV;
  STD(ijk,:)=std;
  
  ijk=ijk+1;
  
  % display the percentage of completion
  if (pk-tstart)/tmax*100>=pcheck
      display(strcat(num2str(pcheck),'%'));
      pcheck=pcheck+10;
  end

end

vel = velocities;

% rotation to camera field of view due to the IDL and Matlab differences 

Rot1=[cos(-pi/2) -sin(-pi/2); sin(-pi/2) cos(-pi/2)];

n = size(velocities); 

for i = 1 : n(1)
    vv = velocities(i,:);
    vv = vv*Rot1;
    velocities(i,:) = vv;

    vv = STD(i,:);
    vv = vv*Rot1;
    STD(i,:) = vv;
end
  
resdir2=strcat(datadir,'/results/');
resfile2=strcat(resdir2, 'velocities.mat');
save(resfile2, 'velocities', 'STD');

% smoothing of the velocities and saving images for the movie

ijk=1;
vel(:,1)=filter([1 1 1 1 1]/5,1,vel(:,1));
vel(:,2)=filter([1 1 1 1 1]/5,1,vel(:,2));

n1=N1;
n2=N2;
resdir=strcat(datadir,'/results/modelled_emissions/ASK2/');

 for pk=tstart:event_setup.step:tstop
  
    sp=num2str(pk,'%05.f'); 
  
    datfile=strcat(datadir,'/results/tomography/ASK2/tomography_ask2_low_', sp,'.mat');
    load(datfile, 'tom');
  
    VVV=vel(ijk,:);
 
    [n11, n22]=tomogrphy_drift(n1,n2,VVV,dt,event_setup.ds,sizh);
               
    [n1, n2]=Oxy_ion_evolution(n11,n22,tom, dt, tau2,tau1, dmin, dmax, sizh);
  
    Op2p_em=n1*1.07*10^(-1)+n2*5.63*10^(-2); % Einstein coefficients
    Op2p_em(isnan(Op2p_em))=0;
     
    [A, B]=projection_fast(projtrace, 2, Op2p_em, event_setup.siz, sizh);

    Aproj=A./B;
    
    resfile = strcat(resdir, 'image_ask2_',sp,'.mat');
    save(resfile, 'Aproj');
    ijk=ijk+1;

 end

end