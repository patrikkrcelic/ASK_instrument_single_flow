function [VVV, STD, n1, n2] = single_flow_1_step(n1, n2, dt,ds, sizh, ...
    projtrace, dmin, dmax, tau1, tau2, tom, im, ml, siz, cutoff)

% Function performes the single flow modling for the single step. 
%
% INPUTS:
%    n1        - Thomographic volume of the J=3/2 state
%    n2        - Thomographic volume of the J=1/2 state
%    dt        - Time step
%    ds        - Horizontal step in meters
%    sizh      - Number of points in volume z component
%    projtrace - Directory of the projection luts
%    dmin      - Lower limit of the tomographic dimension index
%    dmax      - Upper limit of the tomographic dimension index
%    tau1      - "Time of life" profile of the O+(2P) J=3/2 state
%    tau2      - "Time of life" profile of the O+(2p) J=1/2 state
%    tom       - Thomographic volume of the production of O+ (2P) state
%    im        - Observed ASK2 image for the minimisation
%    ml        - Multiplication factor for the velocities, in
%    siz       - ASK image size
%                case of event with velocities larger than 4km/s
%    cutoff    - Paecentege of the maximum intensity used for mean sqare
%                error. Set to 0 to use the full image
%    
%
% OUTPUTS:
%    velovities  - Drift velocities 2D vectors, saved as Nx2 matrix, where
%                  N is the number of the timesteps
%    STD         - Standard deviations of the velocities saved in the same
%                  way as the velocities
%    n1          - New thomographic volume of the J=3/2 state
%    n2          - New thomographic volume of the J=1/2 state
%

  
  
  % degeers of freedom 
  dof = siz*siz - 2;
  % First comparision is for 0 drift 
  count=1;
  
  [n11, n22]=tomogrphy_drift(n1,n2,[0 0],dt,ds,sizh);
               
  [n11, n22]=Oxy_ion_evolution(n11,n22,tom, dt, tau2,tau1, dmin, dmax, sizh);
  
  Op2p_em=n11*1.07*10^(-1)+n22*5.63*10^(-2);

  Op2p_em(isnan(Op2p_em))=0;
  
  [A, B]=projection_fast(projtrace, 2, Op2p_em, siz, sizh);
    
  Aproj=A./B;
  
  im(isnan(im))=0;
  Aproj(isnan(Aproj))=0;
  
  im=im./max(max(im));
  Aproj=Aproj./max(max(Aproj));

  if nargin == 15
      Aproj(Aproj<cutoff)=0;
  end
  
  % Mean square error. Since we dont know the uncertainties of pixel
  % measurements this can be treated as chi-sq wit some additional
  % corrections lated.
  f(2)=sum(sum((im-Aproj).^2))/length(find(im>0));
  
  vel(count,:)=[0 0];
  er(count)=f(2);
  count=count+1;

  % Second set of comparission is with changing step size and changing lead
  % point

  VVV=[0 0];
  VV=[1 0; 0 1; 1 0; 0 1; 1 0; 0 1];
  
  dds=2;
  for ind=1:6
      if ind==3 || ind==4
          dds=1;
      elseif ind==5 || ind==6
          dds=0.5;
      end
     
    % for each dimensipn we testing 2 mean square errors 
    % around the current lead point  
    for lm=1:2
        
        if lm==1
            lm1=1;
            lam=-1*dds;
        else
            lm1=3;
            lam=dds;
        end
   
        vv=VVV+ml*lam*VV(ind,:)/norm(VV(ind,:));
      
        [n11, n22]=tomogrphy_drift(n1,n2,vv,dt,ds,sizh);
               
        [n11, n22]=Oxy_ion_evolution(n11,n22,tom, dt, tau2,tau1, dmin, dmax, sizh);
  
        Op2p_em=n11*1.07*10^(-1)+n22*5.63*10^(-2);
  
        Op2p_em(isnan(Op2p_em))=0;
  
        [A, B]=projection_fast(projtrace, 2, Op2p_em, siz, sizh);
    
        Aproj=A./B;
  
        Aproj(isnan(Aproj))=0;
        Aproj=Aproj./max(max(Aproj));
        if nargin == 13
           Aproj(Aproj<cutoff)=0;
        end
  
        f(lm1)=sum(sum((im-Aproj).^2))/length(find(im>0));
  
        vel(count,:)=vv;
        er(count)=f(lm1);
        count=count+1;
  
     end
     
     % chosing a minimum of mean square error as the new lead point 
     if f(1)==min(f)
         VVV=VVV-ml*dds*VV(ind,:);
         f(2)=f(1);
     elseif f(3)==min(f)
         VVV=VVV+ml*dds*VV(ind,:);
         f(2)=f(3);
     end
    
  end

  % Third set of comparission with array around the lead point. Needed for
  % chi-sq analysis
    
  vgrid=[-0.4 -0.4; -0.4 0; -0.4 0.4; 0 -0.4; 0 0.4; 0.4 -0.4; 0.4 0;...
      0.4 0.4; -0.2 -0.2; -0.2 0; -0.2 0.2; 0 -0.2; 0 0.2; 0.2 -0.2; 0.2 0;...
      0.2 0.2; -0.1 0; 0.1 0; 0 -0.1; 0 0.1];
  
  for ind=1:20
      vv=VVV+vgrid(ind,:);
      
      [n11, n22]=tomogrphy_drift(n1,n2,vv,dt,ds,sizh);
               
      [n11, n22]=Oxy_ion_evolution(n11,n22,tom, dt, tau2,tau1, dmin, dmax, sizh);
  
      Op2p_em=n11*1.07*10^(-1)+n22*5.63*10^(-2);
 
      Op2p_em(isnan(Op2p_em))=0;
  
      [A, B]=projection_fast(projtrace, 2, Op2p_em, siz, sizh);
    
      Aproj=A./B;
      Aproj(isnan(Aproj))=0;
      Aproj=Aproj./max(max(Aproj));
      if nargin == 13
         Aproj(Aproj<cutoff)=0;
      end
  
      vel(count,:)=vv;
      er(count)=sum(sum((im-Aproj).^2))/length(find(im>0));
      count=count+1;
      
  end
  

  % Final velocity form chi_sq
  [VVV, STD]=chi_sq_stat(er, vel, dof);
    
  [n11, n22]=tomogrphy_drift(n1,n2,VVV,dt,ds,sizh);
               
  [n11, n22]=Oxy_ion_evolution(n11,n22,tom, dt, tau2,tau1, dmin, dmax, sizh);
   
  n1=n11;
  n2=n22;
  
  Op2p_em=n11*1.07*10^(-1)+n22*5.63*10^(-2);

  Op2p_em(isnan(Op2p_em))=0;
  
  [A, B]=projection_fast(projtrace, 2, Op2p_em, siz, sizh);
  
  Aproj=A./B;
  Aproj(isnan(Aproj))=0;
  Aproj=Aproj./max(max(Aproj));
  if nargin == 13
     Aproj(Aproj<cutoff)=0;
  end
  
  er(count)=sum(sum((im-Aproj).^2))/length(find(im>0));

  vel(count,:)=VVV;
  
  if er(end)~=min(er)
      klm=find(er==min(er));
      VVV=vel(klm(1),:);
  end

end