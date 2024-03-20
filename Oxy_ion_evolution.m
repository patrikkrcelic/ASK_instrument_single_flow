function [n11 n22]=Oxy_ion_evolution(n11,n22,tom, dt, tau2,tau1, dmin, dmax, sizh)

% Function solves O+(2P) continuity equation.
%
% INPUTS:
%    tom  - O+(2P) production 3D volume
%    dt   - Time resolution
%    tau2 - 'time of life' profile of O+(2P) 1/2 upper state
%    tau2 - 'time of life' profile of O+(2P) 3/2 upper state
%    dmin - Volume minimum horizontal direction index
%    dmax - Volume maximum horizontal direction index
%    sizh - Number of points in volume z component
%
% INPUTS/OUTPUTS:
%    n11  - Initial/Output O+(2P) 1/2 upper state density
%    n22  - Initial/Output O+(2P) 3/2 upper state density
%

  for k=1:sizh
    for i=dmin:dmax
       for j=dmin:dmax
           n11(i,j,k)=n11(i,j,k)+0.7315*tom(i,j,k)*dt;
           if n11(i,j,k)<0
              n11(i,j,k)=0;
           end
           n22(i,j,k)=n22(i,j,k)+0.2685*tom(i,j,k)*dt;
           if n22(i,j,k)<0
               n22(i,j,k)=0;
           end
        end
    end
    n11(:,:,k)=n11(:,:,k)*exp(-(dt)/(tau2(k)));
    n22(:,:,k)=n22(:,:,k)*exp(-(dt)/(tau1(k)));
  end
end