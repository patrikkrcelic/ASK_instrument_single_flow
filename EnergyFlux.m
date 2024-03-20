function [Aen, Aflux, EN] = EnergyFlux (A1, A3, alt, en, N2_em, OI_em, zenith, X1, X2, filtnm)

% Function calculates electron precipitation energy and flux maps from ASK images.
%
% INPUTS:
%    A1     - ASK1 image
%    A3     - ASK3 image
%    alt    - ion chemistry altitudes
%    en     - ion chemistry mono energies
%    N2_em  - N2 emissions from ion chemistry mono run
%    OI_em  - OI + OI_d emissions from ion chemistry mono run
%    zenith - magnetic zenith position
%    X1     - ratio and energy interpolation values
%    X2     - energy and fluxes interpolation values
%    filtnm - optional, filter number for median filter of the images
%
% OUTPUTS:
%    Aen    - energy map same size as ASK1 and ASK3 images
%    Aflux  - flux map same size as ASK1 and ASK3 images
%


    s = size(A1);
    AA=reshape(A1,[s(1)*s(2),1]);
    A1m = prctile(AA, 90);  % it is value of the arc, not the background
    AA=reshape(A3,[s(1)*s(2),1]);
    A3m = prctile(AA, 90);

    EN=ratio2energy(A3m/A1m,X1);
    J=find(en>=EN);
    h3=alt(J(1),max(OI_em(J(1),:))==OI_em(J(1),:));
    h1=alt(J(1),max(N2_em(J(1),:))==N2_em(J(1),:));
    
    A1(A1<=0)=1;
    A3(A3<=0)=1000;

    if nargin == 10
       A3=median_filter(A3, filtnm);
       A1=median_filter(A1, filtnm);
    end

    A = zeros(s(1),s(2));
    Afluxcal = zeros(s(1),s(2));

    for j=1:256
        for k=1:256
            du=round((1-h1/h3)*(zenith(1)-k));
            dv=round((1-h1/h3)*(zenith(2)-j));
            
           if h1<h3
              A(j,k)=A3(j+dv,k+du)/A1(j,k);
           else
              A(j,k)=A3(j,k)/A1(j-dv,k-du);
           end
            
        end
    end
    
    for i=1:s(1)
        for j=1:s(2)
            A(i,j)=ratio2energy(A(i,j),X1);
            Afluxcal(i,j)=ratio2energy(A(i,j),X2);        
        end
    end
  
  Aen=A;
  Aflux=A1./Afluxcal;
 
end
