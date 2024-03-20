function tom = tomography_rest(em, alt, en,  Men, Mflux, Aflux, cameranum, sizh, cutoff)

% Function calculates and fills in tomography volume for ASK2 or ASK3 using
% precipitating energy and flux from tomography volume reference frame
% calculated in tomography_ASK1.
% 
%
% INPUTS:
%    em     - volume emission rates from ion chemistry mono run
%    alt    - ion chemistry altitudes
%    en     - ion chemistry mono energies
%    Men    - energy map in tomography reference frame
%    Mflow  - flux map in tomography reference frame
%    Aflow  - flux map same size as ASK1 and ASK3 images
%    cameranum - ASK camera number, either 2 or 3
%    sizh      - Number of points in volume z component
%    cutoff    - Optional, flux cutoff for ASK2 emission in range between
%                0 and 1.


    tom=ones(150,150,sizh)*0.1;
    s = size(Aflux);
    AA=reshape(Aflux,[s(1)*s(2),1]);
    AA=sort(AA);

    if nargin == 9
       Afmax=AA(round(256*256*cutoff));
    else
        Afmax=AA(round(256*256*0.95));
    end

    for i=1:150
        for j=1:150
            ijk=find(Men(i,j,:)~=0);
            aur=zeros(1,201);
            
            if length(ijk)<1
                continue
            end
            lmn=ijk(end);
            for kl=1:ijk(end)
                b=Men(i,j,kl);
                if isnan(b)
                   continue
                end
                if en(end)<b
                   I=length(en);
                else
                   I=find(en>=b);
                end
                a=alt(I(1),:);
                if cameranum==3
                    aur=aur+(em(I(1),:))*Mflux(i,j,kl);
                elseif cameranum==2
                    Mf=Mflux(i,j,kl);
                  if Mf<Afmax
                      lmn=lmn-1;
                    continue
                  end
                    aur=aur+em(I(1),:)*Mflux(i,j,kl);
                else
                    disp('Something went wrong!');
                    return
                end
            end

            if lmn<1
                continue
            end
            aur=aur/lmn;
        
           for ij=1:length(a)/2
               aa=a(ij);
               a(ij)=a(end-ij+1);
               a(end-ij+1)=aa;
               aa=aur(ij);
               aur(ij)=aur(end-ij+1);
               aur(end-ij+1)=aa;
           end
           
           h=80+(0:sizh-1)*0.2*2050/sizh;
           z = linear_interpolation_extrapolation(aur, a, h);

           tom(i,j,:)=z;  

        end
    end
   
end

            