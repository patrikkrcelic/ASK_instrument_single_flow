function [tom, Men, Mflux] = tomography_ASK1(Aen, Aflux, El, Az, em, alt, en, siz, sizh, zenith, ds)

% Function does tomography for ASK1 emission from the ASK energy and flux 
% images and ion chemistry emision rates profile. Additionaly, function 
% makes energy and flux maps in tomography reference plane needed to fill 
% in tomography volume of ASK2 and ASK3 emissions.
%
% INPUTS:
%    Aen    - ASK energy image
%    Aflux  - ASK flux image
%    El     - Matrix of pixel elevation values
%    Az     - Matrix of pixel azimuth values
%    em     - N2 emissions from ion chemistry mono run
%    alt    - ion chemistry altitudes
%    en     - ion chemistry mono energies
%    siz    - size of the image
%    sizh   - number of points in volume z component
%    zenith - magnetic zenith position
%    ds     - horizontal resolution of 3D volume expresed in km
%
% OUTPUTS:
%    tom    - tomographic volume with size 150x150xsizh
%    Men    - energy map in tomography reference frame
%    Mflux  - flux map in tomography reference frame
%

   X(1)=round(zenith(2));
   X(2)=round(zenith(1));
   [rot, x0, y0, ~] = ASK_rotation_matrix(El, Az, X, siz);
   
   tom=zeros(150,150,sizh);
   Men=zeros(150,150,50);
   Mflux=zeros(150,150,50);
   m=zeros(150,150);

   ds = ds*1000;

   for i=1:siz
       for j=1:siz
           b=Aen(i,j);
           if isnan(b)
               continue
           end

           I=find(en>=b);
           
           if en(end)<b
               I=length(en);
           end
           aur=em(I(1),:);
           a=alt(I(1),:);
           IJ=find(aur==max(aur));
           h=a(IJ(1));
        
           for ij=1:length(a)/2
               aa=a(ij);
               a(ij)=a(end-ij+1);
               a(end-ij+1)=aa;
               aa=aur(ij);
               aur(ij)=aur(end-ij+1);
               aur(end-ij+1)=aa;
           end
           
           el=pi/2-El(i,j);
           az=Az(i,j);

           x=h*1000*tan(el)*sin(az);
           y=h*1000*tan(el)*cos(az);


           Rt=[x y h*1000/cos(el)]*rot;

           x=Rt(1)-x0;
           y=Rt(2)-y0;

           xp=floor(x/ds)+76;
           yp=floor(y/ds)+76;

           if xp>150 || xp<1 || yp>150 || yp<1
              continue
           end

           h=80+(0:sizh-1)*ds/1000*2050/sizh;
           z = linear_interpolation_extrapolation(aur, a, h);
           
           z1 = zeros(1, sizh);

           for k=1:sizh
               z1(k) =tom(xp,yp,k);
           end

           tom(xp,yp,:)=z1+z*Aflux(i,j);
           ijk=find(Men(xp,yp,:)==0);
           if length(ijk)<1
               continue
           end
          
           ijk=ijk(1);
           Men(xp,yp,ijk)=en(I(1));
           Mflux(xp,yp,ijk)=Aflux(i,j);
           m(xp,yp)=m(xp,yp)+1;

       end   
   end

   for i=1:sizh
       tom(:,:,i)=tom(:,:,i)./m;
   end

   tom(isnan(tom))=0;
   tom(tom==0)=0.1;

end