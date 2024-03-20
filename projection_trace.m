function [] = projection_trace(camera, ds, siz, sizh, datdir, mb, zenith)

% Function creates lut-files for faster ASK projection. Projection is
% calculated of the camera defined in the input (1,2 or 3). Data is saved
% in directory datdir + '/results/projection_lut' and luts are saced as 
% 'camera_N_XY_HH.mat', where N is the number of camera and HH is the
% height label in 3D volume. 
%
% INPUTS:
%    camera - ASK camera for projection
%    ds     - horizontal resolution of 3D volume expressed in km
%    siz    - size of the image
%    sizh   - number of points in volume z component
%    datdir - data directory of the event
%    mb     - megablock of the event
%    zenith - magnetic zenith position    
%
%


    X(1)=round(zenith(2));
    X(2)=round(zenith(1));

    ds = ds*1000; % transform from km to m

    A=importdata(strcat(datdir, '/', mb, '_el_az_ask1_lut.txt'));
    Az=A(:,257:end);
    El=A(:,1:256);

    [rot, x0, y0, z0] = ASK_rotation_matrix(El, Az, X, siz);

     
    A=importdata(strcat(datdir, '/', mb, '_el_az_ask', num2str(camera), '_lut.txt'));
    Az=A(:,257:end);
    El=A(:,1:256);

    odir=strcat(datdir, '/results/projection_luts');


    hh=(80+(0:sizh-1)*0.2*2050/sizh+0.1*2050/sizh)*1000;
    
    for k=1:sizh
        h=hh(k);
        imat=zeros(150,150);
        jmat=zeros(150,150);

        rad=zeros(150,150);
       
       for i=1:150
           for j=1:150
            
               x=(i-1)*ds+ds/2-150*ds/2;
               y=(j-1)*ds+ds/2-150*ds/2;
               z=h-z0;

               x=x+x0;
               y=y+y0;
               z=z+z0;

               Xf = [x y z]/rot;
               x=Xf(1);
               y=Xf(2);
               z=Xf(3);
               az=atan(x/y)+pi;
               r=sqrt(x^2+y^2+z^2);
               el=acos(z/r);
               el=pi/2-el;

               rad(i,j)=r;

               Xpos=ea2xy(El,Az,el,az);
               
            
               if Xpos(1)>siz || Xpos(1)<1 || Xpos(2)>siz || Xpos(2)<1
                 continue
               end

               imat(i,j)=Xpos(1);
               jmat(i,j)=Xpos(2);
           end
       end

       save(strcat(odir, '/camera_',num2str(camera), '_XY_', num2str(k), '.mat'), 'imat', 'jmat', 'rad');
    end

end