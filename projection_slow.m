function [A, B] = projection_slow(tom, El, Az, siz, sizh, zenith, ds)

% Function does slow projection from tomography volume to ASK image,
% withouth projection 'look up tables'.
%
% INPUTS:
%    tom    - Volume emission rates 3D volume
%    El     - Matrix of pixel elevation values
%    Az     - Matrix of pixel azimuth values
%    sizh   - Number of points in volume z component
%    zenith - magnetic zenith position
%    ds     - Horizontal spatial resolution of the volume in km
%
% OUTPUTS:
%    A  - Uncalibrated projected image
%    B  - Calibration image
%
 
    
    A=zeros(siz,siz);
    Ap=zeros(siz,siz);
    B=zeros(siz,siz);
    Bp=zeros(siz,siz);

    X(1)=round(zenith(2));
    X(2)=round(zenith(1));
    
    ds=ds*1000;
    
    [rot, x0, y0, z0] = ASK_rotation_matrix(El, Az, X, siz);


    s=0;
    
    for k=1:sizh
        h=(80+(k-1)*0.2*2050/sizh+0.1*2050/sizh)*1000;

        if h>=465000
           if s~=1*siz/256*2
              bjk=sin(pi*(0:2*(s+1))/(2*(s+1))).^2;
              bjk=bjk/sum(bjk);
              Ap=conv2(Ap,bjk,'same');
              Ap=conv2(Ap,bjk','same');
              Bp=conv2(Bp,bjk,'same');
              Bp=conv2(Bp,bjk','same');
              A=A+Ap;
              B=B+Bp;
              Ap=zeros(siz,siz);
              Bp=zeros(siz,siz);
           end
           s=1*siz/256*2;
        elseif h<465000 && h>=235000
           if s~=2*siz/256*2
              bjk=sin(pi*(0:2*(s+1))/(2*(s+1))).^2;
              bjk=bjk/sum(bjk);
              Ap=conv2(Ap,bjk,'same');
              Ap=conv2(Ap,bjk','same');
              Bp=conv2(Bp,bjk,'same');
              Bp=conv2(Bp,bjk','same');
              A=A+Ap;
              Ap=zeros(siz,siz);
              B=B+Bp;
              Bp=zeros(siz,siz);
           end
           s=2*siz/256*2;
        elseif h<235000 && h>=155000
              if s~=3*siz/256*2
                 bjk=sin(pi*(0:2*(s+1))/(2*(s+1))).^2;
                 bjk=bjk/sum(bjk);
                 Ap=conv2(Ap,bjk,'same');
                 Ap=conv2(Ap,bjk','same');
                 Bp=conv2(Bp,bjk,'same');
                 Bp=conv2(Bp,bjk','same');
                 A=A+Ap;
                 Ap=zeros(siz,siz);
                 B=B+Bp;
                 Bp=zeros(siz,siz);
              end
              s=3*siz/256*2;
        elseif h<155000 && h>=115000
              if s~=4*siz/256*2
                 bjk=sin(pi*(0:2*(s+1))/(2*(s+1))).^2;
                 bjk=bjk/sum(bjk);
                 Ap=conv2(Ap,bjk,'same');
                 Ap=conv2(Ap,bjk','same');
                 Bp=conv2(Bp,bjk,'same');
                 Bp=conv2(Bp,bjk','same');
                 A=A+Ap;
                 Ap=zeros(siz,siz);
                 B=B+Bp;
                 Bp=zeros(siz,siz);
              end
              s=4*siz/256*2;
        elseif h<115000 && h>=95000
              if s~=5*siz/256*2
                 bjk=sin(pi*(0:2*(s+1))/(2*(s+1))).^2;
                 bjk=bjk/sum(bjk);
                 Ap=conv2(Ap,bjk,'same');
                 Ap=conv2(Ap,bjk','same');
                 Bp=conv2(Bp,bjk,'same');
                 Bp=conv2(Bp,bjk','same');
                 A=A+Ap;
                 Ap=zeros(siz,siz);
                 B=B+Bp;
                 Bp=zeros(siz,siz);
              end
              s=5*siz/256*2;
        else
              if s~=6*siz/256*2
                 bjk=sin(pi*(0:2*(s+1))/(2*(s+1))).^2;
                 bjk=bjk/sum(bjk);
                 Ap=conv2(Ap,bjk,'same');
                 Ap=conv2(Ap,bjk','same');
                 Bp=conv2(Bp,bjk,'same');
                 Bp=conv2(Bp,bjk','same');
                 A=A+Ap;
                 Ap=zeros(siz,siz);
                 B=B+Bp;
                 Bp=zeros(siz,siz);
              end
              s=6*siz/256*2;
        end
       
       for i=1:150
           for j=1:150
                
               cal=0;
               
               for l=1:sizh
                   
                   cal=cal+(tom(i,j,l)*410000/sizh)...
                       /(((80+(l-1)*0.2*2050/sizh+0.1*2050/sizh)*1000)^2);
               end
               
               if cal==0
                  continue
               end
        

               p=tom(i,j,k)*410000/sizh;
               cal=p/(h^2*cal)*10^(10);

               if p==0 
                  continue
               end
            
               if isnan(p)
                  continue
               end
            
               x=(i-1)*ds+ds/2-15000;
               y=(j-1)*ds+ds/2-15000;
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

               p=p/(r^2);
               cal=cal/(r^2);


               Xpos=ea2xy(El,Az,el,az);
               
            
               if Xpos(1)>256 || Xpos(1)<1 || Xpos(2)>256 || Xpos(2)<1
                 continue
               end

               Ap(Xpos(1),Xpos(2))=Ap(Xpos(1),Xpos(2))+p;
               Bp(Xpos(1),Xpos(2))=Bp(Xpos(1),Xpos(2))+cal;

           end

       end
    end

    bjk=sin(pi*(0:2*(s+1))/(2*(s+1))).^2;
    bjk=bjk/sum(bjk);
    Ap=conv2(Ap,bjk,'same');
    Ap=conv2(Ap,bjk','same');
    A=A+Ap;
    Bp=conv2(Bp,bjk,'same');
    Bp=conv2(Bp,bjk','same');
    B=B+Bp;

end