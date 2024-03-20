function [A, B]=projection_fast(indir, camera, tom, siz, sizh)

% Function does fast projection from tomography volume to ASK image,
% using projection 'look up tables' created with projection_trace function.
%
% INPUTS:
%    indir  - Path to the directory with the projection 'look up tables'
%    camera - ASK camera number
%    tom    - Volume emission rates 3D volume
%    siz    - ASK image size
%    sizh   - Number of points in volume z component
%
% OUTPUTS:
%    A  - Uncalibrated projected image
%    B  - Calibration image
%


    A=zeros(siz,siz);
    Ap=zeros(siz,siz);
    B=zeros(siz,siz);
    Bp=zeros(siz,siz);
    
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

        load(strcat(indir, '/camera_',num2str(camera), '_XY_', num2str(k), '.mat'), 'imat', 'jmat', 'rad');
       
       for i=1:150
           for j=1:150

               if tom(i,j,k)==0
                   continue
               end

               if imat(i,j)==0 || jmat(i,j)==0
                   continue
               end
                
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

               p=p/(rad(i,j)^2);
               cal=cal/(rad(i,j)^2);
            

               Ap(imat(i,j),jmat(i,j))=Ap(imat(i,j),jmat(i,j))+p;
               Bp(imat(i,j),jmat(i,j))=Bp(imat(i,j),jmat(i,j))+cal;

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