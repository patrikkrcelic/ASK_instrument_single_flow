function Xzen = zenith_estimator_model(datdir, event_setup, imnum1, imnum2, imnum3)

% Function estimates the magentic zenith position by modeling the 3 ASK
% images labeled with numbeers imnum1, imnum2, imnum3. 
%
%
% INPUTS:
%    datdir     - data directory of the event
%    event_seup - matadata of the event
%    imnum1     - number of the image 1
%    imnum2     - number of the image 2
%    imnum3     - number of the image 3
%
% OUTPUT:
%    Xzen - magnetic zenith position    
%
%


load(strcat(datdir,'/emissions.mat'), 'N2_em', 'OI_em', 'OI_D_em', 'en', 'alt');
imdir=strcat(datdir,'/ASK1/');
imdir3=strcat(datdir,'/ASK3/');

X1=importdata(strcat(datdir, '/ratio_energy.txt'));
X2=importdata(strcat(datdir, '/energy_fluxint.txt'));

ds=0.2;

siz=event_setup.siz;
sizh=41;

x_start=event_setup.zenith(1);
y_start=event_setup.zenith(2);


p0=[imnum1, imnum2, imnum3];

AA1=importdata(strcat(datdir, '/',  num2str(event_setup.mb) ,'_el_az_ask1_lut.txt'));
El=AA1(:,1:256);
Az=AA1(:,257:end);

AA3=importdata(strcat(datdir, '/', num2str(event_setup.mb), '_el_az_ask3_lut.txt'));
El3=AA3(:,1:256);
Az3=AA3(:,257:end);

Xzenith = zeros(3,2);

for br=1:3
  f=zeros(1,3);
  zen=zeros(29 ,2);
  er=zeros(29 ,1);
    
  p=p0(br);  

  sp=num2str(p,'%05.f');

  imfile=strcat(imdir,'ask1_', sp,'.txt');
  im1=importdata(imfile);
  imfile=strcat(imdir3,'ask3_', sp,'.txt');
  im3=importdata(imfile);

  im1 = image_trimming(im1, event_setup.bck1);
  A1 = im1;
  im1 = im1/max(max(im1));
  
  im3 = image_trimming(im3, event_setup.bck3);
  A3 = im3;
  im3 = im3/max(max(im3));

  count=1;
  dds=20;

  VV=[1 0; 0 1; 1 0; 0 1; 1 0; 0 1];
  X=[x_start; y_start];

  [Aen, Aflux] = EnergyFlux (A1, A3, alt, en, N2_em, OI_em+OI_D_em, X, X1, X2);
  [tom, Men, Mflux]=tomography_ASK1(Aen, Aflux, El, Az, N2_em, alt, en, siz, sizh, X, ds);
  [A, B]=projection_slow(tom, El, Az, siz, sizh, X,ds);
  improj=A./B;
  improj=improj/max(max(improj));
  tom=tomography_rest(OI_em+OI_D_em, alt,en, Men, Mflux, Aflux, 3, sizh);
  [A, B]=projection_slow(tom, El3, Az3, siz, sizh, X,ds);
  improj3=A./B;
  improj3=improj3/max(max(improj3));
  
  f(2)=sum(sum((im1(~isnan(improj))-improj(~isnan(improj))).^2))+...
      sum(sum((im3(~isnan(improj3))-improj3(~isnan(improj3))).^2));


  er(count)=f(2);

  zen(count,:)=X;
  count=count+1;
  for ind=1:6
     for lm=1:2
        
        if lm==1
            lm1=1;
            lam=-1*dds;
        else
            lm1=3;
            lam=dds;
        end
        
        xz=X+lam*VV(ind,:)'/norm(VV(ind,:));
        
        [Aen, Aflux] = EnergyFlux (A1, A3, alt, en, N2_em, OI_em+OI_D_em, xz, X1, X2);
        [tom, Men, Mflux]=tomography_ASK1(Aen, Aflux, El, Az, N2_em, alt, en, siz, sizh, xz, ds);
        [A, B]=projextion_slow(tom, El, Az, siz, sizh, xz,ds);
        improj=A./B;
        improj=improj/max(max(improj));

        tom=tomography_rest(OI_em+OI_D_em, alt,en, Men, Mflux, Aflux, 3, sizh);
        [A, B]=projection_slow(tom, El3, Az3, siz, sizh, xz,ds);
        improj3=A./B;
        improj3=improj3/max(max(improj3));

        f(lm1)=sum(sum((im1(~isnan(improj))-improj(~isnan(improj))).^2))+...
            sum(sum((im3(~isnan(improj3))-improj3(~isnan(improj3))).^2));
        
        er(count)=f(lm1);
        zen(count,:)=xz;
        count=count+1;
        
     end
    
       if f(1)==min(f)
          X=X-dds*VV(ind,:)';
          f(2)=f(1);
       elseif f(3)==min(f)
          X=X+dds*VV(ind,:)';
          f(2)=f(3);
       end
    
       if ind/2==round(ind/2)
          dds=dds/2;
       end
     
  end
  
  J=find(er~=0);
  I=find(min(er(J))==er(J));

  Xzenith(br,:)=zen(I(1),:);

end


zen1=zeros(12 ,2);
er1=zeros(12 ,1);
for br=1:3
    
  p=p0(br);  

  sp=num2str(p,'%05.f');

  imfile=strcat(imdir,'ask1_', sp,'.txt');
  im1=importdata(imfile);
  imfile=strcat(imdir3,'ask3_', sp,'.txt');
  im3=importdata(imfile);

  im1 = image_trimming(im1, event_setup.bck1);
  A1 = im1;
  im1 = im1/max(max(im1));
  
  im3 = image_trimming(im3, event_setup.bck3);
  A3 = im3;
  im3 = im3/max(max(im3));

  for k=1:4
      if k == 4
          X = round([mean(Xzenith(:,1)) mean(Xzenith(:,2))]);
      else
          X = Xzenith(k,:);
      end

      [Aen, Aflux] = EnergyFlux (A1, A3, alt, en, N2_em, OI_em+OI_D_em, X, X1, X2);
      [tom, Men, Mflux]=tomography_ASK1(Aen, Aflux, El, Az, N2_em, alt, en, siz, sizh, X, ds);
      [A, B]=projection_slow(tom, El, Az, siz, sizh, X,ds);
      improj=A./B;
      improj=improj/max(max(improj));

      tom=tomography_rest(OI_em+OI_D_em, alt,en, Men, Mflux, Aflux, 3, sizh);
      [A, B]=projection_slow(tom, El3, Az3, siz, sizh, X,ds);
      improj3=A./B;
      improj3=improj3/max(max(improj3));

      er1(br*4-(k-1)) = sum(sum((im1(~isnan(improj))-improj(~isnan(improj))).^2))+...
         sum(sum((im3(~isnan(improj3))-improj3(~isnan(improj3))).^2));
      zen1(br*4-(k-1),:) = X;
  end 
end

er2(1) = er1(1)+er1(1+4)+er1(1+8);
er2(2) = er1(2)+er1(2+4)+er1(2+8);
er2(2) = er1(3)+er1(4+4)+er1(3+8);
er2(2) = er1(4)+er1(4+4)+er1(4+8);

Xzen = zen1(min(er2)==er2,:);
Xzen = Xzen(1,:);


end


