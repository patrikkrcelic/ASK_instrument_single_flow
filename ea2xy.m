function X=ea2xy(El,Az,el,az)

% Function calculates pixel position from the known "looking" direction, using
% elevation and azimuth matrices crated in idl.
%
% INPUTS:
%    El - Matrix of pixel elevation values
%    Az - Matrix of pixel azimuth values
%    el - elevation of the "looking" direction
%    az - azimuth ot the "looking" direction
%
% OUTPUT:
%    X - 2D pixel position vector 
%    

  if el<min(min(El)) | el >max(max(El))
      X=[-1 -1];
      return
  end

  if az<min(min(Az)) | az >max(max(Az))
      X=[-1 -1];
      return
  end
 
  I=find(El>=el & El<el+0.001 & Az>=az & Az<az+0.001);

  if length(I)<1

      I=find(El>=el & El<el+0.003 & Az>=az & Az<az+0.003);

      if length(I)<1
         X=[-1 -1];
         return
      end
  end

  r=10;
  if length(I)>1
      for i=1:length(I)
          r1=sqrt((Az(I(i))-az)^2+(El(I(i))-el)^2);
          if r1<r
              ij=i;
              r=r1;
          end
      end
      I=I(ij);
  end

  AE=zeros(256,256);
  AE(I)=1;

  for i=1:256
      J=find(AE(i,:)==1); 
      if length(J)>=1
          X(2)=(J(1));
          X(1)=i; 
          return
      end
  end
  
end