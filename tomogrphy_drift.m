function [n11 n22]=tomogrphy_drift(n1,n2,vv,dt,ds,sizh)

% Function calculates ion drift effents on the O+ density volume.
%
% INPUTS:
%    n1   - Initial O+(2P) 1/2 upper state density
%    n2   - Initial O+(2P) 3/2 upper state density
%    vv   - 2D vector drift velocity
%    dt   - Time resolution
%    ds   - Horizontal spatial resolution of the volume in km
%    sizh - Number of points in volume z component
%
% OUTPUTS:
%    n11  - Output O+(2P) 1/2 upper state density
%    n22  - Output O+(2P) 3/2 upper state density
%


  dx=vv(1)*dt/ds;
  dy=vv(2)*dt/ds;
  if isnan(dx)
    dx=0;
  end
  if isnan(dy)
    dy=0;
  end
  s=size(n1);
  ix=1:s(1);
  iy=1:s(2);

  while dx < 0                    
     ix = [1,ix(1:end-1)];
     dx = dx+1;
  end
  while dx > 1                    
     ix = [ix(2:end),ix(end)];
     dx = dx-1;
  end

  while dy < 0
     iy = [1,iy(1:end-1)];
     dy = dy+1;
  end
  while dy > 1
     iy = [iy(2:end),iy(end)];
     dy = dy-1;
  end
      
  for k=1:sizh
    n11(:,:,k) = ( n1(iy(1:end),ix(1:end),k)   *(1-dx) + ...
    n1(iy(1:end),ix([2:end,end]),k) * dx ) * (1-dy) + ...    
    (n1(iy([2:end,end]),ix(1:end),k)    * (1-dx) + ...         
    n1(iy([2:end,end]),ix([2:end,end]),k)*dx    )   * dy;
    
    n22(:,:,k) = ( n2(iy(1:end),ix(1:end),k)   *(1-dx) + ...
    n2(iy(1:end),ix([2:end,end]),k) * dx ) * (1-dy) + ...  
    (n2(iy([2:end,end]),ix(1:end),k)    * (1-dx) + ...         
    n2(iy([2:end,end]),ix([2:end,end]),k)*dx    )   * dy;
    
  end   
end
      
      
      
      
      
      
      
      
      
      
      
      