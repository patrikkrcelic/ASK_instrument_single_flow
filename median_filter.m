function Aout=median_filter(Ain,n)

% Function calculates filtered image using median value of nxn bins for 
% each pixel. 
%
% INPUTS:
%    Ain   - Raw, unfiltered image
%    n     - Dimension of nxn filter 
%
% OUTPUT:
%    Aout  - filtered image
%    

 m=size(Ain);
 Aout=zeros(m(1),m(2));

 for i=n:m(1)-n
     for j=n:m(2)-n
         Aout(i,j)=median(median(Ain(i-n+1:i+n, j-n+1:j+n)));
     end
 end
end
        