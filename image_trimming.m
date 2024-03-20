function Aout = image_trimming(Ain, background)

% Function cuts upper 0.1% and lower 15% from ASK images. Additionaly 
% function substracts background emissions.
%
% INPUTS:
%    Ain         - Original ASK image
%    background  - Value for baskground substraction
%
% OUTPUT:
%    Aout    - Trimmed image
%

   s = size(Ain);
   AA=reshape(Ain,[s(1)*s(2),1]);
   Amax = prctile(AA, 99);
   Amin = prctile(AA, 15);
   Ain(Ain>Amax)=Amax;
   Ain(Ain<Amin)=Amin;

   Aout = Ain - background;

end



