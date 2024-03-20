function en = ratio2energy(r, X)

% Function calulates precipiting electron energy form the energy-ratio 
% curvature, created in idl. Can also be used for flux coeficient
% with flux-energy ratio.
%
% INPUTS:
%    r     - ASK3/ASK1 ratio
%    X     - Energy-ratio curvature values
%
% OUTPUTS:
%    en    - Precipitating electron energy 
%
  
if isnan(r)
    en=NaN;
    return
end

en = linear_interpolation_extrapolation(X(:,2), X(:,1), r);

if en<0
    en=X(end,2);
end

end
