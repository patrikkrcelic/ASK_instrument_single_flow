function y_new = linear_interpolation_extrapolation(y_old, x_old, x_new)

% Function for linear interpolation/extrapolation. For some reason faster
% than the matlab built in function. 
%
% INPUTS:
%    y_old  - Discrete set of known points
%    x_old  - Grid of known points  
%    x_new  - Grid for interpolation/extrapolation. Can be a single value 
%             or a vector
%
% OUTPUTS:
%    y_new  - Interpolation/extrapolation values with same shape and size 
%             as x_new 
%

y_new = zeros(size(x_new));

for i=1:length(x_new)

    if x_new(i)<=x_old(1)
       a=(y_old(2)-y_old(1))/(x_old(2)-x_old(1));
       y_new(i)=y_old(1)-a*(x_old(1)-x_new(i));

    elseif x_new(i)>x_old(end)
       a=(y_old(end)-y_old(end-1))/(x_old(end)-x_old(end-1));
       y_new(i)=y_old(end)+a*(x_new(i)-x_old(end));

    else
       I=find(x_new(i)<=(x_old));
       a=(y_old(I(1))-y_old(I(1)-1))/(x_old(I(1))-x_old(I(1)-1));
       y_new(i)=y_old(I(1)-1)+a*(x_new(i)-x_old(I(1)-1));
    end
    
end

end

