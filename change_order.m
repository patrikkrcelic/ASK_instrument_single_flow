function b = change_order(a)

% Function chages the order of reading the vector
%
% INPUTS
%      a - input vector
%
% OUTPUT:
%      b - output vecotr with reverse reading order
%
   l=length(a);
   
   for i=1:l
       b(i)=a(l-i+1);
   end

end