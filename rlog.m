function y = rlog(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if abs(x) < eps
   y = log1p(x); 
else
   y = log(1+x); 
end


end

