function y = rexp(x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if abs(y) < eps
   y = expm1(y);
else
   y = exp(y)-1;
end

end

