function y = rexp(x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if abs(x) < eps
   y = expm1(x);
else
   y = exp(x)-1;
end

end

