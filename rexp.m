function y = rexp(x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
y = zeros(length(x),1);
if any(abs(x) < eps)
   y(abs(x)<eps) = expm1(x(abs(x)<eps));
end
if any(abs(x)>=eps)
   y(abs(x)>=eps) = exp(x(abs(x)>=eps))-1;
end

end

