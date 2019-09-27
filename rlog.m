function y = rlog(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
y = zeros(length(x),1);
if any(abs(x) < eps)
   y(abs(x)<eps) = log1p(x(abs(x)<eps)); 
end
if any(abs(x)>=eps)
   y(abs(x)>=eps) = log(1+x(abs(x)>=eps)); 
end


end

