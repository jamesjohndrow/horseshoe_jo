function [samp,numdraw] = HSeta1(M)
% implement ordinary rejection sampling idea to sample from the density
% f(x) \propto exp(-Mx)/(1+x), x > 0

n = length(M);
breakloop = false; numdraw = ones(n,1);

xout = zeros(n,1); nid = (1:n)';

while ~breakloop
    v = unifrnd(0,1,[n 1]); x = -rlog(-v)./M;   % x \sim Exp(M)
    u = unifrnd(0,1,[n 1]);
    
    done = u<1./(1+x);
    xout(nid(done)) = x(done);
    numdraw(nid(~done)) = numdraw(nid(~done))+1;
    
    M = M(~done); n = sum(~done); nid = nid(~done);
    
    breakloop = all(done);
    
end
samp = xout;


end

