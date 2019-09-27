function val = neglogpost(M,x)

% plot f(x) = eps*x + log(1+x) which is the negative log-posterior (upto
% constants) of the density proportional to exp(-eps*x)/(1+x) for x > 0

% input: x is a vector of non-negative reals

val = M.*x + rlog(x);

end