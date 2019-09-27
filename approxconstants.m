function [nu, nuvec] = approxconstants(M,a,b)

% calculate the normalizing constant nu of the approximation h_L (see
% notes) to h(x) \propto exp(-eps*x)/(1+x) for x > 0

% Input: 
% eps < 1 determines the density h
% a < 1 determines the first and second segment of the approximation
% b > 1 determines the third and final segment of the approximation

% Output
% nu is the normalizing constant of h_L (see notes)
% nuvec is a 4x1 vector containing the contribution from each segment so
% that nu = sum(nuvec)
% nutrue = exp(eps)*expint(eps) is the normalizing constant of h

% segment 1

nu1 = rlog(a./M);

% segment 2
lambda2 = M.*( neglogpost(M,1./M) - neglogpost(M,a./M) )./(1-a);
nu2 = (1./lambda2).*( exp(-neglogpost(M,a./M)) - exp(-neglogpost(M,1./M)) );

% segment 3
lambda3 = M.*( neglogpost(M,b./M) - neglogpost(M,1./M) )./(b-1);
nu3 = (1./lambda3).*( exp(-neglogpost(M,1./M)) - exp(-neglogpost(M,b./M)) );

% segment 4
nu4 = (1./M).*exp( -neglogpost(M,b./M) );

nu = nu1+nu2+nu3+nu4;
nuvec = [nu1 nu2 nu3 nu4];
%nutrue = exp(M).*expint(M);

end

