using Distributions
include("horseshoe.jl")
SAVE_SAMPLES=true; disp_int = 1000; corX = false; rhoX = 0.9;
delt = 1e-4; nkeep = 100; mh_sigma = true;
ApproxXLX = true; ab = false;
s_sigma = 0.1;

scl = 0.8;

p = 20000; n = 2000;

rnd_seed = 571;
nmc = 20000; # length of Markov chain

is_sim = true;
sim_type = "F"; # Frequentist or Bayesian
plotting = true;

a0 = 1.; b0 = 1.;
BURNIN = 0; MCMC = nmc; thin = 1;

ctr = 0;


println("p: $(p), n: $(n)");
srand(rnd_seed);


# True parameters
if sim_type=="F" # 'Frequentist': BetaTrue is an unknown, deterministic, vector
  SigmaTrue = 2;
  BetaTrue = zeros(p,1);
  BetaTrue[1:23] = 2.^(-range(-2,.25,23));
  TauTrue = 1;
else
  #load('Data/gwas/maize/maize.mat');
  #Xt = X; yt = y;

  #p = size(X,2);
  #n = size(X,1);
  #BetaTrue = [];
  #TauTrue = 1;
end

# Basic Variables
if is_sim
    if corX
        X = rand(Normal(0.,1.),n,p);
        for j=2:p
            X[:,j] = rhoX.*X[:,j-1]+X[:,j];
        end
        X = broadcast(*,X,1./std(X,1));
    else
        X = rand(Normal(),n,p);
    end
    y = X*BetaTrue+SigmaTrue.*rand(Normal(0.,1.),n,1);
    yt = zeros(2,1); Xt = zeros(2,1);
end


Sigma2Est = 1.0; # type: Float64
TauEst = TauTrue;
slice_lambda = true; # whether to update lambda via slice sampling
burn = 0; # number of burn-ins
simtype = sim_type;

bet = horseshoe(y,X,BURNIN,MCMC,scl,SAVE_SAMPLES,a0,b0,BetaTrue,disp_int,plotting,corX,simtype,nkeep,ApproxXLX,is_sim,yt,Xt,delt,rhoX,mh_sigma,s_sigma)

#[pMean,pMedian,pLambda,pSigma,betaout,xiout,sigmaSqout,lambdaout,t]=horseshoe(y,X,BURNIN,MCMC,thin,scl_ub,scl_lb,phasein,SAVE_SAMPLES,a0,b0,BetaTrue,disp_int,...
#          plotting,corX,sim_type,nkeep,ApproxXLX,is_sim,yt,Xt,delt,rhoX,mh_sigma,s_sigma);
