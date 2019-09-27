clear; 
LASTN = maxNumCompThreads(11);
SAVE_SAMPLES=true;disp_int = 1000; corX = false; rhoX = .9;
delt = 1e-4; nkeep = 100; mh_sigma = true;
ApproxXLX = true; ab = false;
s_sigma = .1; 

scl_ub = .8; % scale for Metropolis-Hastings proposals for xi, use .8 except for gwas data use .5
scl_lb = .8; % lb scale, use 0.8 for sims
phasein = 1;

% ps = randsample(5001:50000,20,false);
% ns = zeros(length(ps),1);
% for j=1:length(ps)
%     ns(j) = randsample(1000:5000,1);
% end


ps = 20000.*ones(1,10); ns = 200:200:2000;

%rnd_seed = 571;
nmc = 20000; % length of Markov chain

is_sim = true;
sim_type = 'F'; % Frequentist or Bayesian
plotting = true; 

a0 = 1; b0 = 1;
BURNIN = 1000; MCMC = nmc; thin = 1; 

% 
%delt(1:5000) = 1e-3;
%delt(5001:10000) = 1e-4;
%delt(10001:15000) = 5e-5;
%delt(15001:end) = 1e-5;


ctr = 0;
for p=ps
    ctr = ctr + 1;
    n = ns(ctr);
    disp(['p: ' num2str(p) ' n: ' num2str(n)]);
    
    %rng(rnd_seed);
    
    %p = 10000; % number of parameters   

    % True parameters
    if strcmp(sim_type,'F') % 'Frequentist': BetaTrue is an unknown, deterministic, vector
      SigmaTrue = 2;
      BetaTrue = zeros(p,1);
      BetaTrue(1:23) = 2.^(-(-2:.25:3.5));
      TauTrue = 1;
    elseif strcmp(sim_type,'maize')
        load('Data/maize.mat');
        Xt = X; yt = y;
        
        p = size(X,2);
        n = size(X,1);
        BetaTrue = [];
        TauTrue = 1;
    end

    % Basic Variables
    if is_sim
        if corX
            X = normrnd(0,1,[n p]);
            for j=2:p
                X(:,j) = rhoX.*X(:,j-1)+X(:,j);
            end
            X = bsxfun(@times,X,1./std(X,[],1));
        else
            X = normrnd(0,1,[n p]);
        end
        y = X*BetaTrue+SigmaTrue.*normrnd(0,1,[n 1]);
        yt = []; Xt = [];
        save(strcat('Outputs/sim_vars_',num2str(n),'_',num2str(p),'.mat'),'X','y');
    end

    
    Sigma2Est = 1.0; % type: Float64
    TauEst = TauTrue;
    slice_lambda = true; % whether to update lambda via slice sampling    
    burn = 0; % number of burn-ins

    rng(5171);
    if ab
        [pMean,pMedian,pLambda,pSigma,betaout]=horseshoe_ab_mrg(y,X,BURNIN,MCMC,thin,1,SAVE_SAMPLES,BetaTrue,nkeep,corX,sim_type,rhoX);            
    else
        [pMean,pMedian,pLambda,pSigma,betaout,xiout,sigmaSqout,lambdaout,t]=horseshoe(y,X,BURNIN,MCMC,thin,scl_ub,scl_lb,phasein,SAVE_SAMPLES,a0,b0,BetaTrue,disp_int,...
            plotting,corX,sim_type,nkeep,ApproxXLX,is_sim,yt,Xt,delt,rhoX,mh_sigma,s_sigma);
    end

    
end





