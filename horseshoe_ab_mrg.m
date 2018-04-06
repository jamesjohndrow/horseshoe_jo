

function [pMean,pMedian,pLambda,pSigma,betaout]=horseshoe_ab_mrg(y,X,BURNIN,MCMC,thin,type,SAVE_SAMPLES,BetaTrue,nkeep,corX,simtype,rhoX)
% Function to impelement Horseshoe shrinkage prior (http://faculty.chicagobooth.edu/nicholas.polson/research/papers/Horse.pdf)
% in Bayesian Linear Regression. %%
% Written by Antik Chakraborty (antik@stat.tamu.edu) and Anirban Bhattacharya (anirbanb@stat.tamu.edu)

% Model: y=X\beta+\epslion, \epsilon \sim N(0,\sigma^2) %%
%        \beta_j \sim N(0,\sigma^2 \lambda_j^2 \tau^2) %%
%        \lambda_j \sim Half-Cauchy(0,1), \tau \sim Half-Cauchy (0,1) %%
%        \pi(\sigma^2) \sim 1/\sigma^2 %%

if corX
   simtype = strcat(simtype,'_cor_',strrep(num2str(rhoX*10),'.','_')); 
end




% This function employs the algorithm proposed in "Fast Sampling with Gaussian Scale-Mixture priors in High-Dimensional Regression" by
% Bhattacharya et. al. (2015). The global local scale parameters are updated via a Slice sampling scheme given in the online supplement 
% of "The Bayesian Bridge" by Polson et. al. (2011). Two different algorithms are used to compute posterior samples of the p*1 vector
% of regression coefficients \beta. Type=1 corresponds to the proposed method in Bhattacharya et. al. and Type=2 corresponds to an algorithm provided
% in Rue (2001). We recomend our algorithm when p>n.%%


% Input: y=response, a n*1 vector %%
%        X=matrix of covariates, dimension n*p %%
%        BURNIN= number of burnin MCMC samples %%
%        MCMC= number of posterior draws to be saved %%
%        thin= thinning parameter of the chain %%
%        type= 1 if sampling from posterior of beta is done by Proposed Method %%
%              2 if sampling is done by Rue's algorithm %%
%        SAVE_SAMPLES= binary indicator whethe posterior samples should be saved in a file or not %%


% Output: 
%         pMean= posterior mean of Beta, a p by 1 vector%%
%         pMeadian=posterior median of Beta, a p by 1 vector %%
%         pLambda=posterior mean of local scale parameters, a p by 1 vector %%
%         pSigma=posterior mean of Error variance %% 
%         betaout=posterior samples of beta %%



if SAVE_SAMPLES
   %fprintf('Warning: program will save posterior beta samples to the current working directory \n');
   %dum = input('Enter 1 to proceed or any other key to exit ');
   %if dum~=1
   %   return;
   %end
end


%profile on;


tic;
N=BURNIN+MCMC;
effsamp=(N-BURNIN)/thin;
[n,p]=size(X);

% paramters 
Beta=zeros(p,1); lambda=ones(p,1);
tau=1; sigma_sq=1; Xi = 1./tau.^2;

% output 
betaout=zeros(nkeep,effsamp);
lambdaout=zeros(nkeep,effsamp);
tauout=zeros(effsamp,1);
sigmaSqout=zeros(effsamp,1);
l1out = zeros(effsamp,1);
pexpout = zeros(effsamp,1);

% matrices 
I_n=eye(n); 
l=ones(n,1);
if type==2
   Q_star=X'*X;
end


% start Gibb's sampling 
for i=1:N
    
    LX=bsxfun(@times,(lambda.^2),X');
    XLX = X*LX;
    M = I_n + (1./Xi).*XLX;
    x = M\y;

    ssr = y'*x; ssr = max(ssr,1e-10);
    sigma_sq = 1/gamrnd((n+1)/2,2/(ssr+1));


    % changed prior
    U = (1./Xi).*LX;    
    
    % step 1 %
    u=normrnd(0,tau*lambda);
    v=X*u+normrnd(0,l);

    v_star=(M)\((y./sqrt(sigma_sq))-v);
    Beta=sqrt(sigma_sq)*(u+U*v_star);


    % update lambda_j's in a block using slice sampling
    eta = 1./(lambda.^2); 
    upsi = unifrnd(0,1./(1+eta));
    tempps = Beta.^2/(2*sigma_sq*tau^2); 
    ub = (1-upsi)./upsi;

    % now sample eta from exp(tempv) truncated between 0 & upsi/(1-upsi)
    Fub = 1 - exp(-tempps.*ub); % exp cdf at ub 
    Fub(Fub < (1e-4)) = 1e-4;  % for numerical stability
    up = unifrnd(0,Fub); 
    eta = -log(1-up)./tempps; 
    lambda = 1./sqrt(eta);

    % update tau 
    tempt = sum((Beta./lambda).^2)/(2*sigma_sq); 
    et = 1/tau^2; 
    utau = unifrnd(0,1/(1+et));
    ubt = (1-utau)/utau; 
    Fubt = gamcdf(ubt,(p+1)/2,1/tempt); 
    Fubt = max(Fubt,1e-8); % for numerical stability
    ut = unifrnd(0,Fubt); 
    et = gaminv(ut,(p+1)/2,1/tempt); 
    tau = 1/sqrt(et);
    Xi = 1./tau.^2;

    % update sigma_sq

    if mod(i,1000) == 0
        disp(i)
    end
    
    per_expl = 1-sqrt(sum((Beta-BetaTrue).^2))./sqrt(sum(BetaTrue.^2));
    L1_loss = 1-sum(abs(Beta-BetaTrue))./sum(abs(BetaTrue));
    
    if i > BURNIN && mod(i, thin)== 0
        betaout(:,(i-BURNIN)/thin) = Beta(1:nkeep);
        lambdaout(:,(i-BURNIN)/thin) = lambda(1:nkeep);
        tauout((i-BURNIN)/thin)=tau;
        sigmaSqout((i-BURNIN)/thin)=sigma_sq;
        l1out(i) = L1_loss;
        pexpout(i) = per_expl;
    end
end
pMean=mean(betaout,2);
pMedian=median(betaout,2);
pSigma=mean(sigmaSqout);
pLambda=mean(lambdaout,2);
t=toc;

ci_lo = quantile(betaout(:,5001:end)',.025,1);
ci_hi = quantile(betaout(:,5001:end)',.975,1);
coverage = mean(BetaTrue(1:nkeep)>ci_lo' & BetaTrue(1:nkeep)<ci_hi');
BetaHat = mean(betaout(:,5001:end),2)';
mse = mean((BetaTrue(1:nkeep)-BetaHat').^2);
se = std(betaout(:,5001:end),[],2);

disp(['coverage ' num2str(100*coverage)]);
disp(['mse ' num2str(mse)]);

%fprintf('Execution time of %d Gibbs iteration with (n,p)=(%d,%d)is %f seconds',N,n,p,t)
if SAVE_SAMPLES
    etaout = 1./lambdaout.^2;
    xiout = 1./tauout.^2;
    save(strcat('Outputs/post_reg_horse_ab_mrg_',simtype,'_',num2str(n),'_',num2str(p),'.mat'),'betaout','lambdaout','etaout','tauout','xiout','sigmaSqout','t',...
        'ci_hi','ci_lo','coverage','BetaHat','mse','se','BetaTrue');
end


figure(1);subplot(2,2,1);plot(log(xiout(50:(i-1))),'.');title('log(\xi)');
subplot(2,2,2);plot(log(1./sigmaSqout(50:(i-1))),'.');title('log(\sigma^{-2})');
subplot(2,2,3);plot(betaout(1,50:(i-1)),'.');title('\beta_1');
subplot(2,2,4);plot(etaout(1,50:(i-1)),'.');title('\eta_1');
drawnow;
figure(2); subplot(1,2,1);plot(l1out(1:(i-1)),'.');title('L1');ylim([0 1]);
subplot(1,2,2);plot(pexpout(1:(i-1)),'.');title('pexpl'); ylim([0 1]);
drawnow;
figure(3);
for j0=1:25
    subplot(5,5,j0); plot(betaout(j0,50:(i-1)),'.'); title(strcat('\beta_{',num2str(j0),'}'));
end
drawnow;



%profile viewer;
