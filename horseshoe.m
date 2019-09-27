function[Beta_hat,pMedian,pLambda,pSigma,betaout,xiout,sigmaSqout,lambdaout,t]=horseshoe(y,X,BURNIN,MCMC,thin,scl_ub,scl_lb,phasein,SAVE_SAMPLES,a0,b0,BetaTrue,disp_int,plotting,corX,simtype,nkeep,ApproxXLX,...
    is_sim,yt,Xt,delt,rhoX,mh_sigma,s_sigma)
% Function to impelement Horseshoe shrinkage prior (http://faculty.chicagobooth.edu/nicholas.polson/research/papers/Horse.pdf)
% in Bayesian Linear Regression. %%
% Based on code by Antik Chakraborty (antik@stat.tamu.edu) and Anirban Bhattacharya (anirbanb@stat.tamu.edu)
% Modified by James Johndrow (johndrow@stanford.edu)

% Model: y=X\beta+\epslion, \epsilon \sim N(0,\sigma^2) %%
%        \beta_j \sim N(0,\sigma^2 \lambda_j^2 \tau^2) %%
%        \lambda_j \sim Half-Cauchy(0,1), \tau \sim Half-Cauchy (0,1) %%
%        \pi(\sigma^2) \sim 1/\sigma^2 %%


% This function implements the algorithm proposed in "Scalable MCMC for
% Bayes Shrinkage Priors" by Johndrow and Orenstein (2017). 
% The global local scale parameters are updated via a Slice sampling scheme given in the online supplement 
% of "The Bayesian Bridge" by Polson et. al. (2011). Setting ab = true
% implements the algorithm of Bhattacharya et al. Setting ab=false
% implements the algorith of Johndrow and Orenstein, which uses a block
% update for \tau, \sigma^2, \beta



% Input: y=response, a n*1 vector %%
%        X=matrix of covariates, dimension n*p %%
%        BURNIN= number of burnin MCMC samples %%
%        MCMC= number of posterior draws to be saved %%
%        thin= thinning parameter of the chain %%
%        scl_ub=upper bound on scale for MH proposals
%        scl_lb=lower bound on scale for MH proposals (usually make these
%        equal; 0.8 a good default)
%        phasin=number of iterations over which to transition between upper
%        and lower bound on MH proposals; usually make this 1 and just make
%        scl_ub=scl_lb; only use for particularly challenging cases
%        SAVE_SAMPLES=binary; whether to save samples
%        ab=whether to run the algorithm of Bhattachaya et al.
%        trunc=whether to use the numeric truncations of Bhattacharya et al
%        a0=parameter of gamma prior for sigma2
%        b0=second parameter of gamma prior for sigma2
%        BetaTrue=true beta (for simulations)
%        disp_int=how often to produce output
%        plotting=whether to make plots
%        corX=whether simulations were performed with correlated design



% Output: 
%         pMean= posterior mean of Beta, a p by 1 vector%%
%         pMeadian=posterior median of Beta, a p by 1 vector %%
%         pLambda=posterior mean of local scale parameters, a p by 1 vector %%
%         pSigma=posterior mean of Error variance %% 
%         betaout=posterior samples of beta %%



%profile on;

if ApproxXLX
    if length(delt)==1
        simtype = strcat(simtype,'_approx_',num2str(-log10(delt)));
    else
        simtype = strcat(simtype,'_approx_dseq');
    end
end
if corX
   simtype = strcat(simtype,'_cor_',strrep(num2str(rhoX*10),'.','_')); 
end

tic;
N=BURNIN+MCMC;
if length(delt)==1
   delt = delt.*ones(N,1); 
end
effsamp=(N-BURNIN)/thin;
[n,p]=size(X);
%[U,D,V] = svd(X,'econ');
%Dd = diag(D);
%dV = bsxfun(@times,Dd',V);

% paramters %
Beta=ones(p,1); 
lambda=ones(p,1);
tau=1; sigma_sq=1;

% output %
betaout=zeros(nkeep,effsamp);
lambdaout=zeros(nkeep,effsamp);
etaout = zeros(nkeep,effsamp);
tauout=zeros(effsamp,1);
xiout = zeros(MCMC+BURNIN,1);
sigmaSqout=zeros(effsamp,1);
l1out = zeros(MCMC+BURNIN,1);
pexpout = zeros(MCMC+BURNIN,1);
ACC = zeros(MCMC+BURNIN,1);
ACC_s = zeros(MCMC+BURNIN,1);

% matrices %
if issparse(X)
   I_n = speye(n); 
   %s_max = svds(X,1);
else
    I_n=eye(n); 
end
l=ones(n,1);

%Xy = X'*y;

Xi = tau^(-2);
Eta = lambda.^(-2);
% lambda_old = lambda;

if ~is_sim
   Beta_hat = zeros(p,1); 
end

% start Gibbs sampling %
t = 0;
tic;
for i=1:N      
    if mod(i,disp_int)==0
        disp(num2str(i)); dt = toc; tic;
        t = t+dt;
        disp([num2str(disp_int) ' iterations in ' num2str(dt) ' seconds']);
        %Beta(1:20)'        
    end
    %disp(num2str(i));
    % update tau %
    if i>0
        Eta = lambda.^(-2);
        if i<phasein
            std_MH = (scl_ub.*(phasein-i)+scl_lb.*i)./phasein;
        else
            std_MH = scl_lb;
        end
        prop_Xi = exp(normrnd(log(Xi),std_MH));
        %prop_Xi = exp(trnd(2).*std_MH + log(Xi));
        
        if ApproxXLX
            which_in = lambda.^2.*max(1/Xi,1/prop_Xi)>delt(i);
            id1 = (1:p)'; id1 = id1(which_in);
            rXLX = sum(which_in);
            %disp(num2str(rXLX));
            
            lambda1 = lambda(id1); X1 = X(:,id1);
            
            LX1 = bsxfun(@times,lambda1.^2,X1');           
            LX = LX1;

            if rXLX < n/2
                XX1 = X1'*X1;
            else
                XLX = X1*LX1;
            end
        else
            LX=bsxfun(@times,(lambda.^2),X');
            XLX = X*LX;
            rXLX = n;
        end
        
        
        
        
        if rXLX<n/2 && ApproxXLX
            M_prop = [];
            XL0 = bsxfun(@times,X1,lambda1');
            s = diag(chol(XL0'*XL0));
            %s = svd(full(bsxfun(@times,X1,lambda1')));
            s_prop = s.*(1./sqrt(prop_Xi));            
            cM_prop = sqrt(1+s_prop.^2);            
            M = [];
            s_curr = s.*(1./sqrt(Xi));
            cM = sqrt(1+s_curr.^2);
            x_prop = y - X1*((prop_Xi.*diag(lambda1.^(-2))+XX1)\(X1'*y));
            x = y - X1*((Xi.*diag(lambda1.^(-2))+XX1)\(X1'*y));
        else
            M = I_n + (1./Xi).*XLX;
            M_prop = I_n + (1./prop_Xi).*XLX;
            if issparse(X)
                per = symamd(M_prop);
                cM_prop = diag(chol(M_prop(per,per),'lower'));
                cM = diag(chol(M(per,per),'lower'));
            else
                cM_prop = diag(chol(M_prop));
                cM = diag(chol(M));
            end
            x_prop = M_prop\y;
            x = M\y;
        end
        

        [lrat_prop,~] = lmh_ratio(y,x_prop,prop_Xi,cM_prop,n,a0,b0);
        [lrat_curr,~] = lmh_ratio(y,x,Xi,cM,n,a0,b0);
        log_acc_rat = (lrat_prop-lrat_curr)+(log(prop_Xi)-log(Xi));
              
 
        ACC(i) = (rand < exp(log_acc_rat));
        if ACC(i) % if accepted, update
            Xi = prop_Xi;
            M = M_prop;
            x = x_prop;
        end
        tau = 1./sqrt(Xi);
    end
    
    % update sigma_sq %

        
    % metropolis option
    if mh_sigma
        ssr = y'*x;
        prop_k = exp(normrnd(log(1./sigma_sq),s_sigma));
        %prop_k

        l_prop = (n+a0-2)/2*log(prop_k)-prop_k*(ssr+b0)/2;
        l_curr = (n+a0-2)/2*log(1./sigma_sq)-(1/sigma_sq)*(ssr+b0)/2;
        l_ar = (l_prop-l_curr) + (log(prop_k)-log(1/sigma_sq));
        acc_s = rand<exp(l_ar);
        if acc_s
            sigma_sq = 1/prop_k;
        end
        ACC_s(i) = acc_s;
    else
        ssr = y'*x;
        sigma_sq = 1/gamrnd((n+a0)/2,2/(ssr+b0));
    end
    
    U = (1./Xi).*LX;    
    
    % step 1 %
    u=normrnd(0,tau*lambda);
    v=X*u+normrnd(0,l);
    
    % step 2 %
    
    if ApproxXLX && rXLX<n/2 % woodbury
        tmp = (Xi.*diag(lambda1.^(-2))+XX1)\(X1'*((y./sqrt(sigma_sq))-v));
        v_star = ((y./sqrt(sigma_sq))-v)-X1*tmp;
        
        Beta = sqrt(sigma_sq)*u;
        Beta(id1) = Beta(id1) + sqrt(sigma_sq)*(U*v_star);
        
    elseif ApproxXLX && rXLX>=n/2 % direct solve
         v_star=(M)\((y./sqrt(sigma_sq))-v);        
         Beta = sqrt(sigma_sq)*u;
         Beta(id1) = Beta(id1) + sqrt(sigma_sq)*(U*v_star);        
    else
        v_star=(M)\((y./sqrt(sigma_sq))-v);
        Beta=sqrt(sigma_sq)*(u+U*v_star);
    end        

    if ~is_sim
        if i < floor(BURNIN/2)
            Beta_hat = (i-1)./i.*Beta_hat + 1./i.*Beta;
        elseif i==floor(BURNIN/2)
            Beta_hat = Beta;
            se = Beta.^2;
        else
            idrem = floor(BURNIN/2);
            Beta_hat = (i-idrem-1)./(i-idrem).*Beta_hat + 1./(i-idrem).*Beta;
            se = (i-idrem-1)./(i-idrem).*se + 1./(i-idrem).*Beta.^2;
        end
    end

    %u = unifrnd(0, 1./(Eta+1));
    gamma_rate = (Beta.^2) .* Xi ./ (2.*sigma_sq);
    %Eta = gen_truncated_exp(gamma_rate, (1-u)./u);
    Eta = HSeta_rs_smallM(gamma_rate);
    if any(Eta<=0)
        disp([num2str(sum(Eta<=0)) ' Eta underflowed, replacing = machine epsilon']);
        Eta(Eta<=0) = eps;
    end
    lambda = 1./sqrt(Eta);

    %disp(num2str(i));
    if mod(i,disp_int) == 0 
        disp(num2str(i));    
        mean(ACC_s(1:(i-1)))
        disp(mean(ACC(1:i)));
        disp(['rank of XLX: ' num2str(rXLX)]);
        if is_sim
            Beta(1:20)        
        else
            Bet_ord = sort(Beta_hat,'descend');
            Bet_ord(1:20)
        end
        if plotting && i > BURNIN
            id_upper = i-BURNIN-1;
            figure(1);subplot(2,2,1);plot(log(xiout(50:id_upper)),'.');title('log(\xi)');
            subplot(2,2,2);plot(log(1./sigmaSqout(50:id_upper)),'.');title('log(\sigma^{-2})');
            subplot(2,2,3);plot(betaout(1,50:id_upper),'.');title('\beta_1');
            subplot(2,2,4);plot(etaout(1,50:id_upper),'.');title('\eta_1');
            drawnow;
            figure(2); subplot(1,2,1);plot(l1out(1:id_upper),'.');title('L1');ylim([0 1]);
            subplot(1,2,2);plot(pexpout(1:id_upper),'.');title('pexpl'); ylim([0 1]);
            drawnow;
            figure(3);
            for j0=1:25
                subplot(5,5,j0); plot(betaout(j0,50:id_upper),'.');
            end
            drawnow;
            figure(5); hist(log10((1./Xi).*lambda.^2),100); title('\xi^{-1} \lambda^2_j');drawnow;
            
            if ~is_sim
                figure(4);
                for j0=1:6
                    subplot(2,3,j0); plot(betaout(j0,50:id_upper),'.');
                end
                drawnow;
            end
        end
    end
    
    if is_sim
        per_expl = 1-sqrt(sum((Beta-BetaTrue).^2))./sqrt(sum(BetaTrue.^2));        
        L1_loss = 1-sum(abs(Beta-BetaTrue))./sum(abs(BetaTrue));
    else
        res1 = yt-Xt(:,id1)*Beta_hat(id1);
        per_expl = 1-sqrt(sum(res1.^2))./sqrt(sum(yt.^2));
        L1_loss = 1-sum(abs(res1))./sum(abs(yt));        
    end
    
    if i==BURNIN || (i==1 && BURNIN==0)
        if ~is_sim
            [~,big_id] = sort(abs(Beta_hat),'descend');
            big_id = big_id(1:nkeep/2);
            other_id = setdiff((1:p)',big_id);
            other_id = other_id(1:nkeep/2);
            keep_id = [big_id;other_id];
        else
            keep_id = (1:nkeep)';
        end
    end
    
    if i > BURNIN && mod(i, thin)== 0        
        betaout(:,(i-BURNIN)/thin) = Beta(keep_id);
        lambdaout(:,(i-BURNIN)/thin) = lambda(keep_id);
        etaout(:,(i-BURNIN)/thin) = Eta(keep_id);
        tauout((i-BURNIN)/thin)=tau;
        xiout((i-BURNIN)/thin) = Xi;
        sigmaSqout((i-BURNIN)/thin)=sigma_sq;
        l1out((i-BURNIN)/thin) = L1_loss;
        pexpout((i-BURNIN)/thin) = per_expl;
    end
end
pMean=mean(betaout,2);
pMedian=median(betaout,2);
pSigma=mean(sigmaSqout);
pLambda=mean(lambdaout,2);
dt=toc; t = t+dt;

ci_lo = quantile(betaout(:,5001:end)',.025,1);
ci_hi = quantile(betaout(:,5001:end)',.975,1);

if is_sim
    coverage = mean(BetaTrue(1:nkeep)>ci_lo' & BetaTrue(1:nkeep)<ci_hi');
    BetaHat = mean(betaout(:,5001:end),2)';
    Beta_hat = BetaHat;

    mse = mean((BetaTrue(1:nkeep)-BetaHat').^2);
    se = std(betaout(:,5001:end),[],2);

    disp(['coverage ' num2str(100*coverage)]);
    disp(['mse ' num2str(mse)]);
    keep_id = [];
else
   coverage = []; mse = [];
   BetaHat = Beta_hat;
end


disp([num2str(t) ' seconds elapsed']);
if SAVE_SAMPLES
    save(strcat('Outputs/post_reg_horse_conc_',simtype,'_',num2str(n),'_',num2str(p),'.mat'),'betaout','lambdaout','etaout','tauout','xiout','sigmaSqout','l1out','pexpout','t',...
        'ci_hi','ci_lo','coverage','BetaHat','mse','se','BetaTrue','keep_id');
end


figure(1);subplot(2,2,1);plot(log(xiout(50:(i-BURNIN-1))),'.');title('log(\xi)');
subplot(2,2,2);plot(log(1./sigmaSqout(50:(i-BURNIN-1))),'.');title('log(\sigma^{-2})');
subplot(2,2,3);plot(betaout(1,50:(i-BURNIN-1)),'.');title('\beta_1');
subplot(2,2,4);plot(etaout(1,50:(i-BURNIN-1)),'.');title('\eta_1');
drawnow;
figure(2); subplot(1,2,1);plot(l1out(1:(i-BURNIN-1)),'.');title('L1');ylim([0 1]);
subplot(1,2,2);plot(pexpout(1:(i-BURNIN-1)),'.');title('pexpl'); ylim([0 1]);
drawnow;
figure(3);
for j0=1:25
    subplot(5,5,j0); plot(betaout(j0,50:(i-BURNIN-1)),'.'); title(strcat('\beta_{',num2str(j0),'}'));
end
drawnow;

end



function [lr,x] = lmh_ratio(y,x,Xi,cM,n,a0,b0)
    % marginal of beta, sigma2    
    ssr = y'*x+b0;
    try
        ldetM = 2*sum(log(cM));        
        ll = -.5.*ldetM - ((n+a0)/2).*log(ssr);
        lpr = -log(sqrt(Xi).*(1+Xi));
        lr = ll+lpr;
    catch
        lr = -Inf; warning('proposal was rejected because I+XDX was not positive-definite');
    end
end


function x = gen_truncated_exp(mn,trunc_point)
    r = mn.*trunc_point;
    sml = abs(r)<eps;
    x = zeros(max(length(mn),length(trunc_point)),1);
    tmp = zeros(max(length(mn),length(trunc_point)),1);
    
    if any(sml)
        tmp(sml) = expm1(-r(sml)).*rand(length(mn(sml)),1);
    end
    tmp(~sml) = (exp(-r(~sml))-1).*rand(length(mn(~sml)),1);
    
    sml = abs(tmp)<eps;
    if any(sml)
       x(sml) = -log1p(tmp(sml))./mn(sml); 
    end
    x(~sml) = -log(1+tmp(~sml))./mn(~sml);
    
end

function x = samp_eta(M)
    n = length(M);
    cont = true(n,1);
    x = zeros(n,1);   
    while any(cont)  
        n = sum(cont);
        x(cont) = exprnd(M(cont),[n 1]);
        u = unifrnd(0,1,[n 1]);
        y = (1-u)./u;
        cont(cont) = (y>=x(cont));
    end
end



