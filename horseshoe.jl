
using Distributions
using PyPlot
using CSV
using DataFrames

function lmh_ratio(y::Array{Float64,2},x::Array{Float64,2},Xi::Float64,cM::Base.LinAlg.Cholesky{Float64,Array{Float64,2}},n::Int64,a0::Float64,b0::Float64)
    ssr = y'*x+b0; ssr = ssr[1];
    try
        ldetM = 2.0*logdet(cM);
        ll = -.5.*ldetM - ((n+a0)/2).*log(ssr);
        lpr = -log(sqrt(Xi).*(1+Xi));
        lr = ll+lpr;
    catch
        lr = -Inf; info("proposal was rejected because I+XDX was not positive-definite");
    end
end

function lmh_ratio(y::Array{Float64,2},x::Array{Float64,2},Xi::Float64,cM::Array{Float64,1},n::Int64,a0::Float64,b0::Float64)
    ssr = y'*x+b0; ssr = ssr[1];
    try
        ldetM = 2.0*sum(log(cM));
        ll = -.5.*ldetM - ((n+a0)/2).*log(ssr);
        lpr = -log(sqrt(Xi).*(1+Xi));
        lr = ll+lpr;
    catch
        lr = -Inf; info("proposal was rejected because I+XDX was not positive-definite");
    end
end


function gen_truncated_exp(mn::Array{Float64,2},trunc_point::Array{Float64,2})
    r = mn.*trunc_point;
    sml = abs.(r).<eps(1.0);
    x = zeros(max(length(mn),length(trunc_point)),1);
    tmp = zeros(max(length(mn),length(trunc_point)),1);

    if any(sml)
        tmp[sml] = expm1(-r[sml]).*rand(sum(sml),1);
    end
    tmp[.~(sml)] = (exp(-r[.~(sml)])-1.0).*rand(sum(~sml),1);

    sml = abs(tmp).<eps(1.0);
    if any(sml)
       x[sml] = -log1p(tmp[sml])./mn[sml];
    end
    x[.~(sml)] = -log(1.0+tmp[.~(sml)])./mn[.~(sml)];
end


function horseshoe(y::Array{Float64,2},X::Array{Float64,2},BURNIN::Int64,MCMC::Int64,scl::Float64,SAVE_SAMPLES::Bool,a0::Float64,b0::Float64,BetaTrue::Array{Float64,2},disp_int::Int64,plotting::Bool,
  corX::Bool,simtype::String,nkeep::Int64,ApproxXLX::Bool,is_sim::Bool,yt::Array{Float64,2},Xt::Array{Float64,2},delt::Float64,rhoX::Float64,mh_sigma::Bool,s_sigma::Float64)


if ApproxXLX
    simtype = "$(simtype)_approx_$(-log10(delt))";
end
if corX
   simtype = "$(simtype)_cor_$(rhoX)";
end

N=BURNIN+MCMC;
np0=size(X); n = np0[1]; p = np0[2];


Beta=ones(p,1);
lambda=ones(p,1);
tau=1.; sigma_sq=1.;

betaout=zeros(nkeep,MCMC);
lambdaout=zeros(nkeep,MCMC);
etaout = zeros(nkeep,MCMC);
tauout=zeros(MCMC,1);
xiout = zeros(MCMC+BURNIN,1);
sigmaSqout=zeros(MCMC,1);
l1out = zeros(MCMC+BURNIN,1);
pexpout = zeros(MCMC+BURNIN,1);
ACC = falses(MCMC+BURNIN,1);
ACC_s = falses(MCMC+BURNIN,1);


if issparse(X)
   I_n = sparse(I,n,n);
else
    I_n=Matrix(I,n,n);
end
l=ones(n,1);
u = zeros(p,1);

Xi = tau^(-2);
Eta = lambda.^(-2);

if ~is_sim
   Beta_hat = zeros(p,1);
end


t = 0.;
t1 = time();

for i=1:N
  if mod(i,disp_int)==0
      println("$(i)"); dt = time()-t1; t1=time();
      t = t+dt;
      println("$(disp_int) iterations in $(dt) seconds");
  end

  # update tau
  Eta = lambda.^(-2);
  prop_Xi = exp(rand(Normal(log(Xi),scl)));


  if ApproxXLX
      which_in = (lambda.^2.*max(1/Xi,1/prop_Xi)).>delt;
      id1 = find(which_in);
      rXLX = sum(which_in);

      lambda1 = lambda[id1]; X1 = X[:,id1];

      LX = broadcast(*,lambda1.^2,X1');
      LX1 = broadcast(*,lambda1,X1');

      if rXLX < n/2
          XX1 = X1'*X1;
      else
          XLX = LX1'*LX1;
      end
  else
      LX = broadcast(*,lambda.^2,X');
      LX1 = broadcast(*,lambda,X');
      XLX = LX1'*LX1;
      rXLX = n;
  end




  if rXLX<n/2 && ApproxXLX
      M_prop = [];
      XL0 = broadcast(*,X1,lambda1');
      s = factorize(XL0'*XL0);
      s0 = diag(Matrix(s));
      s_prop = s0.*(1./sqrt(prop_Xi));
      cM_prop = sqrt.(1+s_prop.^2);
      M = [];
      s_curr = s0.*(1./sqrt(Xi));
      cM = sqrt.(1+s_curr.^2);
      X1y = X1'*y;
      cXX_prop = factorize(prop_Xi.*diagm(lambda1.^(-2))+XX1);
      x_prop = y - X1*(cXX_prop\X1y);
      cXX = factorize(Xi.*diagm(lambda1.^(-2))+XX1);
      x = y - X1*(cXX\X1y);
  else
      M = I_n + (1./Xi).*XLX;
      M_prop = I_n + (1./prop_Xi).*XLX;
      cM = factorize(M); cM_prop = factorize(M_prop);
      x_prop = cM_prop\y;
      x = cM\y;
  end
  #
  #
  lrat_prop = lmh_ratio(y,x_prop,prop_Xi,cM_prop,n,a0,b0);
  lrat_curr = lmh_ratio(y,x,Xi,cM,n,a0,b0);
  log_acc_rat = (lrat_prop-lrat_curr)+(log(prop_Xi)-log(Xi));
  #
  #
  ACC[i] = (rand() < exp(log_acc_rat));
  if ACC[i]
     Xi = prop_Xi;
     M = M_prop;
     x = x_prop;
     if rXLX<n/2 && ApproxXLX
       cXX = cXX_prop;
     end
     cM = cM_prop;
  end
  tau = 1./sqrt(Xi);


  # update sigma_sq
  #
  #
  # metropolis option
  if mh_sigma
     ssr = y'*x;
     prop_k = exp(rand(Normal(log(1.0/sigma_sq),s_sigma)));


     l_prop = (n+a0-2.0)/2.0*log(prop_k)-prop_k*(ssr+b0)/2.0;
     l_curr = (n+a0-2.0)/2.0*log(1./sigma_sq)-(1.0/sigma_sq)*(ssr+b0)/2;
     l_ar = (l_prop-l_curr) + (log(prop_k)-log(1.0/sigma_sq));
     l_ar = l_ar[1];
     acc_s = rand()<exp(l_ar);
     if acc_s
         sigma_sq = 1/prop_k;
     end
     ACC_s[i] = acc_s;
  else
     ssr = y'*x;
     sigma_sq = 1/rand(Gamma((n+a0)/2.0,2.0/(ssr+b0)));
  end




  #
  #
  U = (1./Xi).*LX;
  #
  # step 1
  [u[j,1] = rand(Normal(0.0,tau.*lambda[j])) for j in 1:p];
  v=X*u+rand(Normal(0.0,1.0),n,1);
  #
  # step 2
  #
  if ApproxXLX && rXLX<n/2 # woodbury
    tmp = cXX\(X1'*(y./sqrt(sigma_sq)-v));
    v_star = (y./sqrt(sigma_sq)-v)-X1*tmp;
    Beta = sqrt(sigma_sq)*u;
    Beta[id1] = Beta[id1] + sqrt(sigma_sq)*(U*v_star);
  elseif ApproxXLX && rXLX>=n/2 # direct solve
    v_star=cM\(y./sqrt(sigma_sq)-v);
    Beta = sqrt(sigma_sq)*u;
    Beta[id1] = Beta[id1] + sqrt(sigma_sq)*(U*v_star);
  else
    v_star=cM\(y./sqrt(sigma_sq)-v);
    Beta=sqrt(sigma_sq)*(u+U*v_star);
  end

  #
  if ~is_sim
     if i < floor(BURNIN/2)
         Beta_hat = (i-1.0)./i.*Beta_hat + 1./i.*Beta;
     elseif i==floor(BURNIN/2)
         Beta_hat = Beta;
         se = Beta.^2;
     else
         idrem = floor(BURNIN/2);
         Beta_hat = (i-idrem-1.0)./(i-idrem).*Beta_hat + 1./(i-idrem).*Beta;
         se = (i-idrem-1.0)./(i-idrem).*se + 1./(i-idrem).*Beta.^2;
     end
  end
  #
  [u[j,1] = rand(Uniform(0.0, 1./(Eta[j,1]+1))) for j in 1:p];
  gamma_rate = (Beta.^2) .* Xi ./ (2.*sigma_sq);
  Eta_old = Eta;
  Eta = gen_truncated_exp(gamma_rate, (1-u)./u);
  if length(Eta)<p
    println("mystery shortening of Eta occurred");
    Eta = Eta_old;
  end

  if any(Eta.<=0)
     println("$(sum(Eta.<=0)) Eta underflowed, replacing = machine epsilon");
     Eta[Eta.<=0] = eps(1.0);
  end
  lambda = 1./sqrt.(Eta);
  #

  #println("$(i)");
  if mod(i,disp_int) == 0
    i
    mean(ACC_s[1:i])
    mean(ACC[1:i])
    println("rank of XLX $(rXLX)");
    if is_sim
      Beta[1:20]
    else
      Bet_ord = sort(Beta_hat,rev=true);
      Bet_ord[1:20]
    end
    if plotting
#      boxplot(betaout[1:25,1:i]')
    end


  end
  #         if plotting && i > BURNIN
  #             id_upper = i-BURNIN-1;
  #             figure(1);subplot(2,2,1);plot(log(xiout(50:id_upper)),'.');title('log(\xi)');
  #             subplot(2,2,2);plot(log(1./sigmaSqout(50:id_upper)),'.');title('log(\sigma^{-2})');
  #             subplot(2,2,3);plot(betaout(1,50:id_upper),'.');title('\beta_1');
  #             subplot(2,2,4);plot(etaout(1,50:id_upper),'.');title('\eta_1');
  #             drawnow;
  #             figure(2); subplot(1,2,1);plot(l1out(1:id_upper),'.');title('L1');ylim([0 1]);
  #             subplot(1,2,2);plot(pexpout(1:id_upper),'.');title('pexpl'); ylim([0 1]);
  #             drawnow;
  #             figure(3);
  #             for j0=1:25
  #                 subplot(5,5,j0); plot(betaout(j0,50:id_upper),'.');
  #             end
  #             drawnow;
  #             figure(5); hist(log10((1./Xi).*lambda.^2),100); title('\xi^{-1} \lambda^2_j');drawnow;
  #
  #             if ~is_sim
  #                 figure(4);
  #                 for j0=1:6
  #                     subplot(2,3,j0); plot(betaout(j0,50:id_upper),'.');
  #                 end
  #                 drawnow;
  #             end
  #         end
  #     end
  #
  #     if is_sim
  #         per_expl = 1-sqrt(sum((Beta-BetaTrue).^2))./sqrt(sum(BetaTrue.^2));
  #         L1_loss = 1-sum(abs(Beta-BetaTrue))./sum(abs(BetaTrue));
  #     else
  #         res1 = yt-Xt(:,id1)*Beta_hat(id1);
  #         per_expl = 1-sqrt(sum(res1.^2))./sqrt(sum(yt.^2));
  #         L1_loss = 1-sum(abs(res1))./sum(abs(yt));
  #     end
  #
    if i==BURNIN || (i==1 && BURNIN==0)
  #         if ~is_sim
  #             [~,big_id] = sort(abs(Beta_hat),'descend');
  #             big_id = big_id(1:nkeep/2);
  #             other_id = setdiff((1:p)',big_id);
  #             other_id = other_id(1:nkeep/2);
  #             keep_id = [big_id;other_id];
  #         else
      keep_id = range(1,nkeep);
  #         end
    end
  #
    if i > BURNIN
      keep_id = range(1,nkeep);
      betaout[:,i-BURNIN] = Beta[keep_id];
  #         lambdaout(:,(i-BURNIN)/thin) = lambda(keep_id);
  #         etaout(:,(i-BURNIN)/thin) = Eta(keep_id);
  #         tauout((i-BURNIN)/thin)=tau;
  #         xiout((i-BURNIN)/thin) = Xi;
  #         sigmaSqout((i-BURNIN)/thin)=sigma_sq;
  #         l1out((i-BURNIN)/thin) = L1_loss;
  #         pexpout((i-BURNIN)/thin) = per_expl;
    end
end

bo = DataFrame(betaout');
if SAVE_SAMPLES
  CSV.write("/home/james/Documents/GitHub/horseshoe_jo/Outputs/julia_test.csv",DataFrame(betaout'));
end

return bo

end

# pMean=mean(betaout,2);
# pMedian=median(betaout,2);
# pSigma=mean(sigmaSqout);
# pLambda=mean(lambdaout,2);
# dt=toc; t = t+dt;
#
# ci_lo = quantile(betaout(:,5001:end)',.025,1);
# ci_hi = quantile(betaout(:,5001:end)',.975,1);
#
# if is_sim
#     coverage = mean(BetaTrue(1:nkeep)>ci_lo' & BetaTrue(1:nkeep)<ci_hi');
#     BetaHat = mean(betaout(:,5001:end),2)';
#     Beta_hat = BetaHat;
#
#     mse = mean((BetaTrue(1:nkeep)-BetaHat').^2);
#     se = std(betaout(:,5001:end),[],2);
#
#     disp(['coverage ' num2str(100*coverage)]);
#     disp(['mse ' num2str(mse)]);
#     keep_id = [];
# else
#    coverage = []; mse = [];
#    BetaHat = Beta_hat;
# end
#
#
# disp([num2str(t) ' seconds elapsed']);
# if SAVE_SAMPLES
#     save(strcat('Outputs/post_reg_horse_',simtype,'_',num2str(n),'_',num2str(p),'.mat'),'betaout','lambdaout','etaout','tauout','xiout','sigmaSqout','l1out','pexpout','t',...
#         'ci_hi','ci_lo','coverage','BetaHat','mse','se','BetaTrue','keep_id');
# end
#
#
# figure(1);subplot(2,2,1);plot(log(xiout(50:(i-BURNIN-1))),'.');title('log(\xi)');
# subplot(2,2,2);plot(log(1./sigmaSqout(50:(i-BURNIN-1))),'.');title('log(\sigma^{-2})');
# subplot(2,2,3);plot(betaout(1,50:(i-BURNIN-1)),'.');title('\beta_1');
# subplot(2,2,4);plot(etaout(1,50:(i-BURNIN-1)),'.');title('\eta_1');
# drawnow;
# figure(2); subplot(1,2,1);plot(l1out(1:(i-BURNIN-1)),'.');title('L1');ylim([0 1]);
# subplot(1,2,2);plot(pexpout(1:(i-BURNIN-1)),'.');title('pexpl'); ylim([0 1]);
# drawnow;
# figure(3);
# for j0=1:25
#     subplot(5,5,j0); plot(betaout(j0,50:(i-BURNIN-1)),'.'); title(strcat('\beta_{',num2str(j0),'}'));
# end
# drawnow;

#end
