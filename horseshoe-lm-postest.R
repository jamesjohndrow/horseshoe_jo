rm(list=ls(all=T))
setwd('~/Dropbox/Projects/Horseshoe/Code/code_matlab')
require('R.matlab')
require('reshape2')
require('ggplot2')
require('coda')
require('xtable')
require('condmixt')
require('evir')

acf.mcmc <- function(x) {
  return(autocorr(as.mcmc(x),lags = c(1,seq(from=5,to=100,by=5))))
}

es.mcmc <- function(x) {
  return(effectiveSize(as.mcmc(x)))
}


process.acfdat <- function(L,vname,algo,n.sav,n,p) {
  df = data.frame(L)
  names(df) = paste(vname,seq(n.sav),sep='')
  df$lag = c(1,seq(from=5,to=100,by=5))
  df$n = n
  df$p = p
  df$algo = algo
  df = melt(df,id=c('lag','n','p','algo'))
  return(df)
}

process.rawdat <- function(L,vname,algo,n.sav,n,p,nmc) {
  df = data.frame(t(L))
  names(df) = paste(vname,seq(n.sav),sep='')
  df$iter = seq(nmc)
  df$n = n
  df$p = p
  df$algo = algo
  df = melt(df,id=c('iter','n','p','algo'))
  return(df)
}
  
process.acfmat <- function(M,vname,algo,NP) {
  df = data.frame(M)
  names(df) = c(1,seq(from=5,to=100,by=5))
  df$n = NP[,1]
  df$p = NP[,2]
  df$algo = algo
  df$variable = vname
  df = melt(df,id=c('n','p','algo','variable'))
  names(df)[5] = 'lag'
  df = df[,c(5,1:4,6)]
  return(df)
}

process.essmat <- function(M,vname,algo,n.sav,NP) {
  df = data.frame(M)
  names(df) = paste(vname,seq(n.sav),sep='')
  df$n = NP[,1]
  df$p = NP[,2]
  df$algo = algo
  df = melt(df,id=c('n','p','algo'))
  return(df)
}

cov.bet <- function(BET,BetaTrue) {
  ci = t(apply(BET,1,quantile,probs=c(.025,.975)))
  is.in = (BetaTrue>ci[,1]) & (BetaTrue<ci[,2])
  return(mean(is.in))
}


nbet.sav <- 50
shift.eps <- 0
nmc <- 40000
burn <- 1000
burn_fact <- (nmc-burn)/nmc
ns <- c(500,1000)
ps <- c(1000,5000,10000,20000)
nn <- length(ns)
np <- length(ps)
ncs <- 2

BetaTrue <- matrix(0,nbet.sav,1)
BetaTrue[1:5] = 4
BetaTrue[6:15] = 2^(-(seq(from=0,to=4.5,by=.5)));

ES.xi <- matrix(0,nn*np,4)
ESS.xi <- matrix(0,nn*np,4)
ESS.Sigma2 <- matrix(0,nn*np,4)

ESS.bet.block <- matrix(0,nn*np,nbet.sav)
ESS.bet.slice <- matrix(0,nn*np,nbet.sav)
ESS.eta.block <- matrix(0,nn*np,nbet.sav)
ESS.eta.slice <- matrix(0,nn*np,nbet.sav)
ACF.xi.block <- matrix(0,nn*np,21)
ACF.xi.slice <- matrix(0,nn*np,21)
ACF.bet.block <- list()
ACF.bet.slice <- list()
ACF.eta.block <- list()
ACF.eta.slice <- list()

BET.block <- list()
BET.slice <- list()

XI.block <- matrix(0,nn*np,40000)
XI.slice <- matrix(0,nn*np,40000)

NP <- matrix(0,nn*np,2)
T.all <- matrix(0,nn*np,2)

ctr <- 0
for (n in ns) {
  for (p in ps) {
    print(c(n,p))
    ctr <- ctr + 1
    
    NP[ctr,] <- c(n,p)
    
    ACF.bet.block[[ctr]] <- matrix(0,nbet.sav,21)
    ACF.eta.block[[ctr]] <- matrix(0,nbet.sav,21)
    ACF.bet.slice[[ctr]] <- matrix(0,nbet.sav,21)
    ACF.eta.slice[[ctr]] <- matrix(0,nbet.sav,21)
    
    dat <- readMat(paste('Outputs/post_reg_horse_',n,'_',p,'.mat',sep=''))
    BET.block[[ctr]] <- dat$betaout; XI.block[ctr,] <- dat$xiout
    T.all[ctr,1] <- dat$t
    es <- effectiveSize(dat$xiout[burn:nmc]+shift.eps)
    ES.xi[ctr,1:2] <- c(n,p); ESS.xi[ctr,1:2] <- c(n,p)
    ES.xi[ctr,3] <- es/(nmc-burn); ESS.xi[ctr,3] <- es/(dat$t*burn_fact)
    
    ESS.Sigma2[ctr,1:2] <- c(n,p); ESS.Sigma2[ctr,3] <- effectiveSize(dat$sigmaSqout[burn:nmc]+shift.eps)/(dat$t*burn_fact)
    
    ESS.bet.block[ctr,] <- apply(dat$betaout[,burn:nmc]+shift.eps,1,effectiveSize)/(dat$t*burn_fact)
    ESS.eta.block[ctr,] <- apply(dat$etaout[,burn:nmc]+shift.eps,1,effectiveSize)/(dat$t*burn_fact)
    
    ACF.xi.block[ctr,] <- acf.mcmc(dat$xiout[burn:nmc])
    
    ACF.bet.block[[ctr]] <- apply(dat$betaout[,burn:nmc],1,acf.mcmc)
    ACF.eta.block[[ctr]] <- apply(dat$etaout[,burn:nmc],1,acf.mcmc)
    
    dat0 <- readMat(paste('Outputs/post_reg_horse_ab_',n,'_',p,'.mat',sep=''))
    BET.slice[[ctr]] <- dat0$betaout; XI.slice[ctr,] <- dat0$xiout
    
    T.all[ctr,2] <- dat0$t
    es <- effectiveSize(dat0$xiout[burn:nmc]+shift.eps)
    ES.xi[ctr,4] <- es/(nmc-burn); ESS.xi[ctr,4] <- es/(dat0$t*burn_fact)
    ESS.Sigma2[ctr,4] <- effectiveSize(dat$sigmaSqout[burn:nmc]+shift.eps)/(dat0$t*burn_fact)
    
    ESS.bet.slice[ctr,] <- apply(dat0$betaout[,burn:nmc]+shift.eps,1,effectiveSize)/(dat0$t*burn_fact)
    ESS.eta.slice[ctr,] <- apply(dat0$etaout[,burn:nmc]+shift.eps,1,effectiveSize)/(dat0$t*burn_fact)
    
    ACF.xi.slice[ctr,] <- acf.mcmc(dat0$xiout[burn:nmc])
    
    ACF.bet.slice[[ctr]] <- apply(dat0$betaout[,burn:nmc],1,acf.mcmc)
    ACF.eta.slice[[ctr]] <- apply(dat0$etaout[,burn:nmc],1,acf.mcmc)
    
  }
}


# stuff that focuses on xi
ES.xi <- data.frame(ES.xi)
names(ES.xi) <- c('n','p','mh','slice')
rownames(ES.xi) <- NULL
lm1 <- lm(log(1/mh)~log(p)+log(n),data=ES.xi)
lm2 <- lm(log(1/slice)~log(p)+log(n),data=ES.xi)
summary(lm1)
summary(lm2)

ESS.xi <- cbind(ESS.xi,ESS.xi[,3]/ESS.xi[,4])
colnames(ESS.xi) <- c('n','p','block','slice','efficiency ratio')

print.xtable(xtable(ES.xi,digits=0),file='Figures/ess_xi_lm.tex',include.rownames = F)
print.xtable(xtable(ESS.xi,digits=c(0,0,0,4,4,2),caption='Effective samples per second for $\\xi$',label = 'tab:ess_xi'),
             file='Figures/essec_xi_lm.tex',include.rownames = F,caption.placement = 'top')


dat <- readMat('Outputs/post_reg_horse_1000_20000.mat')
dat0 <- readMat('Outputs/post_reg_horse_ab_1000_20000.mat')

XI <- cbind(dat$xiout,dat0$xiout)
XI <- data.frame(XI)
names(XI) <- c('block','slice')
XI$iter <- seq(nrow(XI))
XI <- melt(XI,id='iter')
names(XI)[2] <- 'sampler'

s.xi <- sort(dat$xiout[burn:nmc],decreasing = T)
hill(s.xi[1:1000])
hillest(dat$xiout[burn:nmc],500)


png(filename='Figures/xi_trace.png',width=900,height=600)
ggplot(XI,aes(x=iter,y=value,col=sampler)) + geom_point(size=.5,alpha=.3) + theme(text=element_text(size=32)) +
labs(y=expression(xi))
dev.off()

XI2 <- XI
XI2$value <- log(XI2$value)
XI2$variable <- 'log(xi)'
XI$variable <- 'xi'
XI <- rbind(XI,XI2)
png(filename = 'Figures/xi_dens.png',width=900,height=600)
ggplot(XI,aes(x=value,col=sampler,fill=sampler)) + geom_density(alpha=.25) + 
  theme(text=element_text(size=32),axis.text.x = element_text(angle=90)) + facet_wrap(~variable,scales='free')
dev.off()

bet <- data.frame(cbind(dat$betaout[10,1000:40000],dat0$betaout[10,1000:40000]))
names(bet) <- c('block','slice')
bet$iter <- seq(nrow(bet))
bet <- melt(bet,id='iter')
names(bet)[2] <- 'sampler'

png(filename = 'Figures/trace_bet10.png',width=900,height=600)
ggplot(bet,aes(x=iter,y=value)) + geom_point(size=.2) + facet_wrap(~sampler,scales = 'free') + theme(text=element_text(size=32)) +
  labs(y=expression(beta[10]))
dev.off()


# general summaries for other parameters
for (j in 1:(nn*np)) {
  n <- NP[j,1]; p <- NP[j,2]
  print(c(n,p))
  if (j==1) {
    df.raw <- process.rawdat(BET.block[[j]],'bet','block',nbet.sav,n,p,nmc)
    df.acf <- process.acfdat(ACF.bet.block[[j]],'bet','block',nbet.sav,n,p)
  } else {
    df <- process.rawdat(BET.block[[j]],'bet','block',nbet.sav,n,p,nmc)
    df.raw <- rbind(df.raw,df)
    
    df <- process.acfdat(ACF.bet.block[[j]],'bet','block',nbet.sav,n,p)
    df.acf <- rbind(df.acf,df)
  }  
    df <- process.acfdat(ACF.bet.slice[[j]],'bet','slice',nbet.sav,n,p)
    df.acf <- rbind(df.acf,df)
    
    df <- process.rawdat(BET.slice[[j]],'bet','slice',nbet.sav,n,p,nmc)
    df.raw <- rbind(df.raw,df) 
    
    df <- process.acfdat(ACF.eta.block[[j]],'eta','block',nbet.sav,n,p)
    df.acf <- rbind(df.acf,df)
    
    df <- process.acfdat(ACF.eta.slice[[j]],'eta','slice',nbet.sav,n,p)
    df.acf <- rbind(df.acf,df)
}

png('Figures/acf_box.png',width=1200,height=800)
ggplot(df.acf[df$n==1000,],aes(x=factor(lag),y=value)) + 
  geom_boxplot(outlier.size = .4) + facet_wrap(~algo+factor(p),ncol=4) + 
  theme(axis.text.x = element_text(angle=90),text=element_text(size=24)) + labs(x='lag',y='autocorrelation')
dev.off()

# ESS for variables other than xi
df.ess <- process.essmat(ESS.bet.block,'bet','block',nbet.sav,NP)
df <- process.essmat(ESS.bet.slice,'bet','slice',nbet.sav,NP)
df.ess <- rbind(df.ess,df)
df <- process.essmat(ESS.eta.block,'eta','block',nbet.sav,NP)
df.ess <- rbind(df.ess,df)
df <- process.essmat(ESS.eta.slice,'eta','slice',nbet.sav,NP)
df.ess <- rbind(df.ess,df)

png('Figures/ess_box.png',width=500,height=600)
ggplot(df.ess,aes(x=factor(p),y=value)) + geom_boxplot() + facet_wrap(~algo+n,ncol=2) + scale_y_log10() +
  theme(text=element_text(size=24),axis.text.x = element_text(angle=90)) + labs(x='p',y=expression(T[e]/t))
dev.off()

ess.Sigma2 <- ESS.Sigma2
colnames(ess.Sigma2) <- c('n','p','block','slice')
print.xtable(xtable(ess.Sigma2,digits = c(0,0,0,3,3),caption='Effective sample
             size per second for $\\sigma^2$',label='tab:ess_sigma2'),include.rownames = F,file='Figures/sigma2ess.tex',caption.placement='top')


df.acf.xi <- process.acfmat(ACF.xi.block,'xi','block',NP)
df <- process.acfmat(ACF.xi.slice,'xi','slice',NP)
df.acf.xi <- rbind(df.acf.xi,df)


png('Figures/acf_xi.png',width=900,height=900)
ggplot(df.acf.xi,aes(x=factor(lag),y=value,shape=factor(p))) + geom_point() + facet_wrap(~algo+n,ncol=2) +
  theme(axis.text.x = element_text(angle=90),text=element_text(size=24)) + labs(y="autocorrelation",x="lag")
dev.off()


# estimating complexity using all runs
ps <- c(seq(from=1000,to=10000,by=1000),15000,20000)
ns <- c(500,1000)

nn <- length(ns)
np <- length(ps)
ncs <- 2
shift.eps <- 0


ES.xi <- matrix(0,nn*np,4)
ESS.xi <- matrix(0,nn*np,4)
AC1.xi <- matrix(0,nn*np,4)


ctr <- 0
for (n in ns) {
  for (p in ps) {
    print(c(n,p))
    ctr <- ctr + 1
    
    dat <- readMat(paste('Outputs/post_reg_horse_',n,'_',p,'.mat',sep=''))
    es <- effectiveSize(dat$xiout[burn:nmc]+shift.eps)
    ES.xi[ctr,1:2] <- c(n,p); ESS.xi[ctr,1:2] <- c(n,p)
    ES.xi[ctr,3] <- es/(nmc-burn); ESS.xi[ctr,3] <- es/(dat$t*burn_fact)
    AC1.xi[ctr,1:2] <- c(n,p)
    AC1.xi[ctr,3] <- autocorr(as.mcmc(dat$xiout[burn:nmc]),lags=1)
    
    dat0 <- readMat(paste('Outputs/post_reg_horse_ab_',n,'_',p,'.mat',sep=''))
    es <- effectiveSize(dat0$xiout[burn:nmc]+shift.eps)
    ES.xi[ctr,4] <- es/(nmc-burn); ESS.xi[ctr,4] <- es/(dat0$t*burn_fact)
    AC1.xi[ctr,4] <- autocorr(as.mcmc(dat0$xiout),lags=1)
    
  }
}

ES.df <- data.frame(ES.xi)
names(ES.df) <- c('n','p','block','slice')
lm.es.block <- lm(log(1/block)~log(n)+log(p),data=ES.df)
lm.es.slice <- lm(log(1/slice)~log(n)+log(p),data=ES.df)

summary(lm.es.block)
summary(lm.es.slice)

ES.df <- melt(ES.df,id=c('n','p'))
names(ES.df) <- c('n','p','sampler','Te')

png(filename = 'Figures/ES_np.png',width=800,height=500)
ggplot(ES.df,aes(x=log(p),y=log(1/Te),shape=factor(n))) + geom_point(size=3) + facet_wrap(~sampler) + 
  theme(text=element_text(size=24)) + labs(y=expression(log(1/T[e])))
dev.off()

ggplot(ES.df,aes(x=log(p/n),y=log(1/Te))) + geom_point(size=3) + facet_wrap(~variable) + theme(text=element_text(size=24))


dat1 <- readMat('Outputs/post_reg_horse_ab_trunc_500_10000.mat')
xi.trunc <- dat1$xiout
xi.trunc <- data.frame(xi.trunc)
names(xi.trunc) <- 'xi'
xi.trunc$iter <- seq(nrow(xi.trunc))

png(filename = 'Figures/xi_trace_trunc.png',width=800,height=500)
ggplot(xi.trunc,aes(x=iter,y=log(xi))) + geom_point(size=.5) + theme(text=element_text(size=24))
dev.off()






