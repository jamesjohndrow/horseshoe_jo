rm(list=ls(all=T))
setwd('/data1/GitHub/horseshoe/Code/')

if (!require(pacman)) {install.packages('pacman')}
pacman::p_load(R.matlab,coda,h5,ggplot2,gridExtra,ggpubr,xtable,tidyr,dplyr,mcmcse,readr,parallel,ARTP2)

f.es <- function(x){ess(x,method='obm',size='cuberoot')}
f.ac <- function(x){
  ac <- autocorr(as.mcmc(x),lags = seq(100))
}

ct.common <- function(t,x,y) {
  c1 <- sum(abs(x)>t)
  c2 <- sum(abs(y)>t)
  c3 <- sum(abs(x)>t & abs(y)>t)
  return(c(c1,c2,c3))
}



ab.vs <- 'ab_mrg'

burn <- 0
nmc <- 2000
#######
# kl, w
#######

dat <- readMat('Outputs/post_reg_horse_F_250_1000.mat')

df <- data.frame(KL=dat$KL[100:nmc],W2 = dat$W2[100:nmc],KL2 = dat$KL2[100:nmc],iter=seq(from=100,to=nmc))
df <- df %>% gather(variable,value,-iter)

png('Figures/KL_W2.png',width = 600,height = 400)
ggplot(df,aes(x=value)) + geom_histogram() + facet_wrap(~variable,scales = 'free') 
dev.off()

df.W <- data.frame(W2 = dat$W2,iter = seq(nmc))

png('Figures/W2.png',width = 600,height = 400)
ggplot(df.W,aes(x=W2)) + geom_histogram() + theme(text=element_text(size=24)) + xlab(expression(W[2]))
dev.off()

burn <- 0
nmc <- 25000

########
# maize
########

dat <- readMat('Outputs/post_reg_horse_maize_approx_4_2267_98385.mat')
df.maize <- read.bed('Data/gwas/maize/imputemaize.bed', 'Data/gwas/maize/imputemaize.bim', 'Data/gwas/maize/imputemaize.fam',encode012 = TRUE)
load('Outputs/lasso_maize.RData')

bet.horse <- dat$BetaHat
lasso.in <- which(bets!=0)

df.bets <- data.frame(beta.horse = bet.horse,beta.lasso=bets,id.num=seq(length(bet.horse)))
threshes <- seq(from=0.1,to=0.0005,by=-0.001)
threshes <- as.matrix(threshes,length(threshes),1)
cts <- apply(threshes,1,ct.common,x=df.bets$beta.horse,y=df.bets$beta.lasso)
cts <- t(cts)
df.cts <- data.frame(cts=cts)
names(df.cts) <- c('horse','lasso','intersect')
df.cts$threshold <- threshes
df.cts <- df.cts %>% gather(method,value,-threshold)


df.cts$method <- factor(df.cts$method)

png('Figures/lasso_horse_compare.png',width=600,height=400)
ggplot(df.cts,aes(x=threshold,y=value+1,col=method)) + geom_line(size=2) + scale_y_log10(breaks=c(5,10,50,100,500,1000)) +
  theme(text=element_text(size=32)) + ylab('1 + # exceeding T')
dev.off()

# ggplot(df.bets,aes(x=beta.horse,y=beta.lasso)) + geom_point(size=.3,alpha=.3)
# ggplot(df.bets[abs(df.bets$beta.lasso)>.01,],aes(x=beta.horse,y=beta.lasso)) + geom_point(size=.3,alpha=.3)
# ggplot(df.bets[abs(df.bets$beta.lasso)>.01,],aes(x=beta.horse)) + stat_ecdf()
# ggplot(df.bets[abs(df.bets$beta.lasso)>.01,],aes(x=beta.horse)) + stat_ecdf()

# X <- as.matrix(df.maize)
# X.sub <- X[,dat$keep.id[1:10]]
# Cr <- cor(X.sub,X[,-dat$keep.id[1:10]])
# s.Cr <- sort(c(Cr),decreasing = T)
# which1 <- apply(as.matrix(Cr),1,function(x){which(x==1)})

bet <- t(dat$betaout)
lambda <- t(dat$lambdaout)
tm <- c(dat$t)

es <- c(apply(bet[burn:nmc,],2,f.es),apply(lambda[burn:nmc,],2,f.es),f.es(log(dat$xiout[burn:nmc])),f.es(-2*log(dat$sigmaSqout[burn:nmc])))
df.es <- data.frame(es=c(es),var=c(paste('beta',seq(200),sep=''),paste('lambda',seq(200),sep=''),'log_xi','-2log(sigma2)'),ess=c(es)/tm)

c(min(df.es$es),min(df.es$ess))
c(median(df.es$es),median(df.es$ess))

p1 <- ggplot(df.es,aes(x=es)) + geom_histogram() + scale_x_log10(limits=c(500,25000),breaks=c(500,1000,5000,10000,20000)) + 
  theme(text=element_text(size=24),axis.text.x=element_text(angle=90)) + xlab('') + ggtitle(expression(n[e]))
p2 <- ggplot(df.es,aes(x=ess)) + geom_histogram() + scale_x_log10(limits=c(0.05,1.1),breaks=c(0.05,0.1,0.5,1)) + 
  theme(text=element_text(size=24),axis.text.x=element_text(angle=90)) + xlab('') + ggtitle(expression(n[e]/t))

png('Figures/ess_maize.png',width=600,height=300)
ggarrange(p1,p2,ncol=2)
dev.off()

keep.id <- dat$keep.id
df.bet <- data.frame(bet)
mn.bet <- apply(bet,2,mean)
rnks.bet <- rank(-abs(mn.bet))
names(df.bet) <- keep.id
df.bet$iter <- seq(nmc)
df.bet <- df.bet %>% gather(variable,value,-iter)
df.bet$variable <- as.numeric(df.bet$variable)

df.bet <- df.bet %>% mutate(is.in=abs(value)>1e-3)
df.mip <- df.bet %>% group_by(variable) %>% summarise(mip=mean(is.in))

png('Figures/MIP-maize.png',width=900,height=600)
ggplot(df.mip,aes(x=mip)) + geom_histogram() + theme(text=element_text(size=32)) + xlab(expression(P(abs(beta[j])>.001)))
dev.off()

bet.size <- apply(bet,1,function(x){sum(abs(x)>1e-3)})
df.size <- data.frame(iter=seq(nmc),size=bet.size)

png('Figures/ModSize.png',width=900,height=600)
ggplot(df.size,aes(x=size)) + geom_histogram(bins = 20) + theme(text=element_text(size=32)) + xlab(expression(abs(j : abs(beta[j])>.001)))
dev.off()

df.keep <- data.frame(id.num=keep.id)
rnks.df <- data.frame(rnk=rnks.bet,id.num=keep.id)
df.bet <- df.bet %>% inner_join(rnks.df,by=c("variable" = "id.num"))

#df.bet <- df.bet %>% inner_join(df.bets,by=c("variable" = "id.num"))
df.bets.keep <- rnks.df %>% inner_join(df.keep,by=c("id.num" = "id.num"))
df.bets.keep <- df.bets.keep %>% inner_join(df.bets,by=c("id.num" = "id.num"))

png('Figures/beta_marginals_maize.png',width=900,height=600)
ggplot(df.bet[df.bet$rnk<10,],aes(x=value)) + geom_density(size=1.2) + facet_wrap(~rnk,scales='free') + 
  theme(text=element_text(size=32),axis.text.x=element_text(angle=90)) + xlab(expression(beta)) +
  geom_vline(data=df.bets.keep[df.bets.keep$rnk<10,],aes(xintercept=beta.horse,colour='red')) + 
  geom_vline(data=df.bets.keep[df.bets.keep$rnk<10,],aes(xintercept=beta.lasso,colour='blue')) + 
  theme(legend.position = 'none')
dev.off()

df.bet.wide <- df.bet %>% filter(rnk>=10 & rnk<20) %>% select(-rnk) %>% spread(variable,value)

names(df.bet.wide)[2:ncol(df.bet.wide)] <- paste('X',names(df.bet.wide)[2:ncol(df.bet.wide)],sep='')

png('Figures/biplot_maize_example.png',width=900,height=600)
ggplot(df.bet.wide,aes(x=X19570,y=X32720)) +   stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  ylim(-.08,.02) + xlim(-.08,.02) + theme(text=element_text(size=24))
dev.off()

#ggplot(df.bet.wide,aes(x=bet72059,y=bet94631)) +   stat_density_2d(aes(fill = ..level..), geom = "polygon") 
#######
# the old application, currently not used in paper
#######

burn <- 5000
nmc <- 20000

load('Outputs/affymetrix_info.RData')
#resps <- c('aco','bactin','gluc10','gpat','pepck','scd1','srebp1','wt10wk')
resps <- c('gluc10','wt10wk')

resp <- 'scd1'
ctr <- 0
for (resp in resps) {
  ctr <- ctr + 1
  print(resp)
  dat <- readMat(paste('Outputs/post_reg_horse_genes_',resp,'_60_45265.mat',sep=''))
  bet <- t(dat$betaout)
  t <- dat$t
  lambda <- t(dat$lambdaout)
  
  es <- c(apply(bet[burn:nmc,],2,f.es),apply(lambda[burn:nmc,],2,f.es),f.es(log(dat$xiout[burn:nmc])),f.es(-2*log(dat$sigmaSqout[burn:nmc])))
  df.es <- data.frame(es=es,var=c(paste('beta',seq(100),sep=''),paste('lambda',seq(100),sep=''),'log_xi','-2log(sigma2)'),ess=es/dat$t)
  df.es$phenotype <- resp
  
  if (ctr==1) {
    ES <- df.es
  } else {
    ES <- rbind(ES,df.es)
  }
  
  ids <- dat$ids[1:2]
  aff_id <- id[ids]
  aff_id <- data.frame(ID=aff_id)
  genedat <- aff %>% inner_join(aff_id,by='ID')
  genedat <- aff[ids,]
  names(genedat) <- gsub(' ','_',names(genedat))
  
  gnames <- genedat$ID

  bet.big <- bet[,1:2]
  bet.big <- data.frame(bet.big)
  names(bet.big) <- gnames
  bet.big$iter <- seq(nmc)
  bet.big <- bet.big %>% gather(gene,value,-iter)
  bet.big$phenotype <- resp
  
  if (ctr==1) {
    B <- bet.big
  } else {
    B <- rbind(B,bet.big)
  }
}

#genedat <- read_tsv('Data/GPL340-26774-affyinfo-short.csv',col_names = TRUE,col_types = cols(.default="c"))
#genedat <- genedat[ids,]

p1 <- ggplot(ES,aes(x=ess)) + geom_histogram() + scale_x_log10() + ggtitle(expression(n[e]/t)) + facet_grid(~phenotype) + 
  theme(text=element_text(size=24)) + xlab('')
p2 <- ggplot(ES,aes(x=es)) + geom_histogram() + scale_x_log10() + ggtitle(expression(n[e])) + facet_grid(~phenotype) + 
  theme(text=element_text(size=24)) + xlab('')

#p1 <- ggplot(ES,aes(x=phenotype,y=ess)) + geom_boxplot() + scale_y_log10() + ggtitle(expression(n[e]/t))
#p2 <- ggplot(ES,aes(x=phenotype,y=es)) + geom_boxplot() + scale_y_log10() + ggtitle(expression(n[e]))

png('Figures/es_ess_gene.png',width=900,height=600)
ggarrange(p1,p2,nrow=2)
dev.off()

png('Figures/beta_marginals_gene.png',width=900,height=600)
ggplot(B[B$iter>burn,],aes(x=value)) + geom_histogram() + facet_wrap(~phenotype+gene) +
  theme(text=element_text(size=24)) + xlab(expression(beta))
dev.off()

print(min(ES$es))
print(min(ES$ess))
print(median(ES$es))
print(median(ES$ess))

####################
# the 2000-20000 examples
####################

ab <- readMat(paste('Outputs/post_reg_horse_',ab.vs,'_F_2000_20000.mat',sep=''))
ex <- readMat('Outputs/post_reg_horse_F_2000_20000.mat')
ap <- readMat('Outputs/post_reg_horse_F_approx_4_2000_20000.mat')

ess.ab <- c(apply(ab$betaout[,burn:nmc],1,f.es),apply(ab$etaout,1,f.es),f.es(log(ab$xiout)),f.es(log(1/ab$sigmaSqout)))/c(ab$t)
ess.ex <- c(apply(ex$betaout[,burn:nmc],1,f.es),apply(ex$etaout,1,f.es),f.es(log(ex$xiout)),f.es(log(1/ex$sigmaSqout)))/(200*119.7)
ess.ap <- c(apply(ap$betaout[,burn:nmc],1,f.es),apply(ap$etaout,1,f.es),f.es(log(ap$xiout)),f.es(log(1/ap$sigmaSqout)))/c(ap$t)

nm.es <- c(paste('bet',seq(100),sep=''),paste('eta',seq(100),sep=''),'xi','sigma2')
nm.es <- rep(nm.es,3)

df.es <- data.frame(var=nm.es,es=c(ess.ab,ess.ex,ess.ap),algo=rep(c('old','new','approximate'),each=202))
df.es$algo <- factor(df.es$algo,levels=c('old','new','approximate'))

df.mn <- df.es %>% group_by(algo) %>% summarise(md=median(es))
df.mn$lab <- paste('median: ',as.character(round(df.mn$md,2)))
df.mn$es <- c(10,10,0.4)
df.mn$y <- c(40,40,40)

df.mx <- df.es %>% group_by(algo) %>% summarise(mx=max(es))
df.mx$lab <- paste('max: ',as.character(round(df.mx$mx,2)))
df.mx$es <- c(50,50,.05)
df.mx$y <- c(38,38,38)

df.mi <- df.es %>% group_by(algo) %>% summarise(mn=min(es))
df.mi$lab <- paste('min: ',as.character(round(df.mi$mn,5)))
df.mi$es <- c(50,50,.05)
df.mi$y <- c(35,35,35)

png('Figures/hist_ess_compare.png',width=800,height=400)
ggplot(df.es,aes(x=es)) + geom_histogram() + facet_grid(~algo) + scale_x_log10() +
  geom_text(data=df.mn,aes(x=es,y=y,label=lab),size=9) +
  geom_text(data=df.mi,aes(x=es,y=y,label=lab),size=9) +
  theme(text=element_text(size=32),axis.text.x = element_text(angle=90)) + xlab(expression(n[e]))
  #+ 
  #geom_text(data=df.mx,aes(x=es,y=y,label=lab)) + 
  #geom_text(data=df.mi,aes(x=es,y=y,label=lab)) 
dev.off()

df.xi <- data.frame(old=c(f.ac(ab$xiout)),new=c(f.ac(ex$xiout)),approximate=c(f.ac(ap$xiout)))
df.xi$lag <- seq(100)

df.xi <- df.xi %>% gather(algo,ac,-lag)
df.xi$algo <- factor(df.xi$algo,levels=c('old','new','approximate'))

png('Figures/acs_compare_noapprox.png',width=800,height=400)
ggplot(df.xi %>% filter(algo!='approximate'),aes(x=lag,y=ac,lty=algo)) + geom_line(size=2) +
  theme(text=element_text(size=32),legend.key.width = unit(3,"cm")) + ylim(0,1)
dev.off()


png('Figures/acs_compare.png',width=800,height=400)
ggplot(df.xi,aes(x=lag,y=ac,lty=algo)) + geom_line(size=2) +
  theme(text=element_text(size=32),legend.key.width = unit(3,"cm")) 
dev.off()

tmp <- df.es[df.es$var=='xi',]
tmp <- tmp %>% spread(algo,es)
tmp <- tmp[,-1]
rownames(tmp) <- 'ess'
print.xtable(xtable(tmp,digits = 4,caption = 'Effective samples per second, $\\log(\\xi)$',label='tab:ess_xi'),file='Figures/ess_xi.tex')

#####################################
# smaller sizes: the 200-2000 examples
#####################################

ab <- readMat(paste('Outputs/post_reg_horse_',ab.vs,'_F_200_2000.mat',sep=''))
ex <- readMat('Outputs/post_reg_horse_F_200_2000.mat')
ap <- readMat('Outputs/post_reg_horse_F_approx_4_200_2000.mat')

samps <- rbind(t(ab$betaout),t(ex$betaout),t(ap$betaout))
samps <- cbind(samps,rbind(t(ab$etaout),t(ex$etaout),t(ap$etaout)))
samps <- cbind(samps,c(c(ab$xiout),c(ex$xiout),c(ap$xiout)))
samps <- cbind(samps,c(c(ab$sigmaSqout),c(ex$sigmaSqout),c(ap$sigmaSqout)))
samps <- data.frame(samps) 
names(samps) <- c(paste('bet',seq(100),sep=''),paste('eta',seq(100),sep=''),'xi','sigma2')
samps$iter <- rep(seq(20000),3)
samps$algo <- rep(c('old','new','approximate'),each=20000)

samps <- samps %>% gather(variable,value,-iter,-algo)
samps$algo <- factor(samps$algo,levels=c('old','new','approximate'))

p1 <- ggplot(samps[samps$variable=='xi',],aes(x=iter,y=log(value))) + geom_point(size=.3) + facet_grid(~algo) +
  theme(text=element_text(size=32),axis.text.x = element_text(angle=90)) + ylab(expression(log(xi)))

p2 <- ggplot(samps[samps$variable=='sigma2',],aes(x=iter,y=log(1/value))) + geom_point(size=.3) + facet_grid(~algo) + 
  geom_hline(yintercept = log(1/4),col='red') + theme(text=element_text(size=32),axis.text.x = element_text(angle=90)) + ylab(expression(-2*log(sigma)))

png('Figures/log_sigma2_traces.png',width = 900,height=500)
ggarrange(p1,p2,nrow=2)
dev.off()


p1 <- ggplot(samps[samps$variable=='bet10',],aes(x=iter,y=value)) + geom_point(size=.3) + facet_grid(~algo) + 
  geom_hline(yintercept = 2^(-.75),col='red') + theme(text=element_text(size=32),axis.text.x = element_text(angle=90)) +
  ylim(c(-1,2)) + ylab(expression(beta[10]))

p2 <- ggplot(samps[samps$variable=='bet10',],aes(x=value)) + geom_density() + facet_wrap(~algo,scales='free') + 
  theme(text=element_text(size=32),axis.text.x = element_text(angle=90)) + xlim(c(-1,2)) + xlab(expression(beta[10]))

png('Figures/beta10_traces.png',width=900,height=500)
ggarrange(p1,p2,nrow=2)
dev.off()


ggplot(samps[samps$variable %in% c('eta1','eta10','eta20','eta30') & samps$algo=='approximate',],aes(x=iter,y=log(value))) + geom_point(size=.3) +
  facet_wrap(~variable)

# computational complexity of the approximate algorithm with independent X

covs <- matrix(0,21,2)
mses <- matrix(0,21,2)

nkeep <- 100
fl <- list.files(path='~/Documents/GitHub/horseshoe/Code/Outputs/Approx_ind',pattern='post_reg_.*')

tms <- matrix(0,20,3)
ctr <- 0
for (f in fl) {
    print(f)
    ctr <- ctr + 1
    dat <- readMat(paste('Outputs/Approx_ind/',f,sep=''))
    if (ctr==1) {
      df.bet <- data.frame(t(dat$betaout))
      names(df.bet) <- c(seq(nkeep))
      df.bet$iter <- seq(nrow(df.bet))
      n <- as.numeric(substr(f,27,30))
      df.bet$n <- n
      p <- as.numeric(substr(f,32,nchar(f)-4))
      df.bet$p <- p
      df.bet$xi <- dat$xiout
      df.bet$sigma <- sqrt(dat$sigmaSqout)
    } else {
      tmp <- data.frame(t(dat$betaout))
      names(tmp) <- c(seq(nkeep))
      tmp$iter <- seq(nrow(tmp))
      if (grepl('_',substr(f,27,30))) {
        n <- as.numeric(substr(f,27,29))
      } else {
        n <- as.numeric(substr(f,27,30))  
      }
      tmp$n <- n
      p <- as.numeric(substr(f,32,nchar(f)-4))
      tmp$p <- p
      tmp$xi <- dat$xiout
      tmp$sigma <- sqrt(dat$sigmaSqout)
      df.bet <- rbind(df.bet,tmp)
    }
    covs[ctr,1] <- dat$coverage
    mses[ctr,1] <- dat$mse
    tms[ctr,] <- c(dat$t,n,p)
}

df.tm <- data.frame(tms)
names(df.tm) <- c('t','n','p')
fit.tm <- lm(log(t)~log(n)+log(p),data=df.tm)

tab.fit <- summary(fit.tm)

print.xtable(xtable(tab.fit,digits = 4,caption = 'Regression of $\\log(t)$ on $\\log(N) + \\log(p)$ for approximate algorithm',
                    label='tab:t_approx_reg_big'),file='Figures/t_approx_reg_big.tex')


df.bet <- df.bet %>% gather(variable,value,-n,-p,-iter)

# some res variance plot
df.bet$np <- paste('N=',df.bet$n,', p=',df.bet$p,sep='')

png('Figures/sigma_marginals.png',width=900,height=600)
ggplot(df.bet[df.bet$variable=='sigma' & df.bet$iter>=5000,],aes(x=np,y=value)) + geom_violin(draw_quantiles = c(0.025,0.975)) +
  theme(text=element_text(size=24),axis.text.x=element_text(angle=90)) + xlab('(N,p)') + ylab(expression(sigma))
dev.off()

# current best : obm, cuberoot

es.bet <- df.bet[df.bet$iter>burn,] %>% group_by(variable,n,p) %>% summarise(es=ess(value,method='obm',size='cuberoot'))
es.bet$np <- paste('n=',as.character(es.bet$n),', p=',as.character(es.bet$p),sep='')
es.bet$var <- 1/es.bet$es


box.p.ind <- ggplot(es.bet,aes(x=factor(p),y=es)) + geom_boxplot() + scale_y_log10() +
  theme(axis.text.x=element_text(angle=90),text=element_text(size=32)) + ggtitle('independent design') + xlab('p') + ylab(expression(n[e]))

box.n.ind <- ggplot(es.bet,aes(x=factor(n),y=es)) + geom_boxplot() + scale_y_log10() +
  theme(axis.text.x=element_text(angle=90),text=element_text(size=32)) + ggtitle('independent design') + xlab('N') + ylab(expression(n[e]))

fit.xi <- lm(log(var)~log(n)+log(p),data=es.bet[es.bet$variable=='xi',])

#ggplot(df.bet[df.bet$variable=='xi' & df.bet$iter>1000,],aes(x=iter,y=log(value))) + geom_point(size=.3) + facet_wrap(~n+p)
#ggplot(df.bet[df.bet$variable=='18' & df.bet$iter>1000,],aes(x=iter,y=value)) + geom_point(size=.3) + facet_wrap(~n+p)
#ggplot(df.bet[df.bet$variable=='19' & df.bet$iter>1000,],aes(x=iter,y=value)) + geom_point(size=.3) + facet_wrap(~n+p)

es.bet$variable <- factor(es.bet$variable)
es.bet$N <- es.bet$n
fit <- lm(log(var)~log(N)+log(p)+variable,data=es.bet)

fit.tab <- summary(fit)$coefficients[1:3,]

print.xtable(xtable(fit.tab,caption = 'estimated parameters from regression of $-\\log(n_e)$ on $\\log(N)+\\log(p)$, independent design',
                    label='tab:reg_independent'),file = 'Figures/reg_independent.tex')


# computational cost when X is not independent

fl <- list.files(path='~/Documents/GitHub/horseshoe/Code/Outputs/Approx_cor',pattern='post_reg_.*')
tms <- matrix(0,20,3)

ctr <- 0
for (f in fl) {
    print(f)
    ctr <- ctr + 1
    dat <- readMat(paste('Outputs/Approx_cor/',f,sep=''))
    if (ctr==1) {
      df.bet <- data.frame(t(dat$betaout))
      names(df.bet) <- c(seq(nkeep))
      df.bet$iter <- seq(nrow(df.bet))
      if (grepl('_',substr(f,33,36))) {
        df.bet$n <- as.numeric(substr(f,33,35))  
        df.bet$p <- as.numeric(substr(f,37,nchar(f)-4))
      } else {
        df.bet$n <- as.numeric(substr(f,33,36))
        df.bet$p <- as.numeric(substr(f,38,nchar(f)-4))
      }
      df.bet$xi <- dat$xiout
      df.bet$sigma <- sqrt(dat$sigmaSqout)
    } else {
      tmp <- data.frame(t(dat$betaout))
      names(tmp) <- c(seq(nkeep))
      tmp$iter <- seq(nrow(tmp))
      if (grepl('_',substr(f,33,36))) {
        tmp$n <- as.numeric(substr(f,33,35))
        tmp$p <- as.numeric(substr(f,37,nchar(f)-4))
      } else {
        tmp$n <- as.numeric(substr(f,33,36))  
        tmp$p <- as.numeric(substr(f,38,nchar(f)-4))
      }
      tmp$xi <- dat$xiout
      tmp$sigma <- sqrt(dat$sigmaSqout)
      df.bet <- rbind(df.bet,tmp)
    }
    mses[ctr,2] <- dat$mse
    covs[ctr,2] <- dat$coverage
}

df.bet <- df.bet %>% gather(variable,value,-n,-p,-iter)

es.bet <- df.bet[df.bet$iter>burn,] %>% group_by(variable,n,p) %>% summarise(es=ess(value,method='obm',size='cuberoot'))
es.bet$var <- 1/es.bet$es

box.p.cor <- ggplot(es.bet,aes(x=factor(p),y=es)) + geom_boxplot() + scale_y_log10() +
  theme(axis.text.x=element_text(angle=90),text=element_text(size=32)) + ggtitle('correlated design') + xlab('p') + ylab(expression(n[e]))

box.n.cor <- ggplot(es.bet,aes(x=factor(n),y=es)) + geom_boxplot() + scale_y_log10() +
  theme(axis.text.x=element_text(angle=90),text=element_text(size=32)) + ggtitle('correlated design') + xlab('N') + ylab(expression(n[e]))

png('Figures/box_p_ess.png',width=900,height=600)
ggarrange(box.p.ind,box.p.cor,nrow=2)
dev.off()

png('Figures/box_n_ess.png',width=900,height=600)
ggarrange(box.n.ind,box.n.cor,nrow=2)
dev.off()

# ggplot(df.bet[df.bet$variable=='xi' & df.bet$iter>1000,],aes(x=iter,y=log(value))) + geom_point(size=.3) + facet_wrap(~n+p)
# ggplot(df.bet[df.bet$variable=='18' & df.bet$iter>1000,],aes(x=iter,y=value)) + geom_point(size=.3) + facet_wrap(~n+p)
# ggplot(df.bet[df.bet$variable=='19' & df.bet$iter>1000,],aes(x=iter,y=value)) + geom_point(size=.3) + facet_wrap(~n+p)

es.bet$variable <- factor(es.bet$variable)
es.bet$N <- es.bet$n
fit <- lm(log(var)~log(N)+log(p)+variable,data=es.bet)

fit.tab <- summary(fit)$coefficients[1:3,]

print.xtable(xtable(fit.tab,caption = 'estimated parameters from regression of $-\\log(n_e)$ on $\\log(N)+\\log(p)$, dependent design',
                    label='tab:reg_dependent'),file = 'Figures/reg_dependent.tex')


# performance in estimation

fl <- list.files(path='/data1/GitHub/horseshoe/Code/Outputs/Compare',pattern='post_reg_.*')

cov.comp <- matrix(0,20,3)
mse.comp <- matrix(0,20,3)
tm.total <- matrix(0,60,1)
np <- matrix(0,60,2)
algo <- vector(mode='character')

ctr <- matrix(0,3,1)
ctr0 <- 0
for (f in fl) {
  print(f)
  ctr0 <- ctr0 + 1
  
  if (!grepl('approx',f)) {
    n <- substr(f,regexpr('_F_',f)[1]+3,regexpr('_F_',f)[1]+6)
    n <- sub('_','',n)
    len.n <- nchar(n)
    print(n)
    n <- as.numeric(n)
    p <- substr(f,regexpr('_F_',f)[1]+len.n+4,regexpr('_F_',f)[1]+len.n+7)
    print(p)
    p <- as.numeric(p)
  } else {
    n <- substr(f,regexpr('_4_',f)[1]+3,regexpr('_4_',f)[1]+6)
    n <- sub('_','',n)
    len.n <- nchar(n)
    print(n)
    n <- as.numeric(n)
    p <- substr(f,regexpr('_4_',f)[1]+len.n+4,regexpr('_4_',f)[1]+len.n+7)
    print(p)
    p <- as.numeric(p)
  }
    
  np[ctr0,] <- c(n,p) 
  
  dat <- readMat(paste('Outputs/Compare/',f,sep=''))
  tm.total[ctr0] <- dat$t 
  
  if (grepl('_ab_',f)) {
    ctr[1] <- ctr[1]+1
    cov.comp[ctr[1],1] <- dat$coverage
    mse.comp[ctr[1],1] <- dat$mse
    algo[ctr0] <- 'old'
  }
  if (!grepl('_ab_',f) & !grepl('_approx_',f)) {
    ctr[2] <- ctr[2]+1
    cov.comp[ctr[2],2] <- dat$coverage
    mse.comp[ctr[2],2] <- dat$mse
    algo[ctr0] <- 'new'
  }
  if (grepl('_approx_',f)) {
    ctr[3] <- ctr[3]+1
    cov.comp[ctr[3],3] <- dat$coverage
    mse.comp[ctr[3],3] <- dat$mse
    algo[ctr0] <- 'new_approx'
  }
}

df.time <- data.frame(algo=algo,n=np[,1],p=np[,2],t=tm.total)
fit1 <- lm(log(t)~log(n)+log(p),data=df.time[df.time$algo=='old',])
fit2 <- lm(log(t)~log(n)+log(p),data=df.time[df.time$algo=='new',])
fit3 <- lm(log(t)~log(n)+log(p),data=df.time[df.time$algo=='new_approx',])

print.xtable(xtable(summary(fit1),digits = 4,caption = 'Regression of $\\log(t)$ on $\\log(N) + \\log(p), old algorithm$',
                    label='tab:t_old_reg'),file='Figures/t_old_reg.tex')

print.xtable(xtable(summary(fit2),digits = 4,caption = 'Regression of $\\log(t)$ on $\\log(N) + \\log(p), new algorithm$',
                    label='tab:t_new_reg'),file='Figures/t_new_reg.tex')

print.xtable(xtable(summary(fit3),digits = 4,caption = 'Regression of $\\log(t)$ on $\\log(N) + \\log(p), approximate algorithm$',
                    label='tab:t_approx_reg'),file='Figures/t_approx_reg.tex')

df.tm <- data.frame(dimension=c('$\\log(N)$','$\\log(p)$'),old=c(coef(fit1)[2:3]),new=c(coef(fit2)[2:3]),approximate=c(coef(fit3)[2:3]))
print.xtable(xtable(df.tm,digits = 4,caption='Estimates from regression of $\\log(t)$ on $\\log(N)+\\log(p)$.',label='tab:t_reg'),
      include.rownames = F,file='Figures/t_reg.tex',sanitize.text.function = function(x){x})

df.cov <- data.frame(cov.comp)
names(df.cov) <- c('old','new','approximate')
df.cov$scenario <- seq(20)
df.cov <- df.cov %>% gather(algo,coverage,-scenario)
df.cov$algo <- factor(df.cov$algo,levels=c('old','new','approximate'))

png('Figures/coverage.png',width=600,height=400)
ggplot(df.cov,aes(x=algo,y=coverage)) + geom_boxplot() +
  theme(text=element_text(size=32)) + ylim(0.85,1.0)
dev.off()

df.mse <- data.frame(mse.comp)
names(df.mse) <- c('old','new','approximate')
df.mse$scenario <- seq(20)
df.mse <- df.mse %>% gather(algo,mse,-scenario)
df.mse$algo <- factor(df.mse$algo,levels=c('old','new','approximate'))

png('Figures/mse.png',width=600,height=400)
ggplot(df.mse,aes(x=algo,y=mse)) + geom_boxplot() +
  theme(text=element_text(size=32))
dev.off()


# estimation of beta in the 1000,5000 case

dat <- readMat('Outputs/post_reg_horse_F_approx_4_1000_5000.mat')

bet <- dat$betaout
bet <- t(bet)
bet <- data.frame(bet)
names(bet) <- paste('beta',seq(100),sep='')
bet$iter <- seq(nmc)
bet <- bet %>% gather(variable,value,-iter)
bet$numbet <- as.numeric(substr(bet$variable,5,nchar(bet$variable)))
bet$variable <- factor(bet$variable,levels = paste('beta',seq(100),sep=''))

beta.true <- data.frame(variable=paste('beta',seq(100),sep=''),numbet=seq(100),beta=dat$BetaTrue[1:100])
beta.true$variable <- factor(beta.true$variable,levels = paste('beta',seq(100),sep=''))

png('Figures/beta_marginals_1000_5000.png',width=900,height=900)
ggplot(bet[bet$numbet<=25 & bet$iter>=5000,],aes(x=value)) + geom_density() + facet_wrap(~variable,scales='free') +
  geom_vline(data=beta.true[beta.true$numbet<=25,],aes(xintercept=beta,col='red')) + 
  theme(text=element_text(size=24),axis.text.x=element_text(angle=90),legend.position = 'none') + xlab(expression(beta))
dev.off()


# now for 5000,50000
BetaTrue <- c(2^(-seq(from=-2,to=3.5,by=.25)),rep(0,100-23))

dat <- readMat('Outputs/post_reg_horse_F_approx_4_5000_50000.mat')

bet <- dat$betaout
bet <- t(bet)
bet <- data.frame(bet)
names(bet) <- paste('beta',seq(100),sep='')
bet$iter <- seq(nmc)
bet <- bet %>% gather(variable,value,-iter)
bet$numbet <- as.numeric(substr(bet$variable,5,nchar(bet$variable)))
bet$variable <- factor(bet$variable,levels = paste('beta',seq(100),sep=''))

bet.true <- data.frame(variable=paste('beta',seq(100),sep=''),value=BetaTrue)
bet.true$numbet <- seq(100)
bet.true$variable <- factor(bet.true$variable,levels = paste('beta',seq(100),sep=''))

png('Figures/beta_marginals_5000_50000.png',width=900,height=900)
ggplot(bet[bet$numbet<=25 & bet$iter>=5000,],aes(x=value)) + geom_density() + facet_wrap(~variable,scales='free') + 
  geom_vline(data=bet.true[bet.true$numbet<=25,],aes(xintercept = value,col='red')) +
  theme(text=element_text(size=24),axis.text.x=element_text(angle=90),legend.position='none') + xlab(expression(beta))
dev.off()




