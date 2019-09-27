rm(list=ls(all=T))
setwd('/data1/GitHub/horseshoe_jo/')

if (!require('pacman')) {install.packages('pacman')}
pacman::p_load(R.matlab,ggplot2,dplyr,tidyr,xtable)

ks.stat <- function(x,y) {
  s <- ks.test(x,y)
  return(s$statistic)
}

p <- 10000
n <- 1000
cor.name <- 'cor'
cor.val <- 9
delts <- c(2,3,4,5)

nd <- length(delts)
nkeep <- 100
burn <- 5000
nmc <- 20000

exact <- readMat(paste('Outputs/post_reg_horse_F_',cor.name,'_',cor.val,'_',n,'_',p,'.mat',sep=''))
bet0 <- exact$betaout[,(burn+1):nmc]
m0 <- apply(bet0,1,mean)
v0 <- apply(bet0,1,var)
bet0 <- split(bet0,seq(nkeep))

moments.cor <- matrix(0,nd+1,2)



res <- list()
ks <- matrix(0,nkeep,nd+1)
for (j in 1:nd) {
  res[[j]] <- readMat(paste('Outputs/post_reg_horse_F_approx_',delts[j],'_',cor.name,'_',cor.val,'_',n,'_',p,'.mat',sep=''))
  betd <- res[[j]]$betaout[,(burn+1):nmc]
  m <- apply(betd,1,mean)
  v <- apply(betd,1,var)
  moments.cor[j,] <- c(cor(m,m0),cor(v,v0))
  betd <- split(betd,seq(nkeep))
  tmp <- mapply(ks.stat,bet0,betd)
  ks[,j] <- tmp
}

exact.rep <- readMat(paste('Outputs/post_reg_horse_rep_F_',cor.name,'_',cor.val,'_',n,'_',p,'.mat',sep=''))
betd <- exact.rep$betaout[,(burn+1):nmc]
m <- apply(betd,1,mean)
v <- apply(betd,1,var)
moments.cor[nd+1,] <- c(cor(m,m0),cor(v,v0))
betd <- split(betd,seq(nkeep))
tmp <- mapply(ks.stat,bet0,betd)
ks[,nd+1] <- tmp

df.ks <- data.frame(ks)
names(df.ks) <- c(as.character(delts),'exact')
df.ks$idx <- seq(nkeep)
df.ks <- df.ks %>% gather(delta,ks,-idx)

ggplot(df.ks,aes(x=idx,y=log10(ks),col=factor(delta))) + geom_point() 

png(paste('Figures/delta-eval-',cor.name,'.png',sep=''),width=600,height=400)
ggplot(df.ks,aes(x=factor(delta),y=log10(ks))) + geom_violin() + theme(text=element_text(size=24)) + xlab(expression(-log[10](delta))) +
  ylab(expression(log[10](KS)))
dev.off()

df.cors <- data.frame(moments.cor)
names(df.cors) <- c('mean','variance')
df.cors$delt <- c(delts,'exact')

xtable(df.cors)
