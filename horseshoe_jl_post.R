setwd('~/Documents/GitHub/horseshoe_jo/')

if (!require(pacman)) {install.packages('pacman')}
pacman::p_load(readr,tidyr,ggplot2,dplyr)

df <- read_csv('Outputs/julia_test.csv')

names(df) <- seq(100)
df$iter <- seq(20000)
df <- df %>% gather(variable,value,-iter)
df$variable <- factor(df$variable,levels=seq(100))

ggplot(df[as.numeric(df$variable)<=25 & df$iter>5000,],aes(x=variable,y=value)) + geom_boxplot()
ggplot(df[as.numeric(df$variable)<=25 & df$iter>5000,],aes(x=variable,y=value)) + geom_violin()

ggplot(df[as.numeric(df$variable)<=25 & df$iter>5000,],aes(x=value)) + geom_density() + facet_wrap(~variable,scales='free') 
