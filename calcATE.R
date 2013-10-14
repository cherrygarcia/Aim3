library(ggplot2)
library(plyr)
library(gridExtra)
library(xtable)

## First: Scenario 1. 
## Checks for true model

expit<-function(p){
  exp(p)/(1+exp(p))
}

set.seed(132)
n=100000

z1<-rnorm(n)
z2<-rnorm(n,2,1)

# make treatment variable
#P(tx ~ 1/3)
prob.trt<-expit(-2.5 + (log(1.5)*z2) + log(1.2)*I(z2^2))
t<-rbinom(n, 1, prob.trt)

trueTxWt<-1/prob.trt

# make selection into larger sample
# Selection model (based on one in Cole & Stuart, 2010)
# Mean P(selection)=.1
beta0 <- -2.3
beta1 <- log(2)
prob.sel <- expit(beta0+beta1*z1)
ipsvyselw<-1/prob.sel

meany0<- -3 + z2 + I(z2^2)
y0<-rnorm(n,meany0,2)
y1<-y0 + 3 + 2*z1 + rnorm(n,0,.5)

y<-ifelse(t==1, y1, y0)

# True average effect
ate <- mean(y1)-mean(y0)
