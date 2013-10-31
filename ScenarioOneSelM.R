library(zoo)
library(parallel)
library(harvestr)
library(tmle)

source("TMLEfunctionsMScen1.R")
source("TMLEfunctionsMestScen1.R")

detectCores()

expit<-function(p){
	exp(p)/(1+exp(p))
}

## Define seeds
trials<-500
seed.temp <- gather(trials,seed=123)
Seed <- matrix(nrow=trials,ncol=6)
for(i in 1:trials){
 Seed[i,] <- seed.temp[[i]][2:7]
}

nsims <- 1000

#dataset with n individuals. 
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
#print("P(selection) for z1==0 and z1==1")
#print(mean(prob.sel[z1<0]))
#print(mean(prob.sel[z1>0]))
ipsvyselw<-1/prob.sel

meany0<- -3 + z2 + I(z2^2)
y0<-rnorm(n,meany0,2)
y1<-y0 + 3 + 2*z1 + rnorm(n,0,.5)

y<-ifelse(t==1, y1, y0)

# True average effect
ate <- mean(y1)-mean(y0)
#print(ate)

#unadjusted/biased estimate of ate
mean(y[t==1]) - mean(y[t==0])

# Iterate drawing samples and estimating treatment effects
effects <- bootvar <- coverage.ate <- matrix(NA, ncol=6, nrow=nsims)
colnames(effects) <- colnames(bootvar) <- colnames(coverage.ate) <- c("IPWS1",  "TMLE2S1", "TMLE3S1", "IPWS2", "TMLE2S2", "TMLE3S2")

bias.ate <- effects - ate
#print("ATE bias")
pcent.bias.ate <- (effects - ate)/ate
#print("ATE % bias")

# bounded y
yb1<-(y-min(y))/(max(y)-min(y))
yb2<-ifelse(yb1<0.001, 0.001, yb1)
yb<-ifelse(yb2>0.999, 0.999, yb2)
	#map y values onto [0,1)
	#f<-ecdf(y)
	#ytransf<-f(y)

for (i in 1:nsims) {
	print(i)
	
	ss<-rbinom(n,1, prob.sel)
	sipsw<-(sum(ss)/n)/prob.sel

	pop <- data.frame(t, y0,y1, y, yb, ipsvyselw, sipsw, z1, z2, ss, prob.trt, prob.sel)
	colnames(pop)<-c("t","y0","y1", "y", "yb", "ipsvyselwt", "sipsw", "z1", "z2", "insample", "prob.trt", "prob.svysel")
	pop$id<-index(pop)

	svysamp<-pop[ss==1,]
	
	# Selection model for selection into smaller sample 
	# Mean P(selection)=.5
	beta0 <- -2
	beta1 <- log(1.2)
	beta2 <- log(2)
	beta3 <- log(1.2)

	svysamp$prob.subsel <- expit(beta0+ beta1*svysamp$z1 + beta3*I(svysamp$z1^2) + beta2*svysamp$z2 )
	
	svysamp$ipsubselw<-1/svysamp$prob.subsel

	svysamp$insubsample<-rbinom(nrow(svysamp), 1, svysamp$prob.subsel)

	svysub <- svysamp[svysamp$insubsample==1,]

###IPW	
effects[i, "IPWS1"] <-ipwS1(svysamp)
res.ipw<-unlist(mclapply(1:trials, boot.ipwS1, mc.cores=6))
bootvar[i, "IPWS1"] <-var(res.ipw)
cov<-quantile(res.ipw, probs=c(.025, .975))
coverage.ate[i,"IPWS1"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

effects[i, "IPWS2"] <-ipwS2(svysamp)
res.ipw<-unlist(mclapply(1:trials, boot.ipwS2, mc.cores=6))
bootvar[i, "IPWS2"] <-var(res.ipw)
cov<-quantile(res.ipw, probs=c(.025, .975))
coverage.ate[i,"IPWS2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

###TMLE

effects[i, "TMLE2S1"] <-joffeS1(svysamp)
res.joffe<-unlist(mclapply(1:trials, boot.joffeS1, mc.cores=6))
bootvar[i, "TMLE2S1"] <-var(res.joffe)
cov<-quantile(res.joffe, probs=c(.025, .975))
coverage.ate[i,"TMLE2S1"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

effects[i, "TMLE3S1"] <-tmle3S1(svysamp)
res.tmle3<-unlist(mclapply(1:trials, boot.tmle3S1, mc.cores=6))
bootvar[i, "TMLE3S1"] <-var(res.tmle3)
cov<-quantile(res.tmle3, probs=c(.025, .975))
coverage.ate[i,"TMLE3S1"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

effects[i, "TMLE2S2"] <-joffeS2(svysamp)
res.joffe<-unlist(mclapply(1:trials, boot.joffeS2, mc.cores=6))
bootvar[i, "TMLE2S2"] <-var(res.joffe)
cov<-quantile(res.joffe, probs=c(.025, .975))
coverage.ate[i,"TMLE2S2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

effects[i, "TMLE3S2"] <-tmle3S2(svysamp)
res.tmle3<-unlist(mclapply(1:trials, boot.tmle3S2, mc.cores=6))
bootvar[i, "TMLE3S2"] <-var(res.tmle3)
cov<-quantile(res.tmle3, probs=c(.025, .975))
coverage.ate[i,"TMLE3S2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

}

write.table(effects, file="EffectsS1MMSel.csv", row.names=FALSE, col.names=TRUE, sep=",")
write.table(bootvar, file="VarS1MMSel.csv", row.names=FALSE, col.names=TRUE, sep=",")
write.table(coverage.ate, file="CoverageATES1MMSel.csv", row.names=FALSE, col.names=TRUE, sep=",")