library(zoo)
library(harvestr)
library(parallel)
source("TMLEfunctions.R")

detectCores()

expit<-function(p){
	exp(p)/(1+exp(p))
}

trials<-500
nsims <- 250

## Define seeds
seed.temp <- gather(trials,seed=123)
Seed <- matrix(nrow=trials,ncol=6)
for(i in 1:trials){
 Seed[i,] <- seed.temp[[i]][2:7]
}

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

meany0<-  2*z2
y0<-rnorm(n,meany0,2)
y1<-y0 + 2 + 2*z1 + I(z1^2) + rnorm(n,0,.5)

y<-ifelse(t==1, y1, y0)

# True average effect
ate <- mean(y1)-mean(y0)
#print(ate)

#unadjusted/biased estimate of ate
mean(y[t==1]) - mean(y[t==0])

# Iterate drawing samples and estimating treatment effects
effects <- coverage.ate <- matrix(NA, ncol=3, nrow=nsims)
colnames(effects) <- colnames(coverage.ate) <- c("Naive", "trueTMLE1", "trueTMLE2", "trueTMLE3", "TMLE1T1", "TMLE2T1", "TMLE3T1","TMLE1T2", "TMLE2T2", "TMLE3T2", "TMLE1S1", "TMLE2S1", "TMLE3S1","TMLE1S2", "TMLE2S2", "TMLE3S2","TMLE1O1", "TMLE2O1", "TMLE3O1","TMLE1O2", "TMLE2O2", "TMLE3O2","TMLE1O3", "TMLE2O3", "TMLE3O3")

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

trueTx<-"t ~ poly(z2,2)"
MisTxMod<-"t ~ z2"
MisTxMaj<-"t ~ z1"

trueSel<-"insubsample ~ poly(z1,2) + z2"
MisSelMod<-"insubsample ~ z1 + z2"
MisSelMaj<-"insubsample ~ z1"

trueOut<-"y ~t + z2 + t:poly(z1,2)"
MisOut1<-"y ~t + z2 + t:z1"
MisOut2<-"y ~t + z2 + poly(z1,2)"
MisOut3<-"y ~t + z2 + z1"

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

	# First naive estimate
	temp <- lm(y ~ t, data=svysub)
	effects[i,"Naive"] <- summary(temp)$coef[2,1]	
	sd <- summary(temp)$coef[2,2]
	coverage.ate[i,"Naive"] <- ifelse(effects[i,"Naive"]-2*sd < ate & effects[i,"Naive"] +2*sd > ate, 1, 0)

## True model	
	effects[i, "trueTMLE1"] <- TMLE1(data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=trueOut)
	boot.tmle1<-unlist(mclapply(1:trials, bootTMLE1(trials=trials, data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle1, probs=c(.025, .975))
	coverage.ate[i,"trueTMLE1"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "trueTMLE2"] <- TMLE2(data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=trueOut)
	boot.tmle2<-unlist(mclapply(1:trials, bootTMLE2(trials=trials, data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle2, probs=c(.025, .975))
	coverage.ate[i,"trueTMLE2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "trueTMLE3"] <- TMLE3(data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=trueOut)
	boot.tmle3<-unlist(mclapply(1:trials, bootTMLE3(trials=trials, data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle3, probs=c(.025, .975))
	coverage.ate[i,"trueTMLE3"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

## MisTxMod
	effects[i, "TMLE1T1"] <- TMLE1(data=svysamp, txmodel=MisTxMod, selmodel=trueSel, outmodel=trueOut)
	boot.tmle1<-unlist(mclapply(1:trials, bootTMLE1(trials=trials, data=svysamp, txmodel=MisTxMod, selmodel=trueSel, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle1, probs=c(.025, .975))
	coverage.ate[i,"TMLE1T1"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE2T1"] <- TMLE2(data=svysamp, txmodel=MisTxMod, selmodel=trueSel, outmodel=trueOut)
	boot.tmle2<-unlist(mclapply(1:trials, bootTMLE2(trials=trials, data=svysamp, txmodel=MisTxMod, selmodel=trueSel, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle2, probs=c(.025, .975))
	coverage.ate[i,"TMLE2T1"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE3T1"] <- TMLE3(data=svysamp, txmodel=MisTxMod, selmodel=trueSel, outmodel=trueOut)
	boot.tmle3<-unlist(mclapply(1:trials, bootTMLE3(trials=trials, data=svysamp, txmodel=MisTxMod, selmodel=trueSel, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle3, probs=c(.025, .975))
	coverage.ate[i,"TMLE3T1"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	## MisTxMaj
	effects[i, "TMLE1T2"] <- TMLE1(data=svysamp, txmodel=MisTxMaj, selmodel=trueSel, outmodel=trueOut)
	boot.tmle1<-unlist(mclapply(1:trials, bootTMLE1(trials=trials, data=svysamp, txmodel=MisTxMaj, selmodel=trueSel, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle1, probs=c(.025, .975))
	coverage.ate[i,"TMLE1T2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE2T2"] <- TMLE2(data=svysamp, txmodel=MisTxMaj, selmodel=trueSel, outmodel=trueOut)
	boot.tmle2<-unlist(mclapply(1:trials, bootTMLE2(trials=trials, data=svysamp, txmodel=MisTxMaj, selmodel=trueSel, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle2, probs=c(.025, .975))
	coverage.ate[i,"TMLE2T2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE3T2"] <- TMLE3(data=svysamp, txmodel=MisTxMaj, selmodel=trueSel, outmodel=trueOut)
	boot.tmle3<-unlist(mclapply(1:trials, bootTMLE3(trials=trials, data=svysamp, txmodel=MisTxMaj, selmodel=trueSel, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle3, probs=c(.025, .975))
	coverage.ate[i,"TMLE3T2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	## MisSelMod
	effects[i, "TMLE1S1"] <- TMLE1(data=svysamp, txmodel=trueTx, selmodel=MisSelMod, outmodel=trueOut)
	boot.tmle1<-unlist(mclapply(1:trials, bootTMLE1(trials=trials, data=svysamp, txmodel=trueTx, selmodel=MisSelMod, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle1, probs=c(.025, .975))
	coverage.ate[i,"TMLE1S1"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE2S1"] <- TMLE2(data=svysamp, txmodel=trueTx, selmodel=MisSelMod, outmodel=trueOut)
	boot.tmle2<-unlist(mclapply(1:trials, bootTMLE2(trials=trials, data=svysamp, txmodel=trueTx, selmodel=MisSelMod, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle2, probs=c(.025, .975))
	coverage.ate[i,"TMLE2S1"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE3S1"] <- TMLE3(data=svysamp, txmodel=trueTx, selmodel=MisSelMod, outmodel=trueOut)
	boot.tmle3<-unlist(mclapply(1:trials, bootTMLE3(trials=trials, data=svysamp, txmodel=trueTx, selmodel=MisSelMod, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle3, probs=c(.025, .975))
	coverage.ate[i,"TMLE3S1"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	## MisSelMaj
	effects[i, "TMLE1S2"] <- TMLE1(data=svysamp, txmodel=trueTx, selmodel=MisSelMaj, outmodel=trueOut)
	boot.tmle1<-unlist(mclapply(1:trials, bootTMLE1(trials=trials, data=svysamp, txmodel=trueTx, selmodel=MisSelMaj, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle1, probs=c(.025, .975))
	coverage.ate[i,"TMLE1S2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE2S2"] <- TMLE2(data=svysamp, txmodel=trueTx, selmodel=MisSelMaj, outmodel=trueOut)
	boot.tmle2<-unlist(mclapply(1:trials, bootTMLE2(trials=trials, data=svysamp, txmodel=trueTx, selmodel=MisSelMaj, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle2, probs=c(.025, .975))
	coverage.ate[i,"TMLE2S2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE3S2"] <- TMLE3(data=svysamp, txmodel=trueTx, selmodel=MisSelMaj, outmodel=trueOut)
	boot.tmle3<-unlist(mclapply(1:trials, bootTMLE3(trials=trials, data=svysamp, txmodel=trueTx, selmodel=MisSelMaj, outmodel=trueOut), mc.cores=6))
	cov<-quantile(boot.tmle3, probs=c(.025, .975))
	coverage.ate[i,"TMLE3S2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)
	
	## MisOut1 model	
	effects[i, "TMLE1O1"] <- TMLE1(data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut1)
	boot.tmle1<-unlist(mclapply(1:trials, bootTMLE1(trials=trials, data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut1), mc.cores=6))
	cov<-quantile(boot.tmle1, probs=c(.025, .975))
	coverage.ate[i,"TMLE1O1"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE2O1"] <- TMLE2(data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut1)
	boot.tmle2<-unlist(mclapply(1:trials, bootTMLE2(trials=trials, data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut1), mc.cores=6))
	cov<-quantile(boot.tmle2, probs=c(.025, .975))
	coverage.ate[i,"TMLE2O1"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE3O1"] <- TMLE3(data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut1)
	boot.tmle3<-unlist(mclapply(1:trials, bootTMLE3(trials=trials, data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut1), mc.cores=6))
	cov<-quantile(boot.tmle3, probs=c(.025, .975))
	coverage.ate[i,"TMLE3O1"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	## MisOut2 model	
	effects[i, "TMLE1O2"] <- TMLE1(data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut2)
	boot.tmle1<-unlist(mclapply(1:trials, bootTMLE1(trials=trials, data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut2), mc.cores=6))
	cov<-quantile(boot.tmle1, probs=c(.025, .975))
	coverage.ate[i,"TMLE1O2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE2O2"] <- TMLE2(data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut2)
	boot.tmle2<-unlist(mclapply(1:trials, bootTMLE2(trials=trials, data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut2), mc.cores=6))
	cov<-quantile(boot.tmle2, probs=c(.025, .975))
	coverage.ate[i,"TMLE2O2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE3O2"] <- TMLE3(data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut2)
	boot.tmle3<-unlist(mclapply(1:trials, bootTMLE3(trials=trials, data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut2), mc.cores=6))
	cov<-quantile(boot.tmle3, probs=c(.025, .975))
	coverage.ate[i,"TMLE3O2"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	## MisOut3 model	
	effects[i, "TMLE1O3"] <- TMLE1(data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut3)
	boot.tmle1<-unlist(mclapply(1:trials, bootTMLE1(trials=trials, data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut3), mc.cores=6))
	cov<-quantile(boot.tmle1, probs=c(.025, .975))
	coverage.ate[i,"TMLE1O3"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE2O3"] <- TMLE2(data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut3)
	boot.tmle2<-unlist(mclapply(1:trials, bootTMLE2(trials=trials, data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut3), mc.cores=6))
	cov<-quantile(boot.tmle2, probs=c(.025, .975))
	coverage.ate[i,"TMLE2O3"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

	effects[i, "TMLE3O3"] <- TMLE3(data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut3)
	boot.tmle3<-unlist(mclapply(1:trials, bootTMLE3(trials=trials, data=svysamp, txmodel=trueTx, selmodel=trueSel, outmodel=MisOut3), mc.cores=6))
	cov<-quantile(boot.tmle3, probs=c(.025, .975))
	coverage.ate[i,"TMLE3O3"] <-  ifelse(cov[1] < ate & cov[2] > ate, 1, 0)

}

write.table(effects, file="EffectsS2MM.csv", row.names=FALSE, col.names=TRUE, sep=",")
write.table(coverage.ate, file="CoverageATES2MM.csv", row.names=FALSE, col.names=TRUE, sep=",")