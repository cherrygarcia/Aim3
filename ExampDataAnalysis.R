library(harvestr)
library(parallel)
library(zoo)
library(ggplot2)
library(tmle)

detectCores()
source("exampleDataFunctions.R")

expit<-function(p){
	exp(p)/(1+exp(p))
}

trials <- 1000

## Define seeds
seed.temp <- gather(trials,seed=123)
Seed <- matrix(nrow=trials,ncol=6)
for(i in 1:trials){
 Seed[i,] <- seed.temp[[i]][2:7]
}

## First, get treatment, selection, and outcome model formulas
## for the full model, models are chosen by AIC in a stepwise algorithm
seldat<-svydat[,c(2:4, 6:12, 16:19, 21:24, 30:34,39:40,42,46,50:57,61,62,64, 68)]
txdat<-svydat[,c(2:4, 6:8, 10:12, 16,33,34,39,40,61,62,64)]
fulltxmod<-step(glm(tertscore ~  SEXF*age_cent + factor(meducat) + moth + fath + urbancat + suburb + factor(imgen) + citizen + Language + factor(racecat) + factor(region) + cmage + factor(season) + weekend + begtime, family=binomial(), data=txdat))
#k=log(n) = BIC
less1txmod<-step(glm(tertscore ~ SEXF*age_cent + factor(meducat) + moth + fath + urbancat + suburb + factor(imgen) + citizen + Language + factor(racecat) + factor(region) + cmage + factor(season) + weekend + begtime, family=binomial(), data=txdat), k=9)
less2txmod<-step(glm(tertscore ~ SEXF*age_cent + factor(meducat) + moth + fath + urbancat + suburb + factor(imgen) + citizen + Language + factor(racecat) + factor(region) + cmage + factor(season) + weekend + begtime, family=binomial(), data=txdat), k=20)
fulltxmod$formula
#the following two are the same
less1txmod$formula
less2txmod$formula

fullselmod<-step(glm(insample ~  SEXF*age_cent + factor(meducat) + moth + fath + urbancat + suburb + factor(imgen) + citizen + Language + factor(racecat) + factor(region) + cmage + factor(season) + weekend + begtime + CH33 + pc_psych_minor + d_anxiety12_NIMH2 + d_mdddys12_NIMH2 + pc_pa_minor + pc_pa_severe + pp_pa_minor + pp_pa_severe+ mothwkdichot:mothwork + fathwork + cinc + cp_CdOddh12_NIMH2 + curremp + numparenttrauma + numrx + hrslpwkndnt+ hrslpwknt + hrbdwkndmod + hrbdwkmod + smallgestage, family=binomial(), data=seldat))
#k=log(n) = BIC
less1selmod<-step(glm(insample ~  SEXF*age_cent + factor(meducat) + moth + fath + urbancat + suburb + factor(imgen) + citizen + Language + factor(racecat) + factor(region) + cmage + factor(season) + weekend + begtime + CH33 + pc_psych_minor + d_anxiety12_NIMH2 + d_mdddys12_NIMH2 + pc_pa_minor + pc_pa_severe + pp_pa_minor + pp_pa_severe+ mothwkdichot:mothwork + fathwork + cinc + cp_CdOddh12_NIMH2 + curremp + numparenttrauma + numrx + hrslpwkndnt+ hrslpwknt + hrbdwkndmod + hrbdwkmod + smallgestage, family=binomial(), data=seldat), k=9)
less2selmod<-step(glm(insample ~  SEXF*age_cent + factor(meducat) + moth + fath + urbancat + suburb + factor(imgen) + citizen + Language + factor(racecat) + factor(region) + cmage + factor(season) + weekend + begtime + CH33 + pc_psych_minor + d_anxiety12_NIMH2 + d_mdddys12_NIMH2 + pc_pa_minor + pc_pa_severe + pp_pa_minor + pp_pa_severe+ mothwkdichot:mothwork + fathwork + cinc + cp_CdOddh12_NIMH2 + curremp + numparenttrauma + numrx + hrslpwkndnt+ hrslpwknt + hrbdwkndmod + hrbdwkmod + smallgestage, family=binomial(), data=seldat), k=20)
fullselmod$formula
#I'll use the top one
less1selmod$formula
less2selmod$formula

set.seed(132)

#load data
setwd("/Users/kararudolph/Documents/PhD/NIMH/NCSA/cortisol/simulationproject")
load("exampsvydatOct2013.Rdata")
load("exampsvysubdatOct2013.Rdata")
svydat<-b
svysubdat<-svysub[!is.na(svysub$begtime),]

#1/svywts
summary(1/svydat$final_weight)
#I(A=a)(1/txwts)
svydat$pscore.svy <- predict(glm(fulltxmod$formula , data=svydat, family="binomial"), type="response")
summary(svydat$pscore.svy[svydat$tertscore==1])
summary(1-svydat$pscore.svy[svydat$tertscore==0])
#I(A=a, sub=1)(1/subwts)
svydat$pscore.subsamp1[svydat$tertscore==1]<-predict(glm(fullselmod$formula, data=svydat[svydat$tertscore==1,], family="binomial"), type="response")
summary(svydat$pscore.subsamp1[svydat$tertscore==1 & svydat$insample==1])
svydat$pscore.subsamp1[svydat$tertscore==0]<-predict(glm(fullselmod$formula, data=svydat[svydat$tertscore==0,], family="binomial"), type="response")
summary(svydat$pscore.subsamp1[svydat$tertscore==0 & svydat$insample==1])

#txselprobs
svydat$wts[svydat$tertscore==1 & svydat$insample==1]<-(1/svydat$pscore.svy[svydat$tertscore==1 & svydat$insample==1])*(1/svydat$pscore.subsamp1[svydat$tertscore==1 & svydat$insample==1])
summary(1/svydat$wts[svydat$tertscore==1 & svydat$insample==1])
svydat$wts[svydat$tertscore==0 & svydat$insample==1]<-(1/(1-svydat$pscore.svy[svydat$tertscore==0 & svydat$insample==1]))*(1/svydat$pscore.subsamp1[svydat$tertscore==0 & svydat$insample==1])
summary(1/svydat$wts[svydat$tertscore==0 & svydat$insample==1])

#txselprobs
svydat$wts[svydat$tertscore==1 & svydat$insample==1]<-svydat$final_weight[svydat$tertscore==1 & svydat$insample==1]*(1/svydat$pscore.svy[svydat$tertscore==1 & svydat$insample==1])*(1/svydat$pscore.subsamp1[svydat$tertscore==1 & svydat$insample==1])
summary(1/svydat$wts[svydat$tertscore==1 & svydat$insample==1])
svydat$wts[svydat$tertscore==0 & svydat$insample==1]<-svydat$final_weight[svydat$tertscore==0 & svydat$insample==1]*(1/(1-svydat$pscore.svy[svydat$tertscore==0 & svydat$insample==1]))*(1/svydat$pscore.subsamp1[svydat$tertscore==0 & svydat$insample==1])
summary(1/svydat$wts[svydat$tertscore==0 & svydat$insample==1])

exampmatrix<-matrix(NA, ncol=4, nrow=16)
exampmatrix[,1] <- c("Naive", "IPTW", "IPTSvyW", "IPTSubW", "IPW", "DRWLS", "TMLE", "IPWSelLess", "DRWLSSelLess", "TMLESelLess", "IPWTxLess", "DRWLSTxLess","TMLETxLess",   "IPWTxSelLess", "DRWLSTxSelLess" ,"TMLETxSelLess")
	
# Naive estimate
temp<-lm(cortrate ~ tertscore, data=svysubdat)
exampmatrix[1,2]<-summary(temp)$coef[2,1]
sd <- summary(temp)$coef[2,2]
exampmatrix[1,3] <- summary(temp)$coef[2,1]-2*sd 
exampmatrix[1,4] <- summary(temp)$coef[2,1]+2*sd

res<-unlist(mclapply(1:trials, boot.iptw , mc.cores=2))
exampmatrix[2,2]<-iptw(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[2,3] <- cov[1]
exampmatrix[2,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.iptsvyw , mc.cores=2))
exampmatrix[3,2]<-iptsvyw(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[3,3] <- cov[1]
exampmatrix[3,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.iptsubw , mc.cores=2))
exampmatrix[4,2]<-iptsubw(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[4,3] <- cov[1]
exampmatrix[4,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.ipw , mc.cores=2))
exampmatrix[5,2]<-ipw(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[5,3] <- cov[1]
exampmatrix[5,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.joffe , mc.cores=2))
exampmatrix[6,2]<-joffe(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[6,3] <- cov[1]
exampmatrix[6,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.tmle , mc.cores=2))
exampmatrix[7,2]<-tmlefunc(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[7,3] <- cov[1]
exampmatrix[7,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.ipw.selless1 , mc.cores=2))
exampmatrix[8,2]<-ipw.selless1(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[8,3] <- cov[1]
exampmatrix[8,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.joffe.selless1 , mc.cores=2))
exampmatrix[9,2]<-joffe.selless1(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[9,3] <- cov[1]
exampmatrix[9,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.tmle.selless1 , mc.cores=2))
exampmatrix[10,2]<-tmlefunc.selless1(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[10,3] <- cov[1]
exampmatrix[10,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.ipw.txless1 , mc.cores=2))
exampmatrix[11,2]<-ipw.txless1(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[11,3] <- cov[1]
exampmatrix[11,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.joffe.txless1 , mc.cores=2))
exampmatrix[12,2]<-joffe.txless1(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[12,3] <- cov[1]
exampmatrix[12,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.tmle.txless1 , mc.cores=2))
exampmatrix[13,2]<-tmlefunc.txless1(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[13,3] <- cov[1]
exampmatrix[13,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.ipw.txless1selless1 , mc.cores=2))
exampmatrix[14,2]<-ipw.txless1selless1(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[14,3] <- cov[1]
exampmatrix[14,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.joffe.txless1selless1 , mc.cores=2))
exampmatrix[15,2]<-joffe.txless1selless1(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[15,3] <- cov[1]
exampmatrix[15,4] <- cov[2]

res<-unlist(mclapply(1:trials, boot.tmle.txless1selless1 , mc.cores=2))
exampmatrix[16,2]<-tmlefunc.txless1selless1(svydat)
cov<-quantile(res, probs=c(.025, .975))
exampmatrix[16,3] <- cov[1]
exampmatrix[16,4] <- cov[2]

write.table(exampmatrix, file="ExampleDataEffectsNov1.csv", row.names=FALSE, col.names=TRUE, sep=",")