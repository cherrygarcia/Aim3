## These functions were used to explore why the TMLE CI are larger than IPW in the example data when we didn't find this to be true in the simulation. Evidence points to sensitivity due to positivity/variable weights.

library(harvestr)
library(parallel)
library(zoo)
library(ggplot2)
library(tmle)

detectCores()

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

## load data, get model formulas

#this excludes individuals with small tx-sel probabilities
boot.tmle.excl<- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(fullselmod$formula, data=sampled.data, family="binomial")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	sampled.data$probs<-ifelse(sampled.data$tertscore==1, sampled.data$pscore.svy*sampled.data$a1s1.ps , (1-sampled.data$pscore.svy)*sampled.data$a0s1.ps )
	bounds<-quantile(sampled.data$probs, probs=c(.025, 0.975))

	subsetdat<-sampled.data[sampled.data$probs>bounds[1] & sampled.data$probs<bounds[2],]

	Wa<-data.frame(cbind(subsetdat$begtime, subsetdat$SEXF, subsetdat$age_cent, subsetdat$urbancat, subsetdat$suburb, subsetdat$cmage, subsetdat$meducat, subsetdat$Language, subsetdat$imgen, subsetdat$citizen, subsetdat$region, subsetdat$season, subsetdat$racecat, subsetdat$weekend))
	colnames(Wa)<-c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10","W11", "W12", "W13", "W14")
	Wb<-as.matrix(Wa)
	W<-apply(Wb, c(1,2), function(x) as.numeric(x))

	W1<-W[,1]
	W2<-W[,2]
	W3<-W[,3]
	W4<-W[,4]
	W5<-W[,5]
	W6<-W[,6]
	W7<-W[,7]
	W8<-W[,8]
	W9<-W[,9]
	W10<-W[,10]
	W11<-W[,11]
	W12<-W[,12]
	W13<-W[,13]
	W14<-W[,14]

	fit.tmle<-tmle(Y=subsetdat$cortrate, A=subsetdat$tertscore, W=W, Delta=subsetdat$insample, Qform=Y ~A + W1 + W2*W3 + W4 + W5 + W6 + W7 + W8 + W9 + W10 + factor(W11) + factor(W12) + factor(W13) + factor(W14) , g1W=subsetdat$pscore.svy, pDelta1=(1/subsetdat$final_weight)*cbind(subsetdat$a0s1.ps, subsetdat$a1s1.ps), gbound=0.000001)

	subsetdat$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	subsetdat$num<-subsetdat$dif*subsetdat$final_weight

	ate.tmle.fluc.svy<-sum(subsetdat$num)/sum(subsetdat$final_weight)
	return(ate.tmle.fluc.svy)
	}

boot.ipw.excl <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(fullselmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(fullselmod$formula, data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$probs<-ifelse(sampled.data$tertscore==1, sampled.data$pscore.svy*sampled.data$pscore.subsamp1 , (1-sampled.data$pscore.svy)*sampled.data$pscore.subsamp1 )
	bounds<-quantile(sampled.data$probs, probs=c(.025, 0.975))

	subsetdat<-sampled.data[sampled.data$probs>bounds[1] & sampled.data$probs<bounds[2],]

	subsetdat$subsvytrtwt<-ifelse(subsetdat$insample==1, subsetdat$t.wt*subsetdat$subsamp.wt*subsetdat$final_weight, 0)
			
	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=subsetdat[subsetdat$insample==1,]$subsvytrtwt, data=subsetdat[subsetdat$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)

	}

boot.tmle.trunc <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(fullselmod$formula, data=sampled.data, family="binomial")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	sampled.data$probs<-ifelse(sampled.data$tertscore==1, sampled.data$pscore.svy*sampled.data$a1s1.ps , (1-sampled.data$pscore.svy)*sampled.data$a0s1.ps )
	bounds<-quantile(sampled.data$probs, probs=c(.025, 0.975))

	sampled.data$pscore.svy.trunc[sampled.data$pscore.svy>bounds[1]]<-bounds[1]
    sampled.data$pscore.svy.trunc[sampled.data$pscore.svy<bounds[2]]<-bounds[2]
    sampled.data$a1s1.ps.trunc[sampled.data$a1s1.ps>bounds[1]]<-bounds[1]
    sampled.data$a1s1.ps.trunc[sampled.data$a1s1.ps<bounds[2]]<-bounds[2]
	sampled.data$a0s1.ps.trunc[sampled.data$a0s1.ps>bounds[1]]<-bounds[1]
    sampled.data$a0s1.ps.trunc[sampled.data$a0s1.ps<bounds[2]]<-bounds[2]

	Wa<-data.frame(cbind(sampled.data$begtime, sampled.data$SEXF, sampled.data$age_cent, sampled.data$urbancat, sampled.data$suburb, sampled.data$cmage, sampled.data$meducat, sampled.data$Language, sampled.data$imgen, sampled.data$citizen, sampled.data$region, sampled.data$season, sampled.data$racecat, sampled.data$weekend))
	colnames(Wa)<-c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10","W11", "W12", "W13", "W14")
	Wb<-as.matrix(Wa)
	W<-apply(Wb, c(1,2), function(x) as.numeric(x))

	W1<-W[,1]
	W2<-W[,2]
	W3<-W[,3]
	W4<-W[,4]
	W5<-W[,5]
	W6<-W[,6]
	W7<-W[,7]
	W8<-W[,8]
	W9<-W[,9]
	W10<-W[,10]
	W11<-W[,11]
	W12<-W[,12]
	W13<-W[,13]
	W14<-W[,14]

	fit.tmle<-tmle(Y=sampled.data$cortrate, A=sampled.data$tertscore, W=W, Delta=sampled.data$insample, Qform=Y ~A + W1 + W2*W3 + W4 + W5 + W6 + W7 + W8 + W9 + W10 + factor(W11) + factor(W12) + factor(W13) + factor(W14) , g1W=sampled.data$pscore.svy.trunc, pDelta1=(1/sampled.data$final_weight)*cbind(sampled.data$a0s1.ps.trunc, sampled.data$a1s1.ps.trunc), gbound=0.000001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$final_weight)
	return(ate.tmle.fluc.svy)
	}

boot.tmle.truncwts <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(fullselmod$formula, data=sampled.data, family="binomial")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	sampled.data$probs<-ifelse(sampled.data$tertscore==1, sampled.data$pscore.svy*sampled.data$a1s1.ps*(1/sampled.data$final_weight) , (1-sampled.data$pscore.svy)*sampled.data$a0s1.ps*(1/sampled.data$final_weight) )
	sampled.data$txprobs<-ifelse(sampled.data$tertscore==1, sampled.data$pscore.svy , (1-sampled.data$pscore.svy))
	sampled.data$selprobs<-ifelse(sampled.data$tertscore==1,sampled.data$a1s1.ps*(1/sampled.data$final_weight) ,sampled.data$a0s1.ps*(1/sampled.data$final_weight) )
	bounds<-quantile(sampled.data$selprobs, probs=c(.025, 0.975))
	txbounds<-quantile(sampled.data$txprobs, probs=c(.025, 0.975))

	sampled.data$pscore.svy.trunc[sampled.data$pscore.svy>txbounds[1]]<-txbounds[1]
    sampled.data$pscore.svy.trunc[sampled.data$pscore.svy<txbounds[2]]<-txbounds[2]
    sampled.data$a1s1.ps.trunc[sampled.data$a1s1.ps>bounds[1]]<-bounds[1]
    sampled.data$a1s1.ps.trunc[sampled.data$a1s1.ps<bounds[2]]<-bounds[2]
	sampled.data$a0s1.ps.trunc[sampled.data$a0s1.ps>bounds[1]]<-bounds[1]
    sampled.data$a0s1.ps.trunc[sampled.data$a0s1.ps<bounds[2]]<-bounds[2]

	Wa<-data.frame(cbind(sampled.data$begtime, sampled.data$SEXF, sampled.data$age_cent, sampled.data$urbancat, sampled.data$suburb, sampled.data$cmage, sampled.data$meducat, sampled.data$Language, sampled.data$imgen, sampled.data$citizen, sampled.data$region, sampled.data$season, sampled.data$racecat, sampled.data$weekend))
	colnames(Wa)<-c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10","W11", "W12", "W13", "W14")
	Wb<-as.matrix(Wa)
	W<-apply(Wb, c(1,2), function(x) as.numeric(x))

	W1<-W[,1]
	W2<-W[,2]
	W3<-W[,3]
	W4<-W[,4]
	W5<-W[,5]
	W6<-W[,6]
	W7<-W[,7]
	W8<-W[,8]
	W9<-W[,9]
	W10<-W[,10]
	W11<-W[,11]
	W12<-W[,12]
	W13<-W[,13]
	W14<-W[,14]

	fit.tmle<-tmle(Y=sampled.data$cortrate, A=sampled.data$tertscore, W=W, Delta=sampled.data$insample, Qform=Y ~A + W1 + W2*W3 + W4 + W5 + W6 + W7 + W8 + W9 + W10 + factor(W11) + factor(W12) + factor(W13) + factor(W14) , g1W=sampled.data$pscore.svy.trunc, pDelta1=cbind(sampled.data$a0s1.ps.trunc, sampled.data$a1s1.ps.trunc), gbound=0.000001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$final_weight)
	return(ate.tmle.fluc.svy)
	}


boot.tmle.superlearner <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(fullselmod$formula, data=sampled.data, family="binomial")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	Wa<-data.frame(cbind(sampled.data$begtime, sampled.data$SEXF, sampled.data$age_cent, sampled.data$urbancat, sampled.data$suburb, sampled.data$cmage, factor(sampled.data$meducat), sampled.data$Language, factor(sampled.data$imgen), sampled.data$citizen, sampled.data$region, factor(sampled.data$season), sampled.data$racecat, sampled.data$weekend))
	colnames(Wa)<-c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10","W11", "W12", "W13", "W14")

	fit.tmle<-tmle(Y=sampled.data$cortrate, A=sampled.data$tertscore, W=Wa, Delta=sampled.data$insample, pDelta1=(1/sampled.data$final_weight)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.000001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$final_weight)
	return(ate.tmle.fluc.svy)
	}

boot.tmle.gboundone <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(fullselmod$formula, data=sampled.data, family="binomial")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	Wa<-data.frame(cbind(sampled.data$begtime, sampled.data$SEXF, sampled.data$age_cent, sampled.data$urbancat, sampled.data$suburb, sampled.data$cmage, sampled.data$meducat, sampled.data$Language, sampled.data$imgen, sampled.data$citizen, sampled.data$region, sampled.data$season, sampled.data$racecat, sampled.data$weekend))
	colnames(Wa)<-c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10","W11", "W12", "W13", "W14")
	Wb<-as.matrix(Wa)
	W<-apply(Wb, c(1,2), function(x) as.numeric(x))

	W1<-W[,1]
	W2<-W[,2]
	W3<-W[,3]
	W4<-W[,4]
	W5<-W[,5]
	W6<-W[,6]
	W7<-W[,7]
	W8<-W[,8]
	W9<-W[,9]
	W10<-W[,10]
	W11<-W[,11]
	W12<-W[,12]
	W13<-W[,13]
	W14<-W[,14]

	fit.tmle<-tmle(Y=sampled.data$cortrate, A=sampled.data$tertscore, W=W, Delta=sampled.data$insample, Qform=Y ~A + W1 + W2*W3 + W4 + W5 + W6 + W7 + W8 + W9 + W10 + factor(W11) + factor(W12) + factor(W13) + factor(W14) , g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$final_weight)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$final_weight)
	return(ate.tmle.fluc.svy)
	}

boot.tmle.gboundtwo<- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(fullselmod$formula, data=sampled.data, family="binomial")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	Wa<-data.frame(cbind(sampled.data$begtime, sampled.data$SEXF, sampled.data$age_cent, sampled.data$urbancat, sampled.data$suburb, sampled.data$cmage, sampled.data$meducat, sampled.data$Language, sampled.data$imgen, sampled.data$citizen, sampled.data$region, sampled.data$season, sampled.data$racecat, sampled.data$weekend))
	colnames(Wa)<-c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10","W11", "W12", "W13", "W14")
	Wb<-as.matrix(Wa)
	W<-apply(Wb, c(1,2), function(x) as.numeric(x))

	W1<-W[,1]
	W2<-W[,2]
	W3<-W[,3]
	W4<-W[,4]
	W5<-W[,5]
	W6<-W[,6]
	W7<-W[,7]
	W8<-W[,8]
	W9<-W[,9]
	W10<-W[,10]
	W11<-W[,11]
	W12<-W[,12]
	W13<-W[,13]
	W14<-W[,14]

	fit.tmle<-tmle(Y=sampled.data$cortrate, A=sampled.data$tertscore, W=W, Delta=sampled.data$insample, Qform=Y ~A + W1 + W2*W3 + W4 + W5 + W6 + W7 + W8 + W9 + W10 + factor(W11) + factor(W12) + factor(W13) + factor(W14) , g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$final_weight)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.00001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$final_weight)
	return(ate.tmle.fluc.svy)
	}

boot.tmle.gboundthree <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(fullselmod$formula, data=sampled.data, family="binomial")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	Wa<-data.frame(cbind(sampled.data$begtime, sampled.data$SEXF, sampled.data$age_cent, sampled.data$urbancat, sampled.data$suburb, sampled.data$cmage, sampled.data$meducat, sampled.data$Language, sampled.data$imgen, sampled.data$citizen, sampled.data$region, sampled.data$season, sampled.data$racecat, sampled.data$weekend))
	colnames(Wa)<-c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10","W11", "W12", "W13", "W14")
	Wb<-as.matrix(Wa)
	W<-apply(Wb, c(1,2), function(x) as.numeric(x))

	W1<-W[,1]
	W2<-W[,2]
	W3<-W[,3]
	W4<-W[,4]
	W5<-W[,5]
	W6<-W[,6]
	W7<-W[,7]
	W8<-W[,8]
	W9<-W[,9]
	W10<-W[,10]
	W11<-W[,11]
	W12<-W[,12]
	W13<-W[,13]
	W14<-W[,14]

	fit.tmle<-tmle(Y=sampled.data$cortrate, A=sampled.data$tertscore, W=W, Delta=sampled.data$insample, Qform=Y ~A + W1 + W2*W3 + W4 + W5 + W6 + W7 + W8 + W9 + W10 + factor(W11) + factor(W12) + factor(W13) + factor(W14) , g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$final_weight)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.0001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$final_weight)
	return(ate.tmle.fluc.svy)
	}

boot.tmle.smallSelModel <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(insample ~ factor(season) + curremp, data=sampled.data, family="binomial")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	Wa<-data.frame(cbind(sampled.data$begtime, sampled.data$SEXF, sampled.data$age_cent, sampled.data$urbancat, sampled.data$suburb, sampled.data$cmage, sampled.data$meducat, sampled.data$Language, sampled.data$imgen, sampled.data$citizen, sampled.data$region, sampled.data$season, sampled.data$racecat, sampled.data$weekend))
	colnames(Wa)<-c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10","W11", "W12", "W13", "W14")
	Wb<-as.matrix(Wa)
	W<-apply(Wb, c(1,2), function(x) as.numeric(x))

	W1<-W[,1]
	W2<-W[,2]
	W3<-W[,3]
	W4<-W[,4]
	W5<-W[,5]
	W6<-W[,6]
	W7<-W[,7]
	W8<-W[,8]
	W9<-W[,9]
	W10<-W[,10]
	W11<-W[,11]
	W12<-W[,12]
	W13<-W[,13]
	W14<-W[,14]

	fit.tmle<-tmle(Y=sampled.data$cortrate, A=sampled.data$tertscore, W=W, Delta=sampled.data$insample, Qform=Y ~A + W1 + W2*W3 + W4 + W5 + W6 + W7 + W8 + W9 + W10 + factor(W11) + factor(W12) + factor(W13) + factor(W14) , g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$final_weight)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.000001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$final_weight)
	return(ate.tmle.fluc.svy)
	}

boot.tmle.smallOutModel <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(fullselmod$formula, data=sampled.data, family="binomial")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	Wa<-data.frame(cbind(sampled.data$begtime, sampled.data$SEXF, sampled.data$age_cent, sampled.data$fath, sampled.data$urbancat, sampled.data$citizen, sampled.data$Language,sampled.data$region))
	colnames(Wa)<-c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8")
	Wb<-as.matrix(Wa)
	W<-apply(Wb, c(1,2), function(x) as.numeric(x))

	W1<-W[,1]
	W2<-W[,2]
	W3<-W[,3]
	W4<-W[,4]
	W5<-W[,5]
	W6<-W[,6]
	W7<-W[,7]
	W8<-W[,8]

	fit.tmle<-tmle(Y=sampled.data$cortrate, A=sampled.data$tertscore, W=W, Delta=sampled.data$insample, Qform=Y ~A + W1 + W2 + W3 + W4 + W5 + W6 + W7 + factor(W8) + A:W5 + A:W6 , g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$final_weight)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.000001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$final_weight)
	return(ate.tmle.fluc.svy)
	}

#use boot function instead of the parallel mclapply function so can get BCa confidence intervals. note that the boot function is much slower than the boot function i wrote through mclapply.
boot.fun.tmle <- function(dat, index){
	sampled.data <- dat[index,]

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(fullselmod$formula, data=sampled.data, family="binomial")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	Wa<-data.frame(cbind(sampled.data$begtime, sampled.data$SEXF, sampled.data$age_cent, sampled.data$urbancat, sampled.data$suburb, sampled.data$cmage, sampled.data$meducat, sampled.data$Language, sampled.data$imgen, sampled.data$citizen, sampled.data$region, sampled.data$season, sampled.data$racecat, sampled.data$weekend))
	colnames(Wa)<-c("W1", "W2", "W3", "W4", "W5", "W6", "W7", "W8", "W9", "W10","W11", "W12", "W13", "W14")
	Wb<-as.matrix(Wa)
	W<-apply(Wb, c(1,2), function(x) as.numeric(x))

	W1<-W[,1]
	W2<-W[,2]
	W3<-W[,3]
	W4<-W[,4]
	W5<-W[,5]
	W6<-W[,6]
	W7<-W[,7]
	W8<-W[,8]
	W9<-W[,9]
	W10<-W[,10]
	W11<-W[,11]
	W12<-W[,12]
	W13<-W[,13]
	W14<-W[,14]

	fit.tmle<-tmle(Y=sampled.data$cortrate, A=sampled.data$tertscore, W=W, Delta=sampled.data$insample, Qform=Y ~A + W1 + W2*W3 + W4 + W5 + W6 + W7 + W8 + W9 + W10 + factor(W11) + factor(W12) + factor(W13) + factor(W14) , g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$final_weight)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.000001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$final_weight)
	return(ate.tmle.fluc.svy)
	}

resboot<-boot(svydat, boot.fun.tmle, R=nrow(svydat))
cov<-boot.ci(resboot, type="bca")

res.tmle.excl<-unlist(mclapply(1:trials, boot.tmle.excl , mc.cores=2))
cov.tmle.excl<-quantile(res.tmle.excl, probs=c(.025, .975))

res.ipw.excl<-unlist(mclapply(1:trials, boot.ipw.excl , mc.cores=2))
cov.ipw.excl<-quantile(res.ipw.excl, probs=c(.025, .975))

res.tmle.trunc<-unlist(mclapply(1:trials, boot.tmle.trunc, mc.cores=2))
cov.tmle.trunc<-quantile(res.tmle.trunc, probs=c(.025, .975))

res.tmle.truncwts<-unlist(mclapply(1:trials, boot.tmle.truncwts, mc.cores=2))
cov.tmle.truncwts<-quantile(res.tmle.truncwts, probs=c(.025, .975))

res.tmle.superlearner<-unlist(mclapply(1:trials, boot.tmle.superlearner, mc.cores=2))
cov.tmle.superlearner<-quantile(res.tmle.superlearner, probs=c(.025, .975))

res.tmle.gboundone<-unlist(mclapply(1:trials, boot.tmle.gboundone, mc.cores=2))
cov.tmle.gboundone<-quantile(res.tmle.gboundone, probs=c(.025, .975))\

res.tmle.gboundtwo<-unlist(mclapply(1:trials, boot.tmle.gboundtwo, mc.cores=2))
cov.tmle.gboundtwo<-quantile(res.tmle.gboundtwo, probs=c(.025, .975))

res.tmle.gboundthree<-unlist(mclapply(1:trials, boot.tmle.gboundthree, mc.cores=2))
cov.tmle.gboundthree<-quantile(res.tmle.gboundthree, probs=c(.025, .975))

res.tmle.smallSelModel<-unlist(mclapply(1:trials, boot.tmle.smallSelModel, mc.cores=2))
cov.tmle.smallSelModel<-quantile(res.tmle.smallSelModel, probs=c(.025, .975))

res.tmle.smallOutModel<-unlist(mclapply(1:trials, boot.tmle.smallOutModel, mc.cores=2))
cov.tmle.smallOutModel<-quantile(res.tmle.smallOutModel, probs=c(.025, .975))