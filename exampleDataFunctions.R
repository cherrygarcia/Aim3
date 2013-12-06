boot.joffe <- function(trials){
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
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
			
	model<-glm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data,  weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$final_weight)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
	}

boot.joffe.selless1 <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(less1selmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(less1selmod$formula,  data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
			
	model<-glm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data,  weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$final_weight)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
	}

boot.joffe.txless1 <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(less2txmod$formula, data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(fullselmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(fullselmod$formula,  data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
			
	model<-glm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data,  weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$final_weight)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
	}

boot.joffe.txless1selless1 <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(less2txmod$formula, data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(less1selmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(less1selmod$formula,  data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
			
	model<-glm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data,  weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$final_weight)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
	}

boot.tmle <- function(trials){
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

	fit.tmle<-tmle(Y=sampled.data$cortrate, A=sampled.data$tertscore, W=W, Delta=sampled.data$insample, Qform=Y ~A + W1 + W2*W3 + W4 + W5 + W6 + W7 + W8 + W9 + W10 + factor(W11) + factor(W12) + factor(W13) + factor(W14) , g1W=sampled.data$pscore.svy, pDelta1=(1/sampled.data$final_weight)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.000001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$final_weight)
	return(ate.tmle.fluc.svy)
	}

boot.tmle.selless1 <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula, data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(less1selmod$formula, data=sampled.data, family="binomial")

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

boot.tmle.txless1 <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(less2txmod$formula, data=sampled.data, family="binomial"), type="response")
	
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

boot.tmle.txless1selless1 <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(less2txmod$formula, data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(less1selmod$formula, data=sampled.data, family="binomial")

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

boot.iptw <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$t.wt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)

	}

boot.iptsvyw <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
		
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$final_weight, 0)
			
	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)

	}

boot.iptsubw <- function(trials){
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
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt, 0)
			
	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)

	}

boot.ipw <- function(trials){
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
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
			
	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)

	}


boot.ipw.selless1 <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(less1selmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(less1selmod$formula,  data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
	}

boot.ipw.txless1 <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(less2txmod$formula, data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(fullselmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(fullselmod$formula,  data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)

	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
	}

boot.ipw.txless1selless1 <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(less2txmod$formula, data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(less1selmod$formula , data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(less1selmod$formula,  data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
			
	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
	}

joffe <- function(svydat){
		sampled.data<-svydat 

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(fullselmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(fullselmod$formula, data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
			
	model<-glm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data,  weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$final_weight)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
	}

joffe.selless1 <- function(svydat){
		sampled.data<-svydat 

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula, data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(less1selmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(less1selmod$formula,  data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
			
	model<-glm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data,  weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$final_weight)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
	}

joffe.txless1 <- function(svydat){
		sampled.data<-svydat

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(less2txmod$formula , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(fullselmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(fullselmod$formula,  data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
			
	model<-glm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data,  weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$final_weight)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
	}

joffe.txless1selless1 <- function(svydat){
		sampled.data<-svydat

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(less2txmod$formula , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(less1selmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(less1selmod$formula,  data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
			
	model<-glm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data,  weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$final_weight)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
	}

tmlefunc <- function(svydat){
		sampled.data<-svydat

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

	fit.tmle<-tmle(Y=sampled.data$cortrate, A=sampled.data$tertscore, W=W, Delta=sampled.data$insample, Qform=Y ~A + W1 + W2*W3 + W4 + W5 + W6 + W7 + W8 + W9 + W10 + factor(W11) + factor(W12) + factor(W13) + factor(W14) , g1W=sampled.data$pscore.svy,pDelta1=(1/sampled.data$final_weight)*cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps), gbound=0.0001)

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$final_weight)
	return(ate.tmle.fluc.svy)
	}

tmlefunc.selless1 <- function(svydat){
		sampled.data<-svydat

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(less1selmod$formula, data=sampled.data, family="binomial")

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

tmlefunc.txless1 <- function(svydat){
		sampled.data<-svydat
	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(less2txmod$formula, data=sampled.data, family="binomial"), type="response")
	
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

tmlefunc.txless1selless1 <- function(svydat){
		sampled.data<-svydat

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(less2txmod$formula, data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(less1selmod$formula, data=sampled.data, family="binomial")

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

iptw <- function(svydat){
	sampled.data<-svydat

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$t.wt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)

	}

iptsvyw <- function(svydat){
	sampled.data<-svydat

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
		
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$final_weight, 0)
			
	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)

	}

iptsubw <- function(svydat){
	sampled.data<-svydat

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(fullselmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(fullselmod$formula, data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt, 0)
			
	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)

	}

ipw <- function(svydat){
		sampled.data<-svydat

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(fullselmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(fullselmod$formula, data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
			
	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)

	}

ipw.selless1 <-function(svydat){
		sampled.data<-svydat 

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(fulltxmod$formula, data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(less1selmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(less1selmod$formula,  data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
	}

ipw.txless1 <- function(svydat){
		sampled.data<-svydat

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(less2txmod$formula, data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(fullselmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(fullselmod$formula,  data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)

	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
	}

ipw.txless1selless1 <- function(svydat){
		sampled.data<-svydat

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(less2txmod$formula, data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(less1selmod$formula, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(less1selmod$formula,  data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$final_weight, 0)
			
	lm.subsvytrt.boot <- lm(cortrate ~ tertscore, weights=sampled.data[sampled.data$insample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
	}