library(harvestr)
library(parallel)
library(zoo)
library(ggplot2)

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

boot.fun.tmle.txsel.sel.orig.svy1 <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	#counterfactuals
	#estimate treatment probabilities from survey sample
	#counterfactuals
	m<-glm(tertscore ~ SEXF*age_cent + urbancat + suburb + poly(cmage,2) +meducat + Language + factor(imgen) + citizen + 	factor(region) + factor(season) + weekend + factor(racecat) + begtime , family="binomial", data=sampled.data)
	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1
	
	sampled.data$a1.ps <- predict(m, newdata=data_new1, type="response")
	sampled.data$a0.ps <- predict(m, newdata=data_new0, type="response")
	
	#estimate subsample probabilities from survey sample, 
  	#separately for those treated and untreated
	modelc<-glm(insample ~ SEXF + age_cent + urbancat + suburb + racecat +  fath + fath:SEXF + moth + moth:SEXF + poly(cmage,2) + poly(cinc,2)+ meducat+ numrx +  fathwork + fathwork:SEXF + mothwkdichot:mothwork + SEXF:mothwkdichot:mothwork + Language + factor(imgen) + citizen +  factor(region) + curremp+ hrbdwkndmod  + hrbdwkmod + smallgestage + d_mdddys12_NIMH2 + d_anxiety12_NIMH2 +  cp_CdOddh12_NIMH2, data=sampled.data, family="binomial")
	
	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	sampled.data$subsvya1wt<-ifelse(sampled.data$insample==1, (1/sampled.data$a1.ps)*(1/sampled.data$a1s1.ps)*sampled.data$final_weight, 0)
	sampled.data$subsvya0wt<-ifelse(sampled.data$insample==1, (1/(1-sampled.data$a0.ps))*(1/sampled.data$a0s1.ps)*sampled.data$final_weight,0)

	modela0<-glm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data,  weights=sampled.data$subsvya0wt,  
              family="gaussian")
	modela1<-glm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data,  weights=sampled.data$subsvya1wt,  
              family="gaussian")
	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$y1hat<-predict(modela1, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(modela0, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight
		
	my.ate.tmle.txsel.sel.orig.svy <-sum(sampled.data$num)/sum(sampled.data$final_weight)
		
	return(my.ate.tmle.txsel.sel.orig.svy)
	}

tmle.txsel.sel.orig.svy1 <- function(svydat){
	sampled.data<-svydat 

	#estimate treatment probabilities from survey sample
	#counterfactuals
	m<-glm(tertscore ~ SEXF*age_cent + urbancat + suburb + poly(cmage,2) +meducat + Language + factor(imgen) + citizen + 	factor(region) + factor(season) + weekend + factor(racecat) + begtime , family="binomial", data=sampled.data)
	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1
	
	sampled.data$a1.ps <- predict(m, newdata=data_new1, type="response")
	sampled.data$a0.ps <- predict(m, newdata=data_new0, type="response")
	
	#estimate subsample probabilities from survey sample, 
  	#separately for those treated and untreated
	modelc<-glm(insample ~ SEXF + age_cent + urbancat + suburb + racecat +  fath + fath:SEXF + moth + moth:SEXF + poly(cmage,2) + poly(cinc,2) + meducat+ numrx +  fathwork + fathwork:SEXF + mothwkdichot:mothwork + SEXF:mothwkdichot:mothwork + Language + factor(imgen) + citizen +  factor(region) + curremp+ hrbdwkndmod  + hrbdwkmod + smallgestage + d_mdddys12_NIMH2 + d_anxiety12_NIMH2 +  cp_CdOddh12_NIMH2, data=sampled.data, family="binomial")
	
	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	sampled.data$subsvya1wt<-ifelse(sampled.data$insample==1, (1/sampled.data$a1.ps)*(1/sampled.data$a1s1.ps)*sampled.data$final_weight, 0)
	sampled.data$subsvya0wt<-ifelse(sampled.data$insample==1, (1/(1-sampled.data$a0.ps))*(1/sampled.data$a0s1.ps)*sampled.data$final_weight,0)

	modela0<-glm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data,  weights=sampled.data$subsvya0wt,  
              family="gaussian")
	modela1<-glm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data,  weights=sampled.data$subsvya1wt,  
              family="gaussian")
	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$y1hat<-predict(modela1, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(modela0, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight
		
	my.ate.tmle.txsel.sel.orig.svy <-sum(sampled.data$num)/sum(sampled.data$final_weight)
		
	return(my.ate.tmle.txsel.sel.orig.svy)
	}

boot.fun.tmle.txsel.sel.orig.svy2 <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]  

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(tertscore ~ SEXF*age_cent + urbancat + suburb + poly(cmage,2) +meducat + Language + factor(imgen) + citizen + 	factor(region) + factor(season) + weekend + factor(racecat) + begtime , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(insample ~ SEXF + age_cent + urbancat + suburb + racecat +  fath + fath:SEXF + moth + moth:SEXF + poly(cmage,2) + poly(cinc,2) + meducat+ numrx +  fathwork + fathwork:SEXF + mothwkdichot:mothwork + SEXF:mothwkdichot:mothwork + Language + factor(imgen) + citizen +  factor(region) + curremp+ hrbdwkndmod  + hrbdwkmod + smallgestage + d_mdddys12_NIMH2 + d_anxiety12_NIMH2 +  cp_CdOddh12_NIMH2, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(insample ~ SEXF + age_cent + urbancat + suburb + racecat +  fath + fath:SEXF + moth + moth:SEXF + poly(cmage,2) + poly(cinc,2) + meducat+ numrx +  fathwork + fathwork:SEXF + mothwkdichot:mothwork + SEXF:mothwkdichot:mothwork + Language + factor(imgen) + citizen +  factor(region) + curremp+ hrbdwkndmod  + hrbdwkmod + smallgestage + d_mdddys12_NIMH2 + d_anxiety12_NIMH2 +  cp_CdOddh12_NIMH2, data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
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

boot.fun.tmle.txsel.sel.orig.svy3 <- function(trials){
seed.run <- Seed[trials,]
	 	set.seed(seed.run, "L'Ecuyer-CMRG")

		ind<-sample(nrow(svydat), nrow(svydat), replace=TRUE)
		sampled.data<-svydat[ind,]
  
	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy<- predict(glm(tertscore ~ SEXF*age_cent + urbancat + suburb + poly(cmage,2) +meducat + Language + factor(imgen) + citizen + 	factor(region) + factor(season) + weekend + factor(racecat) + begtime , data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample
#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(insample ~ SEXF + age_cent + urbancat + suburb + racecat +  fath + fath:SEXF + moth + moth:SEXF + poly(cmage,2) + poly(cinc,2)+ meducat+ numrx +  fathwork + fathwork:SEXF + mothwkdichot:mothwork + SEXF:mothwkdichot:mothwork + Language + factor(imgen) + citizen +  factor(region) + curremp+ hrbdwkndmod  + hrbdwkmod + smallgestage + d_mdddys12_NIMH2 + d_anxiety12_NIMH2 +  cp_CdOddh12_NIMH2, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(insample ~ SEXF + age_cent + urbancat + suburb + racecat +  fath + fath:SEXF + moth + moth:SEXF + poly(cmage,2) + poly(cinc,2)+ meducat+ numrx +  fathwork + fathwork:SEXF + mothwkdichot:mothwork + SEXF:mothwkdichot:mothwork + Language + factor(imgen) + citizen +  factor(region) + curremp+ hrbdwkndmod  + hrbdwkmod + smallgestage + d_mdddys12_NIMH2 + d_anxiety12_NIMH2 +  cp_CdOddh12_NIMH2, data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
modelsamp<-glm(insample ~ SEXF + age_cent + urbancat + suburb + racecat +  fath + fath:SEXF + moth + moth:SEXF + poly(cmage,2) + poly(cinc,2)+ meducat+ numrx +  fathwork + fathwork:SEXF + mothwkdichot:mothwork + SEXF:mothwkdichot:mothwork + Language + factor(imgen) + citizen +  factor(region) + curremp+ hrbdwkndmod  + hrbdwkmod + smallgestage + d_mdddys12_NIMH2 + d_anxiety12_NIMH2 +  cp_CdOddh12_NIMH2, data=sampled.data, family="binomial")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$a1s1.ps<-predict(modelsamp, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelsamp, data=data_new0, type="response")
		
	#combine weights
	h<-cbind( ( (sampled.data$tertscore/sampled.data$pscore.svy) - (1-sampled.data$tertscore)/(1-sampled.data$pscore.svy) )*sampled.data$subsamp.wt*sampled.data$final_weight , 1/sampled.data$pscore.svy*1/sampled.data$a1s1.ps*sampled.data$final_weight, -1/(1-sampled.data$pscore.svy)*1/sampled.data$a0s1.ps*sampled.data$final_weight)
			
	model<-lm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data)

	q<-cbind(predict(model, newdata=sampled.data),predict(model, newdata=data_new1), predict(model, newdata=data_new0))

	epsilon<-coef(lm(cortrate ~ -1 + h[,1], data=sampled.data, offset=q[,1]))

	qstar <- q + epsilon*h

	sampled.data$dif<-qstar[,2] - qstar[,3]
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight
		
	my.ate.tmle.txsel.sel.orig.svy3 <-sum(sampled.data$num)/sum(sampled.data$final_weight)
		
	return(my.ate.tmle.txsel.sel.orig.svy3)	}

tmle.txsel.sel.orig.svy2 <- function(svydat){
sampled.data<-svydat

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(tertscore ~ SEXF*age_cent + urbancat + suburb + poly(cmage,2) +meducat + Language + factor(imgen) + citizen + 	factor(region) + factor(season) + weekend + factor(racecat) + begtime , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$tertscore==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(insample ~ SEXF + age_cent + urbancat + suburb + racecat +  fath + fath:SEXF + moth + moth:SEXF + poly(cmage,2) + poly(cinc,2)+ meducat+ numrx +  fathwork + fathwork:SEXF + mothwkdichot:mothwork + SEXF:mothwkdichot:mothwork + Language + factor(imgen) + citizen +  factor(region) + curremp+ hrbdwkndmod  + hrbdwkmod + smallgestage + d_mdddys12_NIMH2 + d_anxiety12_NIMH2 +  cp_CdOddh12_NIMH2, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")

	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(insample ~ SEXF + age_cent + urbancat + suburb + racecat +  fath + fath:SEXF + moth + moth:SEXF + poly(cmage,2) + poly(cinc,2) + meducat+ numrx +  fathwork + fathwork:SEXF + mothwkdichot:mothwork + SEXF:mothwkdichot:mothwork + Language + factor(imgen) + citizen +  factor(region) + curremp+ hrbdwkndmod  + hrbdwkmod + smallgestage + d_mdddys12_NIMH2 + d_anxiety12_NIMH2 +  cp_CdOddh12_NIMH2, data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
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

tmle.txsel.sel.orig.svy3 <- function(svydat){
sampled.data<-svydat
  
	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy<- predict(glm(tertscore ~ SEXF*age_cent + urbancat + suburb + poly(cmage,2) +meducat + Language + factor(imgen) + citizen + 	factor(region) + factor(season) + weekend + factor(racecat) + begtime , data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample
#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$tertscore==1]<-predict(glm(insample ~ SEXF + age_cent + urbancat + suburb + racecat +  fath + fath:SEXF + moth + moth:SEXF + poly(cmage,2) + poly(cinc,2) + meducat+ numrx +  fathwork + fathwork:SEXF + mothwkdichot:mothwork + SEXF:mothwkdichot:mothwork + Language + factor(imgen) + citizen +  factor(region) + curremp+ hrbdwkndmod  + hrbdwkmod + smallgestage + d_mdddys12_NIMH2 + d_anxiety12_NIMH2 +  cp_CdOddh12_NIMH2, data=sampled.data[sampled.data$tertscore==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==1]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==1]
	
	sampled.data$pscore.subsamp1[sampled.data$tertscore==0]<-predict(glm(insample ~ SEXF + age_cent + urbancat + suburb + racecat +  fath + fath:SEXF + moth + moth:SEXF + poly(cmage,2) + poly(cinc,2)+ meducat+ numrx +  fathwork + fathwork:SEXF + mothwkdichot:mothwork + SEXF:mothwkdichot:mothwork + Language + factor(imgen) + citizen +  factor(region) + curremp+ hrbdwkndmod  + hrbdwkmod + smallgestage + d_mdddys12_NIMH2 + d_anxiety12_NIMH2 +  cp_CdOddh12_NIMH2, data=sampled.data[sampled.data$tertscore==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$tertscore==0]<-1/sampled.data$pscore.subsamp1[sampled.data$tertscore==0]
	
modelsamp<-glm(insample ~ SEXF + age_cent + urbancat + suburb + racecat +  fath + fath:SEXF + moth + moth:SEXF + poly(cmage,2) + poly(cinc,2) + meducat+ numrx +  fathwork + fathwork:SEXF + mothwkdichot:mothwork + SEXF:mothwkdichot:mothwork + Language + factor(imgen) + citizen +  factor(region) + curremp+ hrbdwkndmod  + hrbdwkmod + smallgestage + d_mdddys12_NIMH2 + d_anxiety12_NIMH2 +  cp_CdOddh12_NIMH2, data=sampled.data, family="binomial")

	data_new0<-sampled.data
	data_new0$tertscore<-0
	data_new1<-sampled.data
	data_new1$tertscore<-1

	sampled.data$a1s1.ps<-predict(modelsamp, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelsamp, data=data_new0, type="response")
		
	#combine weights
	h<-cbind( ( (sampled.data$tertscore/sampled.data$pscore.svy) - (1-sampled.data$tertscore)/(1-sampled.data$pscore.svy) )*sampled.data$subsamp.wt*sampled.data$final_weight , 1/sampled.data$pscore.svy*1/sampled.data$a1s1.ps*sampled.data$final_weight, -1/(1-sampled.data$pscore.svy)*1/sampled.data$a0s1.ps*sampled.data$final_weight)
			
	model<-lm(cortrate ~ tertscore + begtime + SEXF*age_cent+ urbancat + suburb + cmage+ meducat + Language + imgen+ citizen + as.factor(region) + as.factor(season) + as.factor(racecat) + weekend , data=sampled.data)

	q<-cbind(predict(model, newdata=sampled.data),predict(model, newdata=data_new1), predict(model, newdata=data_new0))

	epsilon<-coef(lm(cortrate ~ -1 + h[,1], data=sampled.data, offset=q[,1]))

	qstar <- q + epsilon*h

	sampled.data$dif<-qstar[,2] - qstar[,3]
	sampled.data$num<-sampled.data$dif*sampled.data$final_weight
		
	my.ate.tmle.txsel.sel.orig.svy3 <-sum(sampled.data$num)/sum(sampled.data$final_weight)
		
	return(my.ate.tmle.txsel.sel.orig.svy3)
	}

set.seed(132)
setwd("/Users/kararudolph/Documents/PhD/NIMH/NCSA/cortisol/simulationproject")
load("exampsvydatOct2013.Rdata")
load("exampsvysubdatOct2013.Rdata")
svydat<-b
svysubdat<-svysub[!is.na(svysub$begtime),]

exampmatrix<-matrix(NA, ncol=4, nrow=4)
exampmatrix[,1] <- c("Naive", "TMLEsvy1", "TMLEsvy2", "TMLEsvy3")
	
# Naive estimate
temp<-lm(cortrate ~ tertscore, data=svysubdat)
exampmatrix[1,2]<-summary(temp)$coef[2,1]
sd <- summary(temp)$coef[2,2]
exampmatrix[1,3] <- summary(temp)$coef[2,1]-2*sd 
exampmatrix[1,4] <- summary(temp)$coef[2,1]+2*sd

boot.tmlesvy1<-unlist(mclapply(1:trials, boot.fun.tmle.txsel.sel.orig.svy1 , mc.cores=2))
exampmatrix[2,2]<-tmle.txsel.sel.orig.svy1(svydat)
cov<-quantile(boot.tmlesvy1, probs=c(.025, .975))
exampmatrix[2,3] <- cov[1]
exampmatrix[2,4] <- cov[2]

boot.tmlesvy2<-unlist(mclapply(1:trials, boot.fun.tmle.txsel.sel.orig.svy2 , mc.cores=2))
exampmatrix[3,2]<-tmle.txsel.sel.orig.svy2(svydat)
cov<-quantile(boot.tmlesvy2, probs=c(.025, .975))
exampmatrix[3,3] <- cov[1]
exampmatrix[3,4] <- cov[2]

boot.tmlesvy3<-unlist(mclapply(1:trials, boot.fun.tmle.txsel.sel.orig.svy3 , mc.cores=2))
exampmatrix[4,2]<-tmle.txsel.sel.orig.svy3(svydat)
cov<-quantile(boot.tmlesvy3, probs=c(.025, .975))
exampmatrix[4,3] <- cov[1]
exampmatrix[4,4] <- cov[2]

write.table(exampmatrix, file="ExampleDataEffectsOctTMLE123.csv", row.names=FALSE, col.names=TRUE, sep=",")

)

exampdf1<-read.csv("ExampleDataEffectsAug.csv", header=TRUE, stringsAsFactors=FALSE)
exampdf2<-read.csv("ExampleDataEffectsOctTMLE123.csv", header=TRUE, stringsAsFactors=FALSE)

exampdf1a<-rbind(exampdf1[c(1:5),], exampdf2[2:4,])
exampdf1a$method<-factor(c("Naive", "IPTW", "IPSvyW", "IPTSvyW", "IPTSW", "TMLE1", "TMLE2", "TMLE3"), levels=c("Naive", "IPTW", "IPSvyW", "IPTSvyW", "IPTSW", "TMLE1", "TMLE2", "TMLE3"))
colnames(exampdf1a)<-c("methoddes", "mean", "LCI", "UCI", "method")
plotcis<-ggplot(exampdf1a, aes(x=method, y=mean)) + 
    geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1) +
    geom_point() + theme_bw() + ylab("Cortisol Rate Difference (ng/ml/hr)") + xlab("")  + geom_hline(aes(yintercept=0)) 

pdf("ExampEst95Fig.pdf")
plotcis
dev.off()