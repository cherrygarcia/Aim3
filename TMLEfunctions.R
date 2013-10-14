require(zoo)
require(harvestr)
require(parallel)

bootTMLE1<- function(trials, data, txmodel, selmodel, outmodel){
	seed.run <- Seed[trials,] 
	set.seed(seed.run, "L'Ecuyer-CMRG")

	ind<-sample(nrow(data), nrow(data), replace=TRUE)
	sampled.data<-as.data.frame(data[ind,])
	#assign("dat",dat,envir=.GlobalEnv)  # put the dat in the global env
	#estimate treatment probabilities from survey sample
	#counterfactuals
	m<-glm(formula=txmodel, family="binomial", data=sampled.data)
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1
	
	sampled.data$a1.ps <- predict(m, newdata=data_new1, type="response")
	sampled.data$a0.ps <- predict(m, newdata=data_new0, type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(formula=selmodel, data=sampled.data, family="binomial")
	
	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")
	
	#combine weights
	sampled.data$subsvya1wt<- ifelse(sampled.data$insubsample==1, (1/sampled.data$a1.ps)*(1/sampled.data$a1s1.ps)*sampled.data$ipsvyselwt, 0)
	sampled.data$subsvya0wt<- ifelse(sampled.data$insubsample==1, (1/(1-sampled.data$a0.ps))*(1/sampled.data$a0s1.ps)*sampled.data$ipsvyselwt, 0)
	
	modela0<-glm(formula=outmodel, data=sampled.data, weights=sampled.data$subsvya0wt, family="gaussian")
	modela1<-glm(formula=outmodel, data=sampled.data, weights=sampled.data$subsvya1wt, family="gaussian")

	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1
			
	sampled.data$y1hat<-predict(modela1, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(modela0, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt
		
	my.ate.tmle.txsel.sel.orig.svy <-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
		
	return(my.ate.tmle.txsel.sel.orig.svy)

}

TMLE1<- function(data, txmodel, selmodel, outmodel){

	#estimate treatment probabilities from survey sample
	#counterfactuals
	m<-glm(formula=txmodel, family="binomial", data=data)
	data_new0<-data
	data_new0$t<-0
	data_new1<-data
	data_new1$t<-1
	
	data$a1.ps <- predict(m, newdata=data_new1, type="response")
	data$a0.ps <- predict(m, newdata=data_new0, type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	modelc<-glm(formula=selmodel, data=data, family="binomial")
	
	data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	data$a0s1.ps<-predict(modelc, data=data_new0, type="response")
	
	#combine weights
	data$subsvya1wt<- ifelse(data$insubsample==1, (1/data$a1.ps)*(1/data$a1s1.ps)*data$ipsvyselwt, 0)
	data$subsvya0wt<- ifelse(data$insubsample==1, (1/(1-data$a0.ps))*(1/data$a0s1.ps)*data$ipsvyselwt, 0)
	
	modela0<-glm(formula=outmodel, data=data, weights=data$subsvya0wt, family="gaussian")
	modela1<-glm(formula=outmodel, data=data, weights=data$subsvya1wt, family="gaussian")

	data_new0<-data
	data_new0$t<-0
	data_new1<-data
	data_new1$t<-1
			
	data$y1hat<-predict(modela1, newdata=data_new1, type="response")
	data$y0hat<-predict(modela0, newdata=data_new0, type="response")
	data$dif<-data$y1hat - data$y0hat
	data$num<-data$dif*data$ipsvyselwt
		
	my.ate.tmle.txsel.sel.orig.svy <-sum(data$num)/sum(data$ipsvyselwt)
		
	return(my.ate.tmle.txsel.sel.orig.svy)

}

## sample use : TMLE1(data=svysamp, txmodel="t ~ poly(z2,2)", selmodel="insubsample ~ poly(z1,2) + z2", outmodel="y ~t + z2 + t:poly(z1,2)")

bootTMLE2 <- function(trials, data, txmodel, selmodel, outmodel){
	seed.run <- Seed[trials,]
	set.seed(seed.run, "L'Ecuyer-CMRG")

	ind<-sample(nrow(data), nrow(data), replace=TRUE)
	sampled.data<-as.data.frame(data[ind,])

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(formula=txmodel , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(formula=selmodel, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(formula=selmodel, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(formula=outmodel, data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1
			
	sampled.data$y1hat<-predict(model, newdata=data_new1, type="response")
	sampled.data$y0hat<-predict(model, newdata=data_new0, type="response")
	sampled.data$dif<-sampled.data$y1hat - sampled.data$y0hat
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
	}

TMLE2 <- function(data, txmodel, selmodel, outmodel){
	#estimate treatment probabilities from survey sample
	data$pscore.svy <- predict(glm(formula=txmodel , data=data, family="binomial"), type="response")
	data$t.wt <- ifelse(data$t==1, 1/data$pscore.svy, 1/(1-data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	data$pscore.subsamp1[data$t==1]<-predict(glm(formula=selmodel, data=data[data$t==1,], family="binomial"), type="response")
	data$subsamp.wt[data$t==1]<-1/data$pscore.subsamp1[data$t==1]
	
	data$pscore.subsamp1[data$t==0]<-predict(glm(formula=selmodel, data=data[data$t==0,], family="binomial"), type="response")
	data$subsamp.wt[data$t==0]<-1/data$pscore.subsamp1[data$t==0]
	
	#combine weights
	data$subsvytrtwt<-ifelse(data$insubsample==1, data$t.wt*data$subsamp.wt*data$ipsvyselwt, 0)
			
	model<-glm(formula=outmodel, data=data, weights=data$subsvytrtwt, family="gaussian")

	data_new0<-data
	data_new0$t<-0
	data_new1<-data
	data_new1$t<-1
			
	data$y1hat<-predict(model, newdata=data_new1, type="response")
	data$y0hat<-predict(model, newdata=data_new0, type="response")
	data$dif<-data$y1hat - data$y0hat
	data$num<-data$dif*data$ipsvyselwt
		
	my.ate.tmle.txsel.sel.orig.svy2 <-sum(data$num)/sum(data$ipsvyselwt)
		
	return(my.ate.tmle.txsel.sel.orig.svy2)
	}

bootTMLE3 <- function(trials, data, txmodel, selmodel, outmodel){
	seed.run <- Seed[trials,] 
	set.seed(seed.run, "L'Ecuyer-CMRG")

	ind<-sample(nrow(data), nrow(data), replace=TRUE)
	sampled.data<-as.data.frame(data[ind,])

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2)+ z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")
	
	#combine weights
	h<-cbind( ( (sampled.data$t/sampled.data$pscore.svy) - (1-sampled.data$t)/(1-sampled.data$pscore.svy) )*sampled.data$subsamp.wt*sampled.data$ipsvyselwt , 1/sampled.data$pscore.svy*1/sampled.data$a1s1.ps*sampled.data$ipsvyselwt, -1/(1-sampled.data$pscore.svy)*1/sampled.data$a0s1.ps*sampled.data$ipsvyselwt)

	modelb<-lm(y ~ t + z2 + t:z1, data=sampled.data)

	q<-cbind(predict(modelb), predict(modelb, newdata=data_new1), predict(modelb, newdata=data_new0))

	epsilon<-coef(lm(y ~ -1 + h[,1],  data=sampled.data, offset=q[,1]))

	qstar<-q + epsilon*h
		
	sampled.data$dif<-qstar[,2] - qstar[,3]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt
		
	my.ate.tmle.txsel.sel.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(my.ate.tmle.txsel.sel.fluc.svy )
    } 

TMLE3 <- function(data, txmodel, selmodel, outmodel){
	#estimate treatment probabilities from survey sample
	data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	data$pscore.subsamp1[data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=data[data$t==1,], family="binomial"), type="response")
	data$subsamp.wt[data$t==1]<-1/data$pscore.subsamp1[data$t==1]
	
	data$pscore.subsamp1[data$t==0]<-predict(glm(insubsample ~ poly(z1,2)+ z2, data=data[data$t==0,], family="binomial"), type="response")
	data$subsamp.wt[data$t==0]<-1/data$pscore.subsamp1[data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=data, family="binomial")
	
	data_new0<-data
	data_new0$t<-0
	data_new1<-data
	data_new1$t<-1

	data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	data$a0s1.ps<-predict(modelc, data=data_new0, type="response")
	
	#combine weights
	h<-cbind( ( (data$t/data$pscore.svy) - (1-data$t)/(1-data$pscore.svy) )*data$subsamp.wt*data$ipsvyselwt , 1/data$pscore.svy*1/data$a1s1.ps*data$ipsvyselwt, -1/(1-data$pscore.svy)*1/data$a0s1.ps*data$ipsvyselwt)

	modelb<-lm(y ~ t + z2 + t:z1, data=data)

	q<-cbind(predict(modelb), predict(modelb, newdata=data_new1), predict(modelb, newdata=data_new0))

	epsilon<-coef(lm(y ~ -1 + h[,1],  data=data, offset=q[,1]))

	qstar<-q + epsilon*h
		
	data$dif<-qstar[,2] - qstar[,3]
	data$num<-data$dif*data$ipsvyselwt
		
	my.ate.tmle.txsel.sel.fluc.svy<-sum(data$num)/sum(data$ipsvyselwt)
	return(my.ate.tmle.txsel.sel.fluc.svy )
    } 

