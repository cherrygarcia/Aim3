require(zoo)
require(parallel)

joffe <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + z2 + t:poly(z1,2), data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

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

joffeT1 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ z2, data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + z2 + t:poly(z1,2), data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

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

joffeT2 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ z1, data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + z2 + t:poly(z1,2), data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

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

joffeS1 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + z2 + t:poly(z1,2), data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

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

joffeS2 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1 , data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + z2 + t:poly(z1,2), data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

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

joffeO1 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + z2 + t:z1, data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

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

joffeO2 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2 , data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + z2 + poly(z1,2), data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

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

joffeO3 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
			
	model<-glm(y ~t + z2 + z1, data=sampled.data, weights=sampled.data$subsvytrtwt, family="gaussian")

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

tmle3<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2), data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + W2 + A:poly(W1,2), g1W=sampled.data$pscore.svy, pDelta1=as.matrix(cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps)))

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

tmle3T1<-function(svysamp){
	sampled.data<-svysamp
	sampled.data$pscore.svy <- predict(glm(t ~ z2, data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + W2 + A:poly(W1,2), g1W=sampled.data$pscore.svy, pDelta1=as.matrix(cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps)))

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

tmle3T2<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ z1, data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + W2 + A:poly(W1,2), g1W=sampled.data$pscore.svy, pDelta1=as.matrix(cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps)))

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

tmle3S1<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2), data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ z1 + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + W2 + A:poly(W1,2), g1W=sampled.data$pscore.svy, pDelta1=as.matrix(cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps)))

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

tmle3S2<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2), data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ z1, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + W2 + A:poly(W1,2), g1W=sampled.data$pscore.svy, pDelta1=as.matrix(cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps)))

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

tmle3O1<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2), data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + W2 + A:W1, g1W=sampled.data$pscore.svy, pDelta1=as.matrix(cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps)))

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

tmle3O2<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2), data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + W2 + poly(W1,2), g1W=sampled.data$pscore.svy, pDelta1=as.matrix(cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps)))

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

## TAKE TWO--this does not work
icO2<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2), data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))

	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	svysub<-sampled.data[sampled.data$insubsample==1,]

	modeld<-lm(y ~ t + z2 + poly(z1,2), data=svysub)
	
	data_new0<-svysub
	data_new0$t<-0
	data_new1<-svysub
	data_new1$t<-1

## it turns out that the qstar1 and qstar0, psi bit don't matter
	svysub$qstara<-predict(modeld, type="response")
	svysub$qstar0<-predict(modeld, data=data_new0, type="response")
	svysub$qstar1<-predict(modeld, data=data_new1, type="response")
	psi<-mean(svysub$qstar1-svysub$qstar0)

	svysub$dif<-svysub$y - svysub$qstara
	svysub$h<- ifelse(svysub$t==1, (1/(svysub$pscore.subsamp1))*(1/svysub$pscore.svy),  (-1/(svysub$pscore.subsamp1))*(1/svysub$pscore.svy) )
	svysub$ic<- svysub$h*(svysub$y - svysub$qstara) + svysub$qstar1 - svysub$qstar0 - psi
	summary(svysub$ic)
	sum(svysub$ic*svysub$ipsvyselwt)/sum(svysub$ipsvyselwt)

	svysub$h<- ifelse(svysub$t==1, svysub$ipsvyselwt*(1/(svysub$pscore.subsamp1))*(1/svysub$pscore.svy),  (-1/(svysub$pscore.subsamp1))*(1/svysub$pscore.svy)*svysub$ipsvyselwt )
	svysub$ic<- svysub$h*(svysub$y - svysub$qstara) + svysub$qstar1 - svysub$qstar0 - psi
	summary(svysub$ic)

## Try if there is no tx effect heterogeneity
## make data
set.seed(132)
n=100000

z1<-rnorm(n)
z2<-rnorm(n,2,1)

# make treatment variable
#P(tx ~ 1/3)
prob.trt<-expit(-2.5 + (log(1.5)*z2) + log(1.2)*I(z2^2))
t<-rbinom(n, 1, prob.trt)

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

meany0<-  2*z2 + z1
y0<-rnorm(n,meany0,2)
y1<-y0 + 3 + rnorm(n,0,.5)

y<-ifelse(t==1, y1, y0)

# True average effect
ate <- mean(y1)-mean(y0)
#print(ate)
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

## true
sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2), data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))

	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	svysub<-sampled.data[sampled.data$insubsample==1,]

	modeld<-lm(y ~ t + z2 + z1, data=svysub)
	
	data_new0<-svysub
	data_new0$t<-0
	data_new1<-svysub
	data_new1$t<-1

## it turns out that the qstar1 and qstar0, psi bit don't matter
	svysub$qstara<-predict(modeld, type="response")
	svysub$qstar0<-predict(modeld, data=data_new0, type="response")
	svysub$qstar1<-predict(modeld, data=data_new1, type="response")
	psi<-mean(svysub$qstar1-svysub$qstar0)

	svysub$dif<-svysub$y - svysub$qstara
	svysub$h<- ifelse(svysub$t==1, (1/(svysub$pscore.subsamp1))*(1/svysub$pscore.svy),  (-1/(svysub$pscore.subsamp1))*(1/svysub$pscore.svy) )
	svysub$ic<- svysub$h*(svysub$y - svysub$qstara) + svysub$qstar1 - svysub$qstar0 - psi
	summary(svysub$ic)
	sum(svysub$ic*svysub$ipsvyselwt)/sum(svysub$ipsvyselwt)

	svysub$h<- ifelse(svysub$t==1, svysub$ipsvyselwt*(1/(svysub$pscore.subsamp1))*(1/svysub$pscore.svy),  (-1/(svysub$pscore.subsamp1))*(1/svysub$pscore.svy)*svysub$ipsvyselwt )
	svysub$ic<- svysub$h*(svysub$y - svysub$qstara) + svysub$qstar1 - svysub$qstar0 - psi
	summary(svysub$ic)

##Weird, it also turns out that it works if there is no effect heterogeneity.
	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + W2 + W1, g1W=sampled.data$pscore.svy, pDelta1=as.matrix(cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps)))

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)

#wrong
	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + W2 , g1W=sampled.data$pscore.svy, pDelta1=as.matrix(cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps)))

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)

#wrong
	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + W1 , g1W=sampled.data$pscore.svy, pDelta1=as.matrix(cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps)))

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)

	#wrong
	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A , g1W=sampled.data$pscore.svy, pDelta1=as.matrix(cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps)))

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)

	return(ate.tmle.fluc.svy)
}

##trying with modified TMLE3. 
TMLE3 <- function(data, txmodel, selmodel, outmodel){
		data<-svysamp
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
        h<-cbind( ( (data$t/data$pscore.svy) - (1-data$t)/data$pscore.svy )*data$subsamp.wt*data$ipsvyselwt , 1/data$pscore.svy*1/data$a1s1.ps*data$ipsvyselwt, -1/data$pscore.svy*1/data$a0s1.ps*data$ipsvyselwt)

        modelb<-lm(y ~ t + z2 + z1, data=data)

        q<-cbind(predict(modelb), predict(modelb, newdata=data_new1), predict(modelb, newdata=data_new0))

        epsilon<-coef(lm(y ~ -1 + h[,1],  data=data, offset=q[,1]))

        qstar<-q + epsilon*h
        data$dif<-qstar[,2] - qstar[,3]
        data$num<-data$dif*data$ipsvyselwt
                
        my.ate.tmle.txsel.sel.fluc.svy<-sum(data$num)/sum(data$ipsvyselwt)
        return(my.ate.tmle.txsel.sel.fluc.svy )
    } 
#WRONG model
TMLE3 <- function(data, txmodel, selmodel, outmodel){
		data<-svysamp
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
        h<-cbind( ( (data$t/data$pscore.svy) - (1-data$t)/(1-data$pscore.svy ))*data$subsamp.wt*data$ipsvyselwt , 1/data$pscore.svy*1/data$a1s1.ps*data$ipsvyselwt, -1/(1-data$pscore.svy)*1/data$a0s1.ps*data$ipsvyselwt)

        modelb<-lm(y ~ t + z1, data=data)

        q<-cbind(predict(modelb), predict(modelb, newdata=data_new1), predict(modelb, newdata=data_new0))

        epsilon<-coef(lm(y ~ -1 + h[,1],  data=data, offset=q[,1]))

        qstar<-q + epsilon*h
        data$dif<-qstar[,2] - qstar[,3]
        data$num<-data$dif*data$ipsvyselwt
                
        my.ate.tmle.txsel.sel.fluc.svy<-sum(data$num)/sum(data$ipsvyselwt)
        return(my.ate.tmle.txsel.sel.fluc.svy )
    } 

tmle3O3<-function(svysamp){
	sampled.data<-svysamp

	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2), data=sampled.data, family="binomial"), type="response")
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	modelc<-glm(insubsample ~ poly(z1,2) + z2, data=sampled.data, family="binomial")
	
	data_new0<-sampled.data
	data_new0$t<-0
	data_new1<-sampled.data
	data_new1$t<-1

	sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
	sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")

	W<-as.matrix(cbind(sampled.data$z1, sampled.data$z2))
	W1<-W[,1]
	W2<-W[,2]

	fit.tmle<-tmle(Y=sampled.data$y, A=sampled.data$t, W=W, Delta=sampled.data$insubsample, Qform=Y ~A + W2 + W1, g1W=sampled.data$pscore.svy, pDelta1=as.matrix(cbind(sampled.data$a0s1.ps, sampled.data$a1s1.ps)))

	sampled.data$dif<-fit.tmle$Qstar[,2] - fit.tmle$Qstar[,1]
	sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt

	ate.tmle.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
	return(ate.tmle.fluc.svy)
}

iptw <- function(svysamp){
	sampled.data<-svysamp

	#sampled.data<-svysamp
	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$t.wt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 

ipsw <- function(svysamp){
	sampled.data<-svysamp

	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	#combine weights
	sampled.data$subsvywt<-ifelse(sampled.data$insubsample==1, sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
		
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$subsvywt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 

iptsvyw <- function(svysamp){
	sampled.data<-svysamp

	#sampled.data<-svysamp
	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))

	#combine weights
	sampled.data$txsvywt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$ipsvyselwt, 0)
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$txsvywt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    }     

ipw <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
		
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 

ipwT1 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ z2 , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
		
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 

ipwT2 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ z1 , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ poly(z1,2) + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
		
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 

ipwS1 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1+ z2, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1 + z2, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
		
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 

ipwS2 <- function(svysamp){
	sampled.data<-svysamp

	#estimate treatment probabilities from survey sample
	sampled.data$pscore.svy <- predict(glm(t ~ poly(z2,2) , data=sampled.data, family="binomial"), type="response")
	sampled.data$t.wt <- ifelse(sampled.data$t==1, 1/sampled.data$pscore.svy, 1/(1-sampled.data$pscore.svy))
	
	#estimate subsample probabilities from survey sample, separately for those treated and untreated.
	sampled.data$pscore.subsamp1[sampled.data$t==1]<-predict(glm(insubsample ~ z1, data=sampled.data[sampled.data$t==1,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==1]<-1/sampled.data$pscore.subsamp1[sampled.data$t==1]
	
	sampled.data$pscore.subsamp1[sampled.data$t==0]<-predict(glm(insubsample ~ z1, data=sampled.data[sampled.data$t==0,], family="binomial"), type="response")
	sampled.data$subsamp.wt[sampled.data$t==0]<-1/sampled.data$pscore.subsamp1[sampled.data$t==0]
	
	#combine weights
	sampled.data$subsvytrtwt<-ifelse(sampled.data$insubsample==1, sampled.data$t.wt*sampled.data$subsamp.wt*sampled.data$ipsvyselwt, 0)
		
	lm.subsvytrt.boot <- lm(y ~ t, weights=sampled.data[sampled.data$insubsample==1,]$subsvytrtwt, data=sampled.data[sampled.data$insubsample==1,])
	my.ate.ipw.boot <- summary(lm.subsvytrt.boot)$coef[2,1]

	return(my.ate.ipw.boot)
    } 
