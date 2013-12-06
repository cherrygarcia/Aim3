#tmle function without using R package
tmlefluc <- function(data, txmodel, selmodel, outmodel){
    sampled.data<-data

    sampled.data$pscore.svy <- predict(glm(formula=txmodel, data=sampled.data, family="binomial"), type="response")
        
    #estimate subsample probabilities from survey sample, separately for 
    #those treated and untreated
    modelc<-glm(formula=selmodel, data=sampled.data, family="binomial")
        
    data_new0<-sampled.data
    data_new0$t<-0
    data_new1<-sampled.data
    data_new1$t<-1

    sampled.data$a1s1.ps<-predict(modelc, data=data_new1, type="response")
    sampled.data$a0s1.ps<-predict(modelc, data=data_new0, type="response")
    
    g0w<-sampled.data$a0s1.ps*(1-sampled.data$pscore.svy)*(1/sampled.data$ipsvyselwt)
    g1w<-sampled.data$a1s1.ps*sampled.data$pscore.svy*(1/sampled.data$ipsvyselwt)

    #clever covariates
    h0w<-(1- sampled.data$t)/g0w
    h1w<-sampled.data$t/g1w
    
    #map y to bounded logit scale
    sampled.data$yb1<-(sampled.data$y - min(sampled.data$y))/(max(sampled.data$y)-min(sampled.data$y))
    sampled.data$yb2<-ifelse(sampled.data$yb1<0.001, 0.001, sampled.data$yb1)
    sampled.data$yb<-ifelse(sampled.data$yb2>0.999, 0.999, sampled.data$yb2)

    modelb<-glm(formula=boundedyoutmodel ,family="binomial", data=sampled.data, subset=insubsample==1)

    q<-cbind(predict(modelb, type="link", newdata=sampled.data), predict(modelb, type="link", newdata=data_new0), predict(modelb, type="link", newdata=data_new1))

    epsilon<-coef(glm(yb ~ -1 + offset(q[,1]) + h0w + h1w , family="binomial", data=sampled.data, subset=insubsample==1 ))

    qstar<-q + c((epsilon[1]*h0w + epsilon[2]*h1w), epsilon[1]/g0w, epsilon[2]/g1w)

    #map q back to unbounded scale
    qmapped<-(plogis(qstar)*(max(sampled.data$y)-min(sampled.data$y))) + min(sampled.data$y)

    sampled.data$dif<-qmapped[,3] - qmapped[,2]
    sampled.data$num<-sampled.data$dif*sampled.data$ipsvyselwt
        
    my.ate.tmle.txsel.sel.fluc.svy<-sum(sampled.data$num)/sum(sampled.data$ipsvyselwt)
    return(my.ate.tmle.txsel.sel.fluc.svy)
}
