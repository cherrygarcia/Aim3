library(ggplot2)
library(plyr)
library(gridExtra)
library(xtable)

setwd("/Users/kararudolph/Documents/PhD/NIMH/NCSA/cortisol/simulationproject/Output/3level")

### Sensitivity to incorrect correction of selection or treatment non-random mechanism under Scenario 1
## Under-correct. 

# Treatment. Treatment mechanism is ~ poly(z2,2). No correction
datTxunder<- read.csv("EffectsS1MMIncorATx.csv", header=TRUE, stringsAsFactors=FALSE)
coverTxunder<- read.csv("CoverageATES1MMIncorATx.csv", header=TRUE, stringsAsFactors=FALSE)
varTxunder<- read.csv("VarS1MMIncorATx.csv", header=TRUE, stringsAsFactors=FALSE)
datTxunderTMLE<- read.csv("EffectsS1MMIncorATxTMLERD.csv", header=TRUE, stringsAsFactors=FALSE)
coverTxunderTMLE<- read.csv("CoverageATES1MMIncorATxTMLERD.csv", header=TRUE, stringsAsFactors=FALSE)
varTxunderTMLE<- read.csv("VarS1MMIncorATxTMLERD.csv", header=TRUE, stringsAsFactors=FALSE)

# Selection. Sub-sample selection mechanism is ~ poly(z1,2) + z2. No correction
datSelunder<- read.csv("EffectsS1MMIncorA.csv", header=TRUE, stringsAsFactors=FALSE)
coverSelunder<- read.csv("CoverageATES1MMIncorA.csv", header=TRUE, stringsAsFactors=FALSE)
varSelunder<- read.csv("VarS1MMIncorA.csv", header=TRUE, stringsAsFactors=FALSE)
datSelunderTMLE<- read.csv("EffectsS1MMIncorATMLERD.csv", header=TRUE, stringsAsFactors=FALSE)
coverSelunderTMLE<- read.csv("CoverageATES1MMIncorATMLERD.csv", header=TRUE, stringsAsFactors=FALSE)
varSelunderTMLE<- read.csv("VarS1MMIncorATMLERD.csv", header=TRUE, stringsAsFactors=FALSE)

underS1box<-c(dats1$trueIPSW, datTxunder$TMLE2, datTxunderTMLE$TMLE3, dats1$trueIPTSvyW, datSelunder$TMLE2, datSelunderTMLE$TMLE3)
underS1cover<-c(covers1$trueIPSW, coverTxunder$TMLE2, coverTxunderTMLE$TMLE3, covers1$trueIPTSvyW, coverSelunder$TMLE2, coverSelunderTMLE$TMLE3)
underS1var<-c(vars1$trueIPSW, varTxunder$TMLE2, varTxunderTMLE$TMLE3, vars1$trueIPTSvyW, varSelunder$TMLE2, varSelunderTMLE$TMLE3)

labinc<-c( rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000) )
catinc<-c(rep("tx", 3000), rep("sel", 3000) )

underS1box.dat<-data.frame(underS1box, labinc, catinc)
underS1cover.dat<-data.frame(underS1cover, labinc, catinc)
underS1var.dat<-data.frame(underS1var, labinc, catinc)

sum.underS1 <- ddply(underS1box.dat, .(labinc, catinc), summarise, 
  AbsBias = mean(underS1box-ate),
    PcentBias = mean((underS1box-ate)/ate)*100,
     MSE=mean((underS1box - ate)^2))

underS1.box.plot<-ggplot(underS1box.dat, aes(x=labinc, y=underS1box)) + geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", shape=5, size=4) + facet_grid(~ catinc ) +
    coord_flip() + geom_hline(yintercept=ate, colour="red")  + scale_x_discrete(limits=c("IPW", "DRWLS", "TMLE")) + ylab("ATE") +xlab("") + ggtitle("ATE, Scenario 1")+ theme( title=element_text(size=9))

sumcover.underS1<- ddply(underS1cover.dat, .(labinc, catinc), summarise, 
  meancov = mean(underS1cover))
sumvar.underS1<- ddply(underS1var.dat, .(labinc, catinc), summarise, 
  meanvar = mean(underS1var))

underS1.cov.plot<-ggplot(data=underS1cover.dat, aes(x=labinc, y=underS1cover))  +facet_grid(~ catinc ) +
    geom_bar( stat="identity", fill="lightblue") +
scale_x_discrete(limits=c("IPW", "DRWLS", "TMLE")) + geom_text(data=sumcover.underS1, aes(x=labinc, y=meancov*100, label=meancov*100), vjust=0, size=3) +
   ggtitle("Coverage, Scenario 1") + xlab("") + ylab("") + theme(axis.text.x=element_text(size=7), title=element_text(size=9))

grid.arrange(underS1.box.plot, underS1.cov.plot, nrow=1)

sum.underS1
sumcover.underS1
sumvar.underS1

## Over-correct

# Treatment 1. Treatment mechanism is ~ z2. Correction mechanism is: ~poly(z2,2)
datTxover1<- read.csv("EffectsS1MMTxOver1.csv", header=TRUE, stringsAsFactors=FALSE)
coverTxover1<- read.csv("CoverageATES1MMTxOver1.csv", header=TRUE, stringsAsFactors=FALSE)
varTxover1<- read.csv("VarS1MMTxOver1.csv", header=TRUE, stringsAsFactors=FALSE)

datTxover1TMLE<- read.csv("EffectsS1MMTxOver1TMLERD.csv", header=TRUE, stringsAsFactors=FALSE)
coverTxover1TMLE<- read.csv("CoverageATES1MMTxOver1TMLERD.csv", header=TRUE, stringsAsFactors=FALSE)
varTxover1TMLE<- read.csv("VarS1MMTxOver1TMLERD.csv", header=TRUE, stringsAsFactors=FALSE)

# Treatment 2. Treatment mechanism is ~ 1. Correction mechanism is: ~poly(z2,2)
datTxover2<- read.csv("EffectsS1MMTxOver2.csv", header=TRUE, stringsAsFactors=FALSE)
coverTxover2<- read.csv("CoverageATES1MMTxOver2.csv", header=TRUE, stringsAsFactors=FALSE)
varTxover2<- read.csv("VarS1MMTxOver2.csv", header=TRUE, stringsAsFactors=FALSE)

datTxover2TMLE<- read.csv("EffectsS1MMTxOver2TMLERD.csv", header=TRUE, stringsAsFactors=FALSE)
coverTxover2TMLE<- read.csv("CoverageATES1MMTxOver2TMLERD.csv", header=TRUE, stringsAsFactors=FALSE)
varTxover2TMLE<- read.csv("VarS1MMTxOver2TMLERD.csv", header=TRUE, stringsAsFactors=FALSE)

datTxover2ipw<- read.csv("EffectsS1MMTxOver2rep2.csv", header=TRUE, stringsAsFactors=FALSE)
coverTxover2ipw<- read.csv("CoverageATES1MMTxOver2rep2.csv", header=TRUE, stringsAsFactors=FALSE)
varTxover2ipw<- read.csv("VarS1MMTxOver2rep2.csv", header=TRUE, stringsAsFactors=FALSE)

# Selection 1. Sub-sample selection mechanism is ~ poly(z1,2). Correction mechanism is:  ~ poly(z1,2) + z2
datSelover1<- read.csv("EffectsS1MMOver1.csv", header=TRUE, stringsAsFactors=FALSE)
coverSelover1<- read.csv("CoverageATES1MMOver1.csv", header=TRUE, stringsAsFactors=FALSE)
varSelover1<- read.csv("VarS1MMOver1.csv", header=TRUE, stringsAsFactors=FALSE)

datSelover1TMLE<- read.csv("EffectsS1MMOver1TMLERD.csv", header=TRUE, stringsAsFactors=FALSE)
coverSelover1TMLE<- read.csv("CoverageATES1MMOver1TMLERD.csv", header=TRUE, stringsAsFactors=FALSE)
varSelover1TMLE<- read.csv("VarS1MMOver1TMLERD.csv", header=TRUE, stringsAsFactors=FALSE)

# Selection 2. Sub-sample selection mechanism is ~ 1. Correction mechanism is:  ~ poly(z1,2) + z2
datSelover2<- read.csv("EffectsS1MMOver2.csv", header=TRUE, stringsAsFactors=FALSE)
coverSelover2<- read.csv("CoverageATES1MMOver2.csv", header=TRUE, stringsAsFactors=FALSE)
varSelover2<- read.csv("VarS1MMOver2.csv", header=TRUE, stringsAsFactors=FALSE)

datSelover2TMLE<- read.csv("EffectsS1MMOver2TMLERD.csv", header=TRUE, stringsAsFactors=FALSE)
coverSelover2TMLE<- read.csv("CoverageATES1MMOver2TMLERD.csv", header=TRUE, stringsAsFactors=FALSE)
varSelover2TMLE<- read.csv("VarS1MMOver2TMLERD.csv", header=TRUE, stringsAsFactors=FALSE)

txoverS1box<-c(datTxover1$IPW, datTxover1$TMLE2, datTxover1TMLE$TMLE3,
  datTxover1$IPWover, datTxover1$TMLE2over, datTxover1TMLE$TMLE3over,
datTxover2ipw$IPW, datTxover2$TMLE2, datTxover2TMLE$TMLE3,
datTxover2$IPWover, datTxover2$TMLE2over, datTxover2TMLE$TMLE3over)

txoverS1cover<-c(coverTxover1$IPW, coverTxover1$TMLE2, coverTxover1TMLE$TMLE3,
  coverTxover1$IPWover, coverTxover1$TMLE2over, coverTxover1TMLE$TMLE3over,
coverTxover2ipw$IPW, coverTxover2$TMLE2, coverTxover2TMLE$TMLE3,
coverTxover2$IPWover, coverTxover2$TMLE2over, coverTxover2TMLE$TMLE3over)

txoverS1var<-c(varTxover1$IPW, varTxover1$TMLE2, varTxover1TMLE$TMLE3,
  varTxover1$IPWover, varTxover1$TMLE2over, varTxover1TMLE$TMLE3over,
varTxover2ipw$IPW, varTxover2$TMLE2, varTxover2TMLE$TMLE3,
varTxover2$IPWover, varTxover2$TMLE2over, varTxover2TMLE$TMLE3over)

seloverS1box<-c(datSelover1$IPW, datSelover1$TMLE2, datSelover1TMLE$TMLE3,
  datSelover1$IPWover, datSelover1$TMLE2over, datSelover1TMLE$TMLE3over,
datSelover2$IPW, datSelover2$TMLE2, datSelover2TMLE$TMLE3,
datSelover2$IPWover, datSelover2$TMLE2over, datSelover2TMLE$TMLE3over)

seloverS1cover<-c(coverSelover1$IPW, coverSelover1$TMLE2, coverSelover1TMLE$TMLE3,
  coverSelover1$IPWover, coverSelover1$TMLE2over, coverSelover1TMLE$TMLE3over,
coverSelover2$IPW, coverSelover2$TMLE2, coverSelover2TMLE$TMLE3,
coverSelover2$IPWover, coverSelover2$TMLE2over, coverSelover2TMLE$TMLE3over)

seloverS1var<-c(varSelover1$IPW, varSelover1$TMLE2, varSelover1TMLE$TMLE3,
  varSelover1$IPWover, varSelover1$TMLE2over, varSelover1TMLE$TMLE3over,
varSelover2$IPW, varSelover2$TMLE2, varSelover2TMLE$TMLE3,
varSelover2$IPWover, varSelover2$TMLE2over, varSelover2TMLE$TMLE3over)

labincover<-c( rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000) , rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000) )
catincover<-c(rep("trueMod", 3000), rep("ModMis", 3000),  rep("trueNone", 3000), rep("NoneMis", 3000))

txoverS1box.dat<-data.frame(txoverS1box, labincover, catincover)
txoverS1cover.dat<-data.frame(txoverS1cover, labincover, catincover)
txoverS1var.dat<-data.frame(txoverS1var, labincover, catincover)

seloverS1box.dat<-data.frame(seloverS1box, labincover, catincover)
seloverS1cover.dat<-data.frame(seloverS1cover, labincover, catincover)
seloverS1var.dat<-data.frame(seloverS1var, labincover, catincover)

sum.txoverS1 <- ddply(txoverS1box.dat, .(labincover, catincover), summarise, 
  AbsBias = mean(txoverS1box-ate),
    PcentBias = mean((txoverS1box-ate)/ate)*100,
     MSE=mean((txoverS1box - ate)^2))

txoverS1.box.plot<-ggplot(txoverS1box.dat, aes(x=catincover, y=txoverS1box)) + geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", shape=5, size=4) + facet_grid(~ labincover ) +
    coord_flip() + geom_hline(yintercept=ate, colour="red")  + scale_x_discrete(limits=c("trueMod", "ModMis", "trueNone", "NoneMis")) + ylab("ATE") +xlab("") + ggtitle("ATE, Scenario 1, Tx Model")+ theme( title=element_text(size=9))

sumcover.txoverS1<- ddply(txoverS1cover.dat, .(labincover, catincover), summarise, 
  meancov = mean(txoverS1cover)*100)
sumvar.txoverS1<- ddply(txoverS1var.dat, .(labincover, catincover), summarise, 
  meanvar = mean(txoverS1var))

txoverS1.cov.plot<-ggplot(data=txoverS1cover.dat, aes(x=catincover, y=txoverS1cover))  +facet_grid(~ labincover ) +
    geom_bar( stat="identity", fill="lightblue") +
scale_x_discrete(limits=c("trueMod", "ModMis", "trueNone", "NoneMis")) + geom_text(data=sumcover.txoverS1, aes(x=catincover, y=meancov*100, label=meancov*100), vjust=0, size=3) +
   ggtitle("Coverage, Scenario 1, Tx Model") + xlab("") + ylab("") + theme(axis.text.x=element_text(size=7), title=element_text(size=9))

grid.arrange(txoverS1.box.plot, txoverS1.cov.plot, nrow=1)

sum.txoverS1
sumcover.txoverS1
sumvar.txoverS1

sum.seloverS1 <- ddply(seloverS1box.dat, .(labincover, catincover), summarise, 
  AbsBias = mean(seloverS1box-ate),
    PcentBias = mean((seloverS1box-ate)/ate)*100,
     MSE=mean((seloverS1box - ate)^2))

seloverS1.box.plot<-ggplot(seloverS1box.dat, aes(x=catincover, y=seloverS1box)) + geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", shape=5, size=4) + facet_grid(~ labincover ) +
    coord_flip() + geom_hline(yintercept=ate, colour="red")  + scale_x_discrete(limits=c("trueMod", "ModMis", "trueNone", "NoneMis")) + ylab("ATE") +xlab("") + ggtitle("ATE, Scenario 1, Tx Model")+ theme( title=element_text(size=9))

sumcover.seloverS1<- ddply(seloverS1cover.dat, .(labincover, catincover), summarise, 
  meancov = mean(seloverS1cover)*100)
sumvar.seloverS1<- ddply(seloverS1var.dat, .(labincover, catincover), summarise, 
  meanvar = mean(seloverS1var))

seloverS1.cov.plot<-ggplot(data=seloverS1cover.dat, aes(x=catincover, y=seloverS1cover))  +facet_grid(~ labincover ) +
    geom_bar( stat="identity", fill="lightblue") +
scale_x_discrete(limits=c("trueMod", "ModMis", "trueNone", "NoneMis")) + geom_text(data=sumcover.seloverS1, aes(x=catincover, y=meancov*100, label=meancov*100), vjust=0, size=3) +
   ggtitle("Coverage, Scenario 1, Tx Model") + xlab("") + ylab("") + theme(axis.text.x=element_text(size=7), title=element_text(size=9))

grid.arrange(seloverS1.box.plot, seloverS1.cov.plot, nrow=1)

sum.seloverS1
sumcover.seloverS1
sumvar.seloverS1


a<-sum.txoverS1[c(2,4,6,8,10,12),-3]
b<-sumcover.txoverS1[c(2,4,6,8,10,12),]
c<-sumvar.txoverS1[c(2,4,6,8,10,12),]

d<-cbind(a[,c(1:3)], c[,3], b[,3], a[,4])
e<-rbind(d[4,], d[2,] , d[6,], d[3,], d[1,] , d[5,])
f<-t(e)
f1<-apply(f[-c(1,2),], c(1,2), function(x) as.numeric(x))

a<-sum.seloverS1[c(2,4,6,8,10,12),-3]
b<-sumcover.seloverS1[c(2,4,6,8,10,12),]
c<-sumvar.seloverS1[c(2,4,6,8,10,12),]

d<-cbind(a[,c(1:3)], c[,3], b[,3], a[,4])
e<-rbind(d[4,], d[2,] , d[6,], d[3,], d[1,] , d[5,])
f<-t(e)
g<-apply(f[-c(1,2),], c(1,2), function(x) as.numeric(x))

g1<-rbind(f1, g)
h<-as.data.frame(g1)
rownames(h)<-c("%BiasT", "VarT", "CovT", "MSET", "%BiasS", "VarS", "CovS", "MSES")
colnames(h)<-c("IPW", "DRWLS", "TMLE", "IPW", "DRWLS", "TMLE")
matdat<- matrix(c(0, rep(1,6), 0, rep(3, 6), 0, rep(1,6), 0, rep(3,6),
0, rep(1,6), 0, rep(3, 6), 0, rep(1,6), 0, rep(3,6)), nrow=8, ncol=7, byrow=TRUE)
print(xtable(h, digits=matdat, include.rownames=TRUE, include.colnames=TRUE, type="latex"))
