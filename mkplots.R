### First make plots when the models are correctly specified
############################################################
library(ggplot2)
library(plyr)
library(gridExtra)
library(xtable)

setwd("/Users/kararudolph/Documents/PhD/NIMH/NCSA/cortisol/simulationproject/Output/3level")
#setwd("/Users/kararudolph/Documents/PhD/NIMH/NCSA/cortisol/simulationproject/Rcode/3level")
## Scenario 1

dats1<-read.csv("EffectsS1MMtrue.csv", header=TRUE, stringsAsFactors=FALSE)
covers1<-read.csv("CoverageATES1MMtrue.csv", header=TRUE, stringsAsFactors=FALSE)
vars1<-read.csv("VarS1MMtrue.csv", header=TRUE, stringsAsFactors=FALSE)

tmleS1box<-c(dats1$Naive, dats1$trueIPTW, dats1$trueIPSW, dats1$trueIPTSvyW, dats1$trueIPW, dats1$trueTMLE2, dats1$trueTMLE3)
tmleS1cover<-c(covers1$Naive, covers1$trueIPTW, covers1$trueIPSW, covers1$trueIPTSvyW, covers1$trueIPW, covers1$trueTMLE2, covers1$trueTMLE3)
tmleS1var<-c(vars1$Naive, vars1$trueIPTW, vars1$trueIPSW, vars1$trueIPTSvyW, vars1$trueIPW, vars1$trueTMLE2, vars1$trueTMLE3)

lab1<-c(rep("Naive", 1000), rep("trueIPTW", 1000), rep("trueIPSW", 1000), rep("trueIPTSvyW", 1000),  rep("trueIPW", 1000), rep("trueDRWLS", 1000),rep("trueTMLE", 1000) )

tmleS1box.dat<-data.frame(tmleS1box, lab1)
tmleS1cover.dat<-data.frame(tmleS1cover, lab1)
tmleS1var.dat<-data.frame(tmleS1var, lab1)

tmleS1.box.plot<-ggplot(tmleS1box.dat, aes(x=lab1, y=tmleS1box)) + geom_boxplot() +
scale_x_discrete(limits=c("trueTMLE", "trueDRWLS", "trueIPW", "trueIPTSvyW", "trueIPTW", "trueIPSW",  "Naive")) +
    stat_summary(fun.y=mean, geom="point", shape=5, size=4) + 
    coord_flip() + geom_hline(yintercept=ate, colour="red")  + ylab("ATE") +xlab("") + ggtitle("ATE, Scenario 1")+ theme( title=element_text(size=9))

sumcover.tmleS1<- ddply(tmleS1cover.dat, .(lab1), summarise, 
  meancov = mean(tmleS1cover)*100)

tmleS1.cov.plot<-ggplot(data=tmleS1cover.dat, aes(x=lab1, y=tmleS1cover))  +
    geom_bar( stat="identity", fill="lightblue") +
scale_x_discrete(limits=c("trueTMLE", "trueDRWLS", "trueIPW", "trueIPTSvyW", "trueIPTW", "trueIPSW",  "Naive"))+ 
  geom_text(data=sumcover.tmleS1, aes(x=lab1, y=meancov, label=meancov), vjust=0, size=3) + ggtitle("Coverage, Scenario 1") + xlab("") + ylab("") + theme(axis.text.x=element_text(size=7), title=element_text(size=9))

sum.tmleS1 <- ddply(tmleS1box.dat, .(lab1), summarise, 
  AbsBias = mean(tmleS1box-ate),
    PcentBias = (mean((tmleS1box-ate)/ate))*100,
     MSE=mean((tmleS1box - ate)^2) )

sum.varS1 <- ddply(tmleS1var.dat, .(lab1), summarise, 
  meanvar=mean(tmleS1var) )

sum.tmleS1
sumcover.tmleS1
sum.varS1

grid.arrange(tmleS1.box.plot, tmleS1.cov.plot, nrow=1)

a<-cbind(sum.tmleS1[c(2,6,7), -2], sumcover.tmleS1[c(2,6,7), 2], sum.varS1[c(2,6,7),2] )
k<-c(a[2,2], a[2,5], a[2,4], a[2,3], a[3,2], a[3,5], a[3,4], a[3,3], a[1,2], a[1,5], a[1,4], a[1,3])
## Scenario 2

dats2<-read.csv("EffectsS2MMtrue.csv", header=TRUE, stringsAsFactors=FALSE)
covers2<-read.csv("CoverageATES2MMtrue.csv", header=TRUE, stringsAsFactors=FALSE)
vars2<-read.csv("VarS2MMtrue.csv", header=TRUE, stringsAsFactors=FALSE)

tmleS2box<-c(dats2$Naive, dats2$trueIPTW, dats2$trueIPSW, dats2$trueIPTSvyW, dats2$trueIPW, dats2$trueTMLE2, dats2$trueTMLE3)
tmleS2cover<-c(covers2$Naive, covers2$trueIPTW, covers2$trueIPSW, covers2$trueIPTSvyW, covers2$trueIPW, covers2$trueTMLE2, covers2$trueTMLE3)
tmleS2var<-c(vars2$Naive, vars2$trueIPTW, vars2$trueIPSW, vars2$trueIPTSvyW, vars2$trueIPW, vars2$trueTMLE2, vars2$trueTMLE3)

lab1<-c(rep("Naive", 1000), rep("trueIPTW", 1000), rep("trueIPSW", 1000), rep("trueIPTSvyW", 1000),  rep("trueIPW", 1000), rep("trueDRWLS", 1000),rep("trueTMLE", 1000) )

tmleS2box.dat<-data.frame(tmleS2box, lab1)
tmleS2cover.dat<-data.frame(tmleS2cover, lab1)
tmleS2var.dat<-data.frame(tmleS2var, lab1)

tmleS2.box.plot<-ggplot(tmleS2box.dat, aes(x=lab1, y=tmleS2box)) + geom_boxplot() +
scale_x_discrete(limits=c("trueTMLE", "trueDRWLS", "trueIPW", "trueIPTSvyW", "trueIPTW", "trueIPSW",  "Naive")) +
    stat_summary(fun.y=mean, geom="point", shape=5, size=4) + 
    coord_flip() + geom_hline(yintercept=ate, colour="red")  + ylab("ATE") +xlab("") + ggtitle("ATE, Scenario 2")+ theme( title=element_text(size=9))

sumcover.tmleS2<- ddply(tmleS2cover.dat, .(lab1), summarise, 
  meancov = mean(tmleS2cover)*100)
tmleS2.cov.plot<-ggplot(data=tmleS2cover.dat, aes(x=lab1, y=tmleS2cover))  +
    geom_bar( stat="identity", fill="lightblue") +
scale_x_discrete(limits=c("trueTMLE", "trueDRWLS", "trueIPW", "trueIPTSvyW", "trueIPTW", "trueIPSW",  "Naive"))+ 
  geom_text(data=sumcover.tmleS2, aes(x=lab1, y=meancov, label=meancov), vjust=0, size=3) + ggtitle("Coverage, Scenario 2") + xlab("") + ylab("") + theme(axis.text.x=element_text(size=7), title=element_text(size=9))
sum.tmleS2 <- ddply(tmleS2box.dat, .(lab1), summarise, 
  AbsBias = mean(tmleS2box-ate),
    PcentBias = mean((tmleS2box-ate)/ate)*100,
     MSE=mean((tmleS2box - ate)^2))
sum.varS2<- ddply(tmleS2var.dat, .(lab1), summarise, 
  meanvar = mean(tmleS2var))

grid.arrange(tmleS2.box.plot, tmleS2.cov.plot, nrow=1)

sum.tmleS2
sumcover.tmleS2
sum.varS2

a<-cbind(sum.tmleS1[c(2,6,7), -2], sumcover.tmleS1[c(2,6,7), 2], sum.varS1[c(2,6,7),2] )
d<-c(a[2,2], a[2,5], a[2,4], a[2,3], a[3,2], a[3,5], a[3,4], a[3,3], a[1,2], a[1,5], a[1,4], a[1,3])

############################################################
### Make plots when the models are misspecified
############################################################

## Scenario 1
dat.S1mistx<-read.csv("EffectsS1MMTx.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1mistx<-read.csv("CoverageATES1MMTx.csv", header=TRUE, stringsAsFactors=FALSE)
var.S1mistx<-read.csv("VarS1MMTx.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1missel<-read.csv("EffectsS1MMSel.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1missel<-read.csv("CoverageATES1MMSel.csv", header=TRUE, stringsAsFactors=FALSE)
var.S1missel<-read.csv("VarS1MMSel.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1misout<-read.csv("EffectsS1MMOut.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1misout<-read.csv("CoverageATES1MMOut.csv", header=TRUE, stringsAsFactors=FALSE)
var.S1misout<-read.csv("VarS1MMOut.csv", header=TRUE, stringsAsFactors=FALSE)

tmleS1box<-c( dat.S1mistx$IPWT1, dat.S1mistx$TMLE2T1, dat.S1mistx$TMLE3T1,
dat.S1mistx$IPWT2, dat.S1mistx$TMLE2T2, dat.S1mistx$TMLE3T2, 
dat.S1missel$IPWS1, dat.S1missel$TMLE2S1, dat.S1missel$TMLE3S1,
dat.S1missel$IPWS2, dat.S1missel$TMLE2S2, dat.S1missel$TMLE3S2, 
dat.S1misout$TMLE2O1, dat.S1misout$TMLE3O1,
dat.S1misout$TMLE2O2, dat.S1misout$TMLE3O2, 
dat.S1misout$TMLE2O3, dat.S1misout$TMLE3O3 )

tmleS1cover<-c( cover.S1mistx$IPWT1, cover.S1mistx$TMLE2T1, cover.S1mistx$TMLE3T1,
cover.S1mistx$IPWT2, cover.S1mistx$TMLE2T2, cover.S1mistx$TMLE3T2, 
cover.S1missel$IPWS1, cover.S1missel$TMLE2S1, cover.S1missel$TMLE3S1,
cover.S1missel$IPWS2, cover.S1missel$TMLE2S2, cover.S1missel$TMLE3S2,  
cover.S1misout$TMLE2O1, cover.S1misout$TMLE3O1,
cover.S1misout$TMLE2O2, cover.S1misout$TMLE3O2, 
cover.S1misout$TMLE2O3, cover.S1misout$TMLE3O3 )

tmleS1var<-c( var.S1mistx$IPWT1, var.S1mistx$TMLE2T1, var.S1mistx$TMLE3T1,
var.S1mistx$IPWT2, var.S1mistx$TMLE2T2, var.S1mistx$TMLE3T2, 
var.S1missel$IPWS1, var.S1missel$TMLE2S1, var.S1missel$TMLE3S1,
var.S1missel$IPWS2, var.S1missel$TMLE2S2, var.S1missel$TMLE3S2, 
var.S1misout$TMLE2O1, var.S1misout$TMLE3O1,
var.S1misout$TMLE2O2, var.S1misout$TMLE3O2, 
var.S1misout$TMLE2O3, var.S1misout$TMLE3O3 )

labS2<-c(rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), 
  rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), 
   rep("DRWLS", 1000),  rep("TMLE", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000))

catS2<-c(rep("modtx", 3000), rep("majtx", 3000), rep("modsel", 3000), rep("majsel", 3000), rep("modout", 2000), rep("majout1", 2000),rep("majout2", 2000) )

tmleS1box.dat<-data.frame(tmleS1box, labS2, catS2)
tmleS1cover.dat<-data.frame(tmleS1cover, labS2, catS2)
tmleS1var.dat<-data.frame(tmleS1var, labS2, catS2)

sum.tmleS1 <- ddply(tmleS1box.dat, .(labS2, catS2), summarise, 
  AbsBias = mean(tmleS1box-ate),
    PcentBias = mean((tmleS1box-ate)/ate)*100,
     MSE=mean((tmleS1box - ate)^2))

tmleS1.box.plot<-ggplot(tmleS1box.dat, aes(x=catS2, y=tmleS1box)) + geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", shape=5, size=4) + facet_grid(~ labS2 ) +
    coord_flip() + geom_hline(yintercept=ate, colour="red")  + scale_x_discrete(limits=c("modtx", "majtx", "modsel", "majsel","modout", "majout1", "majout2")) + ylab("ATE") +xlab("") + ggtitle("ATE, Scenario 1")+ theme( title=element_text(size=9))

sumcover.tmleS1<- ddply(tmleS1cover.dat, .(labS2, catS2), summarise, 
  meancov = mean(tmleS1cover)*100)
sumvar.tmleS1<- ddply(tmleS1var.dat, .(labS2, catS2), summarise, 
  meanvar = mean(tmleS1var))

tmleS1.cov.plot<-ggplot(data=tmleS1cover.dat, aes(x=catS2, y=tmleS1cover))  +facet_grid(~ labS2 ) +
    geom_bar( stat="identity", fill="lightblue") +
scale_x_discrete(limits=c("modtx", "majtx", "modsel", "majsel","modout", "majout1", "majout2")) + 
   ggtitle("Coverage, Scenario 1") + xlab("") + ylab("") + theme(axis.text.x=element_text(size=7), title=element_text(size=9))

grid.arrange(tmleS1.box.plot, tmleS1.cov.plot, nrow=1)

sumcover.tmleS1
sumvar.tmleS1

a<-cbind(sum.tmleS1[, c(1,2,4)], sumvar.tmleS1[,3], sumcover.tmleS1[, 3], sum.tmleS1[,5] )
b<-rbind(a[c(1:7),], c("IPW", "majout1", 0,0,0,0), c("IPW", "majout2", 0,0,0,0) , a[c(8:9),] , c("IPW", "modout", 0,0,0,0), a[c(10:18),])
b1<-cbind(b[c(8:14),],b[1:7,],b[15:21,])
b2<-b1[,-c(1,7,8, 13,14)]
b3<-rbind(b2[7,],b2[4,], b2[5,], b2[1:2,], b2[6,], b2[3,])
#b3<-rbind(b2["modtx",],b2["majtx",], b2["modout",],b2["majout1",], b2["majout2",], b2["modsel",],b2["majsel",] )
e<-b3[,-1]

## Scenario 2

dat.S2mistx<-read.csv("EffectsS2MMTx.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2mistx<-read.csv("CoverateATES2MMTx.csv", header=TRUE, stringsAsFactors=FALSE)
var.S2mistx<-read.csv("VarS2MMTx.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S2missel<-read.csv("EffectsS2MMSel.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2missel<-read.csv("CoverageATES2MMSel.csv", header=TRUE, stringsAsFactors=FALSE)
var.S2missel<-read.csv("VarS2MMSel.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S2misout<-read.csv("EffectsS2MMOut.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2misout<-read.csv("CoverageATES2MMOut.csv", header=TRUE, stringsAsFactors=FALSE)
var.S2misout<-read.csv("VarS2MMOut.csv", header=TRUE, stringsAsFactors=FALSE)

tmleS2box<-c( dat.S2mistx$IPWT1, dat.S2mistx$TMLE2T1, dat.S2mistx$TMLE3T1,
dat.S2mistx$IPWT2, dat.S2mistx$TMLE2T2, dat.S2mistx$TMLE3T2, 
dat.S2missel$IPWS1, dat.S2missel$TMLE2S1, dat.S2missel$TMLE3S1,
dat.S2missel$IPWS2, dat.S2missel$TMLE2S2, dat.S2missel$TMLE3S2, 
dat.S2misout$TMLE2O1, dat.S2misout$TMLE3O1,
dat.S2misout$TMLE2O2, dat.S2misout$TMLE3O2, 
dat.S2misout$TMLE2O3, dat.S2misout$TMLE3O3 )

tmleS2cover<-c( cover.S2mistx$IPWT1, cover.S2mistx$TMLE2T1, cover.S2mistx$TMLE3T1,
cover.S2mistx$IPWT2, cover.S2mistx$TMLE2T2, cover.S2mistx$TMLE3T2, 
cover.S2missel$IPWS1, cover.S2missel$TMLE2S1, cover.S2missel$TMLE3S1,
cover.S2missel$IPWS2, cover.S2missel$TMLE2S2, cover.S2missel$TMLE3S2,  
cover.S2misout$TMLE2O1, cover.S2misout$TMLE3O1,
cover.S2misout$TMLE2O2, cover.S2misout$TMLE3O2, 
cover.S2misout$TMLE2O3, cover.S2misout$TMLE3O3 )

tmleS2var<-c( var.S2mistx$IPWT1, var.S2mistx$TMLE2T1, var.S2mistx$TMLE3T1,
var.S2mistx$IPWT2, var.S2mistx$TMLE2T2, var.S2mistx$TMLE3T2, 
var.S2missel$IPWS1, var.S2missel$TMLE2S1, var.S2missel$TMLE3S1,
var.S2missel$IPWS2, var.S2missel$TMLE2S2, var.S2missel$TMLE3S2, 
var.S2misout$TMLE2O1, var.S2misout$TMLE3O1,
var.S2misout$TMLE2O2, var.S2misout$TMLE3O2, 
var.S2misout$TMLE2O3, var.S2misout$TMLE3O3 )

labS2<-c(rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), 
  rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), 
   rep("DRWLS", 1000),  rep("TMLE", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000))

catS2<-c(rep("modtx", 3000), rep("majtx", 3000), rep("modsel", 3000), rep("majsel", 3000), rep("modout", 2000), rep("majout1", 2000),rep("majout2", 2000) )

tmleS2box.dat<-data.frame(tmleS2box, labS2, catS2)
tmleS2cover.dat<-data.frame(tmleS2cover, labS2, catS2)
tmleS2var.dat<-data.frame(tmleS2var, labS2, catS2)

sum.tmleS2 <- ddply(tmleS2box.dat, .(labS2, catS2), summarise, 
  AbsBias = mean(tmleS2box-ate),
    PcentBias = mean((tmleS2box-ate)/ate)*100,
     MSE=mean((tmleS2box - ate)^2))

tmleS2.box.plot<-ggplot(tmleS2box.dat, aes(x=catS2, y=tmleS2box)) + geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", shape=5, size=4) + facet_grid(~ labS2 ) +
    coord_flip() + geom_hline(yintercept=ate, colour="red")  + scale_x_discrete(limits=c("modtx", "majtx", "modsel", "majsel","modout", "majout1", "majout2")) + ylab("ATE") +xlab("") + ggtitle("ATE, Scenario 2")+ theme( title=element_text(size=9))

sumcover.tmleS2<- ddply(tmleS2cover.dat, .(labS2, catS2), summarise, 
  meancov = mean(tmleS2cover)*100)
sumvar.tmleS2<- ddply(tmleS2var.dat, .(labS2, catS2), summarise, 
  meanvar = mean(tmleS2var))

tmleS2.cov.plot<-ggplot(data=tmleS2cover.dat, aes(x=catS2, y=tmleS2cover))  +facet_grid(~ labS2 ) +
    geom_bar( stat="identity", fill="lightblue") +
scale_x_discrete(limits=c("modtx", "majtx", "modsel", "majsel","modout", "majout1", "majout2")) + 
   ggtitle("Coverage, Scenario 2") + xlab("") + ylab("") + theme(axis.text.x=element_text(size=7), title=element_text(size=9))

grid.arrange(tmleS2.box.plot, tmleS2.cov.plot, nrow=1)

sum.tmleS2
sumcover.tmleS2
sumvar.tmleS2

a<-cbind(sum.tmleS2[, c(1,2,4)], sumvar.tmleS2[,3], sumcover.tmleS2[, 3], sum.tmleS2[,5] )
b<-rbind(a[c(1:7),], c("IPW", "majout1", 0,0,0,0), c("IPW", "majout2", 0,0,0,0) , a[c(8:9),] , c("IPW", "modout", 0,0,0,0), a[c(10:18),])
b1<-cbind(b[c(8:14),],b[1:7,],b[15:21,])
b2<-b1[,-c(1,7,8, 13,14)]
b3<-rbind(b2[7,],b2[4,], b2[5,], b2[1:2,], b2[6,], b2[3,])
#b3<-rbind(b2["modtx",],b2["majtx",], b2["modout",],b2["majout1",], b2["majout2",], b2["modsel",],b2["majsel",] )
f<-b3[,-1]
names(f)<-names(e)
g<-rbind(k,e,d,f)
g1<-apply(g, c(1,2), function(x) as.numeric(x))
h<-as.data.frame(g1)
colnames(h)<-c("%Bias", "Var", "Cov", "MSE", "%Bias", "Var", "Cov", "MSE", "%Bias", "Var", "Cov", "MSE")
rownames(h)<-c("S1TrueMod", "S1ModMisTx", "S1MajMisTx", "S1ModMisOut1", "S1MajMisOut2", "S1MajMisOut3", "S1ModMisSel", "S1MajMisSel","S2TrueMod", "S2ModMisTx", "S2MajMisTx", "S2ModMisOut1", "S2MajMisOut2", "S2MajMisOut3", "S2ModMisSel", "S2MajMisSel" )

print(xtable(h, digits=c(0,1,3,1,3,1,3,1,3,1,3,1,3)), include.rownames=TRUE, include.colnames=TRUE, type="latex")


### Sensitivity to incorrect correction of selection or treatment non-random mechanism under Scenario 1
## Under-correct. 

# Treatment. Treatment mechanism is ~ poly(z2,2). No correction
datTxunder<- read.csv("EffectsS1MMIncorATx.csv", header=TRUE, stringsAsFactors=FALSE)
coverTxunder<- read.csv("CoverageATES1MMIncorATx.csv", header=TRUE, stringsAsFactors=FALSE)
varTxunder<- read.csv("VarS1MMIncorATx.csv", header=TRUE, stringsAsFactors=FALSE)

# Selection. Sub-sample selection mechanism is ~ poly(z1,2) + z2. No correction
datSelunder<- read.csv("EffectsS1MMIncorA.csv", header=TRUE, stringsAsFactors=FALSE)
coverSelunder<- read.csv("CoverageATES1MMIncorA.csv", header=TRUE, stringsAsFactors=FALSE)
varSelunder<- read.csv("VarS1MMIncorA.csv", header=TRUE, stringsAsFactors=FALSE)

underS1box<-c(dats1$trueIPSW, datTxunder$TMLE2, datTxunder$TMLE3, dats1$trueIPTSvyW, datSelunder$TMLE2, datSelunder$TMLE3)
underS1cover<-c(covers1$trueIPSW, coverTxunder$TMLE2, coverTxunder$TMLE3, covers1$trueIPTSvyW, coverSelunder$TMLE2, coverSelunder$TMLE3)
underS1var<-c(vars1$trueIPSW, varTxunder$TMLE2, varTxunder$TMLE3, vars1$trueIPTSvyW, varSelunder$TMLE2, varSelunder$TMLE3)

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
datTxover1<- read.csv("EffectsS1MMIncorATx.csv", header=TRUE, stringsAsFactors=FALSE)
coverTxover1<- read.csv("CoverageATES1MMIncorATx.csv", header=TRUE, stringsAsFactors=FALSE)
varTxover1<- read.csv("VarS1MMIncorATx.csv", header=TRUE, stringsAsFactors=FALSE)
# Treatment 2. Treatment mechanism is ~ 1. Correction mechanism is: ~poly(z2,2)
datTxover2<- read.csv("EffectsS1MMIncorATx.csv", header=TRUE, stringsAsFactors=FALSE)
coverTxover2<- read.csv("CoverageATES1MMIncorATx.csv", header=TRUE, stringsAsFactors=FALSE)
varTxover2<- read.csv("VarS1MMIncorATx.csv", header=TRUE, stringsAsFactors=FALSE)

# Selection 1. Sub-sample selection mechanism is ~ poly(z1,2). Correction mechanism is:  ~ poly(z1,2) + z2
datSelover1<- read.csv("EffectsS1MMIncorA.csv", header=TRUE, stringsAsFactors=FALSE)
coverSelover1<- read.csv("CoverageATES1MMIncorA.csv", header=TRUE, stringsAsFactors=FALSE)
varSelover1<- read.csv("VarS1MMIncorA.csv", header=TRUE, stringsAsFactors=FALSE)
# Selection 2. Sub-sample selection mechanism is ~ 1. Correction mechanism is:  ~ poly(z1,2) + z2
datSelover2<- read.csv("EffectsS1MMIncorA.csv", header=TRUE, stringsAsFactors=FALSE)
coverSelover2<- read.csv("CoverageATES1MMIncorA.csv", header=TRUE, stringsAsFactors=FALSE)
varSelover2<- read.csv("VarS1MMIncorA.csv", header=TRUE, stringsAsFactors=FALSE)

txoverS1box<-c(datTxover1$IPW, datTxover1$TMLE2, datTxover1$TMLE3,
  datTxover1$IPWover, datTxover1$TMLE2over, datTxover1$TMLE3over,
datTxover2$IPW, datTxover2$TMLE2, datTxover2$TMLE3,
datTxover2$IPWover, datTxover2$TMLE2over, datTxover2$TMLE3over)

txoverS1cover<-c(coverTxover1$IPW, coverTxover1$TMLE2, coverTxover1$TMLE3,
  coverTxover1$IPWover, coverTxover1$TMLE2over, coverTxover1$TMLE3over,
coverTxover2$IPW, coverTxover2$TMLE2, coverTxover2$TMLE3,
coverTxover2$IPWover, coverTxover2$TMLE2over, coverTxover2$TMLE3over)

txoverS1var<-c(varTxover1$IPW, varTxover1$TMLE2, varTxover1$TMLE3,
  varTxover1$IPWover, varTxover1$TMLE2over, varTxover1$TMLE3over,
varTxover2$IPW, varTxover2$TMLE2, varTxover2$TMLE3,
varTxover2$IPWover, varTxover2$TMLE2over, varTxover2$TMLE3over)

seloverS1box<-c(datSelover1$IPW, datSelover1$TMLE2, datSelover1$TMLE3,
  datSelover1$IPWover, datSelover1$TMLE2over, datSelover1$TMLE3over,
datSelover2$IPW, datSelover2$TMLE2, datSelover2$TMLE3,
datSelover2$IPWover, datSelover2$TMLE2over, datSelover2$TMLE3over)

seloverS1cover<-c(coverSelover1$IPW, coverSelover1$TMLE2, coverSelover1$TMLE3,
  coverSelover1$IPWover, coverSelover1$TMLE2over, coverSelover1$TMLE3over,
coverSelover2$IPW, coverSelover2$TMLE2, coverSelover2$TMLE3,
coverSelover2$IPWover, coverSelover2$TMLE2over, coverSelover2$TMLE3over)

seloverS1var<-c(varSelover1$IPW, varSelover1$TMLE2, varSelover1$TMLE3,
  varSelover1$IPWover, varSelover1$TMLE2over, varSelover1$TMLE3over,
varSelover2$IPW, varSelover2$TMLE2, varSelover2$TMLE3,
varSelover2$IPWover, varSelover2$TMLE2over, varSelover2$TMLE3over)

labincover<-c( rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000), rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000) ,
  rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000) 
  rep("IPW", 1000), rep("DRWLS", 1000),  rep("TMLE", 1000) )
catincover<-c(rep("trueMod", 3000), rep("ModMis", 3000),  rep("trueNone", 3000), rep("NoneMis", 3000),)

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
  meancov = mean(txoverS1cover))
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
  meancov = mean(seloverS1cover))
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

