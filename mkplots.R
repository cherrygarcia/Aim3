setwd("/Users/kararudolph/Documents/PhD/NIMH/NCSA/cortisol/simulationproject/Output/3level")

dat<-read.csv("EffectsTrueModelScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover<-read.csv("CoverageATETrueModelScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.btwts<-read.csv("EffectsTrueModelScenOne3levels_bswts.csv", header=TRUE, stringsAsFactors=FALSE)
cover.btwts<-read.csv("CoverageATETrueModelScenOne3levels_bswts.csv", header=TRUE, stringsAsFactors=FALSE)
dat.flucsvy2<-read.csv("EffectsScenOne3levelsflucsvy2.csv", header=TRUE, stringsAsFactors=FALSE)
cover.flucsvy2<-read.csv("CoverageATEScenOne3levelsflucsvy2.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1true2<-read.csv("EffectsScenOne3levelsflucsvymatrixmod.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1true2<-read.csv("CoverageATEScenOne3levelsflucsvymatrixmod.csv", header=TRUE, stringsAsFactors=FALSE)

tmleS1box<-c(dat$Naive, dat$Tx_Sub, dat$SvySub, dat.btwts$BootIPW, dat.flucsvy2$TMLESvy, dat.S1true2$TMLESvy2, dat.S1true2$TMLESvyFluc)
tmleS1cover<-c(cover$Naive, cover$Tx_Sub, cover$SvySub, cover.btwts$BootIPW, cover.flucsvy2$TMLESvy, cover.S1true2$TMLESvy2, cover.S1true2$TMLESvyFluc)

lab1<-c(rep("Naive", 250), rep("IPTW", 250), rep("IPSW", 250), rep("IPTSW", 250),  rep("TMLE1", 250),rep("TMLE2", 250),rep("TMLE3", 250) )

tmleS1box.dat<-data.frame(tmleS1box, lab1)
tmleS1cover.dat<-data.frame(tmleS1cover, lab1)

tmleS1.box.plot<-ggplot(tmleS1box.dat, aes(x=lab1, y=tmleS1box)) + geom_boxplot() +
scale_x_discrete(limits=c("TMLE3", "TMLE2", "TMLE1", "IPTSW", "IPTW", "IPSW",  "Naive")) +
    stat_summary(fun.y=mean, geom="point", shape=5, size=4) + 
    coord_flip() + geom_hline(yintercept=ate, colour="red")  + ylab("ATE") +xlab("") + ggtitle("ATE, Scenario 1")+ theme( title=element_text(size=9))

sumcover.tmleS1<- ddply(tmleS1cover.dat, .(lab1), summarise, 
  meancov = mean(tmleS1cover))
tmleS1.cov.plot<-ggplot(data=tmleS1cover.dat, aes(x=lab1, y=tmleS1cover))  +
    geom_bar( stat="identity", fill="lightblue") +
scale_x_discrete(limits=c("TMLE3", "TMLE2", "TMLE1", "IPTSW", "IPTW", "IPSW",  "Naive"))+ 
  geom_text(data=sumcover.tmleS1, aes(x=lab1, y=meancov*100, label=meancov*100), vjust=0, size=3) + ggtitle("Coverage, Scenario 1") + xlab("") + ylab("") + theme(axis.text.x=element_text(size=7), title=element_text(size=9))

grid.arrange(tmleS1.box.plot, tmleS1.cov.plot, nrow=1)

## Scenario 2

dats2<-read.csv("EffectsTrueModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
covers2<-read.csv("CoverageATETrueModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.flucsvy2<-read.csv("EffectsScenTwo3levelsTruesvy.csv", header=TRUE, stringsAsFactors=FALSE)
cover.flucsvy2<-read.csv("CoverageATEScenTwo3levelsTruesvy.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S2true2<-read.csv("EffectsScenTwo3levelsflucsvymatrixmod.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2true2<-read.csv("CoverageATEScenTwo3levelsflucsvymatrixmod.csv", header=TRUE, stringsAsFactors=FALSE)

tmleS2box<-c(dats2$Naive, dats2$Tx_Sub, dats2$SvySub, dats2$TxSvySub, dat.flucsvy2$TMLESvy, dat.S2true2$TMLESvy2, dat.S2true2$TMLESvyFluc)
tmleS2cover<-c(covers2$Naive, covers2$Tx_Sub, covers2$SvySub, covers2$TxSvySub, cover.flucsvy2$TMLESvy, cover.S2true2$TMLESvy2, cover.S2true2$TMLESvyFluc)

lab1<-c(rep("Naive", 250), rep("IPTW", 250), rep("IPSW", 250), rep("IPTSW", 250),  rep("TMLE1", 250),rep("TMLE2", 250),rep("TMLE3", 250) )

tmleS2box.dat<-data.frame(tmleS2box, lab1)
tmleS2cover.dat<-data.frame(tmleS2cover, lab1)

tmleS2.box.plot<-ggplot(tmleS2box.dat, aes(x=lab1, y=tmleS2box)) + geom_boxplot() +
scale_x_discrete(limits=c("TMLE3", "TMLE2", "TMLE1", "IPTSW", "IPTW", "IPSW",  "Naive")) +
    stat_summary(fun.y=mean, geom="point", shape=5, size=4) + 
    coord_flip() + geom_hline(yintercept=ate, colour="red")  + ylab("ATE") +xlab("") + ggtitle("ATE, Scenario 1")+ theme( title=element_text(size=9))

sumcover.tmleS2<- ddply(tmleS2cover.dat, .(lab1), summarise, 
  meancov = mean(tmleS2cover))
tmleS2.cov.plot<-ggplot(data=tmleS2cover.dat, aes(x=lab1, y=tmleS2cover))  +
    geom_bar( stat="identity", fill="lightblue") +
scale_x_discrete(limits=c("TMLE3", "TMLE2", "TMLE1", "IPTSW", "IPTW", "IPSW",  "Naive"))+ 
  geom_text(data=sumcover.tmleS2, aes(x=lab1, y=meancov*100, label=meancov*100), vjust=0, size=3) + ggtitle("Coverage, Scenario 1") + xlab("") + ylab("") + theme(axis.text.x=element_text(size=7), title=element_text(size=9))

grid.arrange(tmleS2.box.plot, tmleS2.cov.plot, nrow=1)

### Model Misspecification
## Scenario 1
dat.S1mistxmin<-read.csv("EffectsMisspecTx1ModelScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1mistxmin<-read.csv("CoverageATEMisspecTx1ModelScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1mistxmaj<-read.csv("EffectsMisspecTx2ModelScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1mistxmaj<-read.csv("CoverageATEMisspecTx2ModelScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1misselmin<-read.csv("EffectsMisSel1ScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1misselmin<-read.csv("CoverageATEMisSel1ScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1misselmaj<-read.csv("EffectsMisSel2ScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1misselmaj<-read.csv("CoverageATEMisSel2ScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1misoutmin<-read.csv("EffectsMisspecOut1ModelScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1misoutmin<-read.csv("CoverageATEMisspecOut1ScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1misoutmaj1<-read.csv("EffectsMisspecOut2ModelScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1misoutmaj1<-read.csv("CoverageATEMisspecOut2ScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1misoutmaj2<-read.csv("EffectsMisspecOut3ModelScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1misoutmaj2<-read.csv("CoverageATEMisspecOut3ScenOne3levels.csv", header=TRUE, stringsAsFactors=FALSE)

dat.S1mistx2<-read.csv("EffectsScenOne3levelsflucsvyMMTxMis.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1mistx2<-read.csv("CoverageATEScenOne3levelsflucsvyMMTxMis.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1missel2<-read.csv("EffectsScenOne3levelsflucsvyMMSelMis.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1missel2<-read.csv("CoverageATEScenOne3levelsflucsvyMMSelMis.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1misout2<-read.csv("EffectsScenOne3levelsflucsvyMMOutMis.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1misout2<-read.csv("CoverageATEScenOne3levelsflucsvyMMOutMis.csv", header=TRUE, stringsAsFactors=FALSE)

dat.S1mistx3<-read.csv("EffectsScenOne3levelsTxMissvy.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1mistx3<-read.csv("CoverageATEScenOne3levelsTxMissvy.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1missel3<-read.csv("EffectsScenOne3levelsflucsvy3.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1missel3<-read.csv("CoverageATEScenOne3levelsflucsvy3.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1missel4<-read.csv("EffectsScenOne3levelsflucsvy5.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1missel4<-read.csv("CoverageATEScenOne3levelsflucsvy5.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S1misout4<-read.csv("EffectsScenOne3levelsOutMissvy.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S1misout4<-read.csv("CoverageATEScenOne3levelsOutMissvy.csv", header=TRUE, stringsAsFactors=FALSE)

tmleS1box<-c(
  dat.S1mistxmin$TxSvySub, dat.S1mistx3$TMLESvyT1, dat.S1mistx3$TMLESvy2T1, dat.S1mistx2$TMLESvyFlucT1, 
  dat.S1mistxmaj$TxSvySub, dat.S1mistx3$TMLESvyT2, dat.S1mistx3$TMLESvy2T2, dat.S1mistx2$TMLESvyFlucT2, 
  dat.S1misselmin$TxSvySub, dat.S1missel3$TMLESvyM1, dat.S1missel2$TMLESvy2S1, dat.S1missel2$TMLESvyFlucS1, 
  dat.S1misselmaj$TxSvySub, dat.S1missel4$TMLESvyM2, dat.S1missel4$TMLESvy2M2, dat.S1missel2$TMLESvyFlucS2,
  dat.S1misoutmin$TxSvySub, dat.S1misout4$TMLESvyO1, dat.S1misout4$TMLESvy2O1, dat.S1misout2$TMLESvyFlucO1,
  dat.S1misoutmaj1$TxSvySub, dat.S1misout4$TMLESvyO2, dat.S1misout4$TMLESvy2O2, dat.S1misout2$TMLESvyFlucO2,
  dat.S1misoutmaj2$TxSvySub, dat.S1misout4$TMLESvyO3, dat.S1misout4$TMLESvy2O3, dat.S1misout2$TMLESvyFlucO3)
tmleS1cover<-c(
  cover.S1mistxmin$TxSvySub, cover.S1mistx3$TMLESvyT1, cover.S1mistx3$TMLESvy2T1, cover.S1mistx2$TMLESvyFlucT1, 
  cover.S1mistxmaj$TxSvySub, cover.S1mistx3$TMLESvyT2, cover.S1mistx3$TMLESvy2T2, cover.S1mistx2$TMLESvyFlucT2, 
  cover.S1misselmin$TxSvySub, cover.S1missel3$TMLESvyM1, cover.S1missel2$TMLESvy2S1, cover.S1missel2$TMLESvyFlucS1, 
  cover.S1misselmaj$TxSvySub, cover.S1missel4$TMLESvyM2, cover.S1missel4$TMLESvy2M2, cover.S1missel2$TMLESvyFlucS2,
  cover.S1misoutmin$TxSvySub, cover.S1misout4$TMLESvyO1, cover.S1misout4$TMLESvy2O1, cover.S1misout2$TMLESvyFlucO1,
  cover.S1misoutmaj1$TxSvySub, cover.S1misout4$TMLESvyO2, cover.S1misout4$TMLESvy2O2, cover.S1misout2$TMLESvyFlucO2,
  cover.S1misoutmaj2$TxSvySub, cover.S1misout4$TMLESvyO3, cover.S1misout4$TMLESvy2O3, cover.S1misout2$TMLESvyFlucO3)
labS2<-c(rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250), rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250),rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250),rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250),rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250),rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250),rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250))
catS2<-c(rep("modtx", 1000), rep("majtx", 1000), rep("modsel", 1000), rep("majsel", 1000), rep("modout", 1000), rep("majout1", 1000),rep("majout2", 1000) )
tmleS1box.dat<-data.frame(tmleS1box, labS2, catS2)
tmleS1cover.dat<-data.frame(tmleS1cover, labS2, catS2)

sum.tmleS1 <- ddply(tmleS1box.dat, .(labS2, catS2), summarise, 
  AbsBias = mean(tmleS1box-ate),
    PcentBias = mean((tmleS1box-ate)/ate),
     MSE=mean((tmleS1box - ate)^2))

tmleS1.box.plot<-ggplot(tmleS1box.dat, aes(x=catS2, y=tmleS1box)) + geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", shape=5, size=4) + facet_grid(~ labS2 ) +
    coord_flip() + geom_hline(yintercept=ate, colour="red")  + scale_x_discrete(limits=c("modtx", "majtx", "modsel", "majsel","modout", "majout1", "majout2")) + ylab("ATE") +xlab("") + ggtitle("ATE, Scenario 1")+ theme( title=element_text(size=9))

sumcover.tmleS1<- ddply(tmleS1cover.dat, .(labS2, catS2), summarise, 
  meancov = mean(tmleS1cover))
tmleS1.cov.plot<-ggplot(data=tmleS1cover.dat, aes(x=catS2, y=tmleS1cover))  +facet_grid(~ labS2 ) +
    geom_bar( stat="identity", fill="lightblue") +
scale_x_discrete(limits=c("modtx", "majtx", "modsel", "majsel","modout", "majout1", "majout2")) + 
   ggtitle("Coverage, Scenario 1") + xlab("") + ylab("") + theme(axis.text.x=element_text(size=7), title=element_text(size=9))

grid.arrange(tmleS1.box.plot, tmleS1.cov.plot, nrow=1)


## Scenario 2
dat.S2mistxmin<-read.csv("EffectsMisspecTx1ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2mistxmin<-read.csv("CoverageATEMisspecTx1ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S2mistxmaj<-read.csv("EffectsMisspecTx2ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2mistxmaj<-read.csv("CoverageATEMisspecTx2ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S2misselmin<-read.csv("EffectsMisspecSel1ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2misselmin<-read.csv("CoverageATEMisspecSel1ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S2misselmaj<-read.csv("EffectsMisspecSel2ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2misselmaj<-read.csv("CoverageATEMisspecSel2ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S2misoutmin<-read.csv("EffectsOutcome1ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2misoutmin<-read.csv("CoverageATEOutcome1ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S2misoutmaj1<-read.csv("EffectsOutcome2ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2misoutmaj1<-read.csv("CoverageATEOutcome2ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S2misoutmaj2<-read.csv("EffectsOutcome3ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2misoutmaj2<-read.csv("CoverageATEOutcome3ModelScenTwo3levels.csv", header=TRUE, stringsAsFactors=FALSE)

dat.S2mistx2<-read.csv("EffectsScenTwo3levelsflucsvyMMTxMis.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2mistx2<-read.csv("CoverageATEScenTwo3levelsflucsvyMMTxMis.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S2missel2<-read.csv("EffectsScenTwo3levelsflucsvyMMSelMis.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2missel2<-read.csv("CoverageATEScenTwo3levelsflucsvyMMSelMis.csv", header=TRUE, stringsAsFactors=FALSE)
dat.S2misout2<-read.csv("EffectsScenTwo3levelsflucsvyMMOutMis.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2misout2<-read.csv("CoverageATEScenTwo3levelsflucsvyMMOutMis.csv", header=TRUE, stringsAsFactors=FALSE)

dat.S2mistx3<-read.csv("EffectsScenTwo3levelsredoMMTxMis.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2mistx3<-read.csv("CoverageATEScenTwo3levelsredoMMTxMis.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2mistx4<-read.csv("CoverageATEScenTwo3levelsTxMissvy.csv", header=TRUE, stringsAsFactors=FALSE)

dat.S2missel3<-read.csv("EffectsScenTwo3levelsredoMMSelMis.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2missel3<-read.csv("CoverageATEScenTwo3levelsredoMMSelMis.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2missel4<-read.csv("CoverageATEScenTwo3levelsMisSelsvy.csv", header=TRUE, stringsAsFactors=FALSE)

dat.S2misout4<-read.csv("EffectsScenTwo3levelsredoMMOutMis.csv", header=TRUE, stringsAsFactors=FALSE)
cover.S2misout4<-read.csv("CoverageATEScenTwo3levelsredoMMOutMis.csv", header=TRUE, stringsAsFactors=FALSE)

tmleS2box<-c(
  dat.S2mistxmin$TxSvySub, dat.S2mistx3$TMLESvy1T1, dat.S2mistx3$TMLESvy2T1, dat.S2mistx2$TMLESvyFlucT1, 
  dat.S2mistxmaj$TxSvySub, dat.S2mistx3$TMLESvy1T2, dat.S2mistx3$TMLESvy2T2, dat.S2mistx2$TMLESvyFlucT2, 
  dat.S2misselmin$TxSvySub, dat.S2missel3$TMLESvy1S1, dat.S2missel3$TMLESvy2S1, dat.S2missel2$TMLESvyFlucS1, 
  dat.S2misselmaj$TxSvySub, dat.S2missel3$TMLESvy1S2, dat.S2missel3$TMLESvy2S2, dat.S2missel2$TMLESvyFlucS2,
  dat.S2misoutmin$TxSvySub, dat.S2misout4$TMLESvy1O1, dat.S2misout4$TMLESvy2O1, dat.S2misout2$TMLESvyFlucO1,
  dat.S2misoutmaj1$TxSvySub, dat.S2misout4$TMLESvy1O2, dat.S2misout4$TMLESvy2O2, dat.S2misout2$TMLESvyFlucO2,
  dat.S2misoutmaj2$TxSvySub, dat.S2misout4$TMLESvy1O3, dat.S2misout4$TMLESvy2O3, dat.S2misout2$TMLESvyFlucO3)
tmleS2cover<-c(
  cover.S2mistxmin$TxSvySub, cover.S2mistx3$TMLESvy1T1, cover.S2mistx3$TMLESvy2T1, cover.S2mistx2$TMLESvyFlucT1, 
  cover.S2mistxmaj$TxSvySub, cover.S2mistx4$TMLESvyT2, cover.S2mistx3$TMLESvy2T2, cover.S2mistx2$TMLESvyFlucT2, 
  cover.S2misselmin$TxSvySub, cover.S2missel3$TMLESvy1S1, cover.S2missel3$TMLESvy2S1, cover.S2missel2$TMLESvyFlucS1, 
  cover.S2misselmaj$TxSvySub, cover.S2missel4$TMLESvyS2, cover.S2missel3$TMLESvy2S2, cover.S2missel2$TMLESvyFlucS2,
  cover.S2misoutmin$TxSvySub, cover.S2misout4$TMLESvy1O1, cover.S2misout4$TMLESvy2O1, cover.S2misout2$TMLESvyFlucO1,
  cover.S2misoutmaj1$TxSvySub, cover.S2misout4$TMLESvy1O2, cover.S2misout4$TMLESvy2O2, cover.S2misout2$TMLESvyFlucO2,
  cover.S2misoutmaj2$TxSvySub, cover.S2misout4$TMLESvy1O3, cover.S2misout4$TMLESvy2O3, cover.S2misout2$TMLESvyFlucO3)
labS2<-c(rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250), rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250),rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250),rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250),rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250),rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250),rep("IPTSW", 250), rep("TMLE1", 250), rep("TMLE2", 250),  rep("TMLE3", 250))
catS2<-c(rep("modtx", 1000), rep("majtx", 1000), rep("modsel", 1000), rep("majsel", 1000), rep("modout", 1000), rep("majout1", 1000),rep("majout2", 1000) )
tmleS2box.dat<-data.frame(tmleS2box, labS2, catS2)
tmleS2cover.dat<-data.frame(tmleS2cover, labS2, catS2)

sum.tmleS2 <- ddply(tmleS2box.dat, .(labS2, catS2), summarise, 
  AbsBias = mean(tmleS2box-ate),
    PcentBias = mean((tmleS2box-ate)/ate),
     MSE=mean((tmleS2box - ate)^2))

tmleS2.box.plot<-ggplot(tmleS2box.dat, aes(x=catS2, y=tmleS2box)) + geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", shape=5, size=4) + facet_grid(~ labS2 ) +
    coord_flip() + geom_hline(yintercept=ate, colour="red")  + scale_x_discrete(limits=c("modtx", "majtx", "modsel", "majsel","modout", "majout1", "majout2")) + ylab("ATE") +xlab("") + ggtitle("ATE, Scenario 2")+ theme( title=element_text(size=9))

sumcover.tmleS2<- ddply(tmleS2cover.dat, .(labS2, catS2), summarise, 
  meancov = mean(tmleS2cover))
tmleS2.cov.plot<-ggplot(data=tmleS2cover.dat, aes(x=catS2, y=tmleS2cover))  +facet_grid(~ labS2 ) +
    geom_bar( stat="identity", fill="lightblue") +
scale_x_discrete(limits=c("modtx", "majtx", "modsel", "majsel","modout", "majout1", "majout2")) + 
   ggtitle("Coverage, Scenario 2") + xlab("") + ylab("") + theme(axis.text.x=element_text(size=7), title=element_text(size=9))

grid.arrange(tmleS2.box.plot, tmleS2.cov.plot, nrow=1)

@

<<echo=FALSE, fig=TRUE>>=
tmletrue.cov.plot
@

<<echo=FALSE, results=tex>>=
xtable(sum.tmleS2, type="latex", digits=3)
@

<<echo=FALSE, fig=TRUE>>=
tmlemis.box.plot
@

<<echo=FALSE, fig=TRUE>>=
tmlemis.cov.plot
@

<<echo=FALSE, results=tex>>=
xtable(sum.tmleS1, type="latex", digits=3)
@