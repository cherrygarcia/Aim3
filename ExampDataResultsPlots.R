exampmatrix<-read.csv("ExampleDataEffectsNov1.csv", header=TRUE)

a<-data.frame(exampmatrix[-c(1:4),])
colnames(a)<-c("method", "mean", "LCI", "UCI")
a$estim<-factor(c("IPW", "DRWLS", "TMLE", "IPW", "DRWLS", "TMLE","IPW", "DRWLS", "TMLE" ,"IPW", "DRWLS", "TMLE"), levels=c("IPW", "DRWLS", "TMLE"))
a$deg<-factor(c(rep("A", 3), rep("B", 3), rep("C", 3), rep("D", 3)), levels=c("A", "B", "C", "D"))
b<-a
b[,2]<-as.numeric(as.character(b[,2]))
b[,3]<-as.numeric(as.character(b[,3]))
b[,4]<-as.numeric(as.character(b[,4]))
pdf("ExampEst95FigSens.pdf")
ggplot(b, aes(x=estim, y=mean)) + facet_wrap(~deg) + geom_errorbar(aes(ymin=LCI, ymax=UCI)) + geom_point() + theme_bw() + ylab("Cortisol Rate Difference (ng/ml/hr)") + xlab("")+ geom_hline(aes(yintercept=0))
dev.off()

d<-data.frame(exampmatrix[-c(8:16),])
d[,2]<-as.numeric(as.character(d[,2]))
d[,3]<-as.numeric(as.character(d[,3]))
d[,4]<-as.numeric(as.character(d[,4]))
colnames(d)<-c("method", "mean", "LCI", "UCI")
d$method<-factor(c("Naive", "IPTW", "IPTSvyW", "IPTSubW", "IPW", "DRWLS", "TMLE"), levels=c("Naive", "IPTW", "IPTSvyW", "IPTSubW", "IPW", "DRWLS", "TMLE"))
plotcis<-ggplot(d, aes(x=method, y=mean)) + 
    geom_errorbar(aes(ymin=LCI, ymax=UCI), width=.1) +
    geom_point() + theme_bw() + ylab("Cortisol Rate Difference (ng/ml/hr)") + xlab("")  + geom_hline(aes(yintercept=0)) 

pdf("ExampEst95Fig.pdf")
plotcis
dev.off()