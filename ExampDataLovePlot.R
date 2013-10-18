setwd("/Users/kararudolph/Documents/PhD/NIMH/NCSA/cortisol/simulationproject")
load("exampsvysubdatOct2013.Rdata")
# Make plots of standardized mean differences
library(plyr)
svysub$pblack<-ifelse(svysub$racecat==1,1,0)
svysub$phispanic<-ifelse(svysub$racecat==0,1,0)
svysub$pwhite<-ifelse(svysub$racecat==3,1,0)
svysub$pother<-ifelse(svysub$racecat==2,1,0)

svysub$pLTHS<-ifelse(svysub$meducat==0,1,0)
svysub$pHSG<-ifelse(svysub$meducat==1,1,0)
svysub$pSC<-ifelse(svysub$meducat==2,1,0)
svysub$pCG<-ifelse(svysub$meducat==3,1,0)

svysub$northeast<-ifelse(svysub$region==1,1,0)
svysub$south<-ifelse(svysub$region==3,1,0)
svysub$west<-ifelse(svysub$region==4,1,0)
svysub$midwest<-ifelse(svysub$region==2,1,0)

svysub$spring<-ifelse(svysub$season==0,1,0)
svysub$summer<-ifelse(svysub$season==1,1,0)
svysub$fall<-ifelse(svysub$season==2,1,0)
svysub$winter<-ifelse(svysub$season==3,1,0)

svysub$bmi<-703*(svysub$lbs / (svysub$feet*12 + svysub$inch)^2)

svysub$immigrant<-ifelse(svysub$imgen==1,1,0)

datplot<-svysub[,c(10:14, 17,19, 21, 40, 51,53, 64, 67:69, 76,77,80,91, 93,94,107,110:126)]
a<-c()
for(i in 1:length(names(datplot))){
	a[i]<-class(datplot[,i])
}

b<-datplot[,c(1,2,11,12, 15)]
a<-datplot[,c(3:10,13:14, 16:39)]
c<-apply(b, c(1,2), function(x) as.numeric(as.character(x)))
d<-cbind(a,c)

txmeans<-daply(d, .(tertscore), colwise(mean))
txvars<-daply(d, .(tertscore), colwise(var))

stdmeandif<-(unlist(txmeans[2,]) - unlist(txmeans[1,]))/sqrt(unlist(txvars[2,]))
nam2<-c("South","Black", "Hispanic", "ESL", "MomHSGrad", "suburb", "bmi", "numTruama", "summer", "smallgestage", "weekend", "immigrant","wkdySleep", "wkndBedtime",  "MomWork","OtherRace",  "spring","fall", "gender", "numRx",  "West", "age", "wkdyBedtime", "wkndSleep",  "winter", "MomSmCollege",  "livewthMom", "citizen", "timeSamp", "curremp", "urban", "MomAge",  "Northeast", "livewthDad",    "income", "MomColGrad",   "Midwest", "White"        )
orderedpear<-stdmeandif[order(stdmeandif, decreasing=TRUE)]
orderedpear1<-(as.numeric(orderedpear)*100)

load("exampsvydatOct2013.Rdata")

svysamp2<-b
svysamp2$pblack<-ifelse(svysamp2$racecat==1,1,0)
svysamp2$phispanic<-ifelse(svysamp2$racecat==0,1,0)
svysamp2$pwhite<-ifelse(svysamp2$racecat==3,1,0)
svysamp2$pother<-ifelse(svysamp2$racecat==2,1,0)

svysamp2$pLTHS<-ifelse(svysamp2$meducat==0,1,0)
svysamp2$pHSG<-ifelse(svysamp2$meducat==1,1,0)
svysamp2$pSC<-ifelse(svysamp2$meducat==2,1,0)
svysamp2$pCG<-ifelse(svysamp2$meducat==3,1,0)

svysamp2$northeast<-ifelse(svysamp2$region==1,1,0)
svysamp2$south<-ifelse(svysamp2$region==3,1,0)
svysamp2$west<-ifelse(svysamp2$region==4,1,0)
svysamp2$midwest<-ifelse(svysamp2$region==2,1,0)

svysamp2$spring<-ifelse(svysamp2$season==0,1,0)
svysamp2$summer<-ifelse(svysamp2$season==1,1,0)
svysamp2$fall<-ifelse(svysamp2$season==2,1,0)
svysamp2$winter<-ifelse(svysamp2$season==3,1,0)

svysamp2$bmi<-703*(svysamp2$lbs / (svysamp2$feet*12 + svysamp2$inch)^2)

svysamp2$immigrant<-ifelse(svysamp2$imgen==1,1,0)

dat2plot<-svysamp2[,c(3,4,6:8,11,12,16,30,32,34,40,46,50:57,62,64,68,113:130)]

c("South","Black", "Hispanic", "ESL", "MomHSGrad", "suburb", "bmi", "numTruama", "summer", "smallgestage", "weekend", "immigrant","wkdySleep", "wkndBedtime",  "MomWork","OtherRace",  "spring","fall", "gender", "numRx",  "West", "age", "wkdyBedtime", "wkndSleep",  "winter", "MomSmCollege",  "livewthMom", "citizen", "timeSamp", "curremp", "urban", "MomAge",  "Northeast", "livewthDad",    "income", "MomColGrad",   "Midwest", "White"        )

a<-c()
for(i in 1:length(names(dat2plot))){
	a[i]<-class(dat2plot[,i])
}

b1<-dat2plot[,c(1,2,9,20:22)]
a1<-dat2plot[,c(3:8,10:19,23:42)]

c1<-apply(b1, c(1,2), function(x) as.numeric(as.character(x)))
d1<-cbind(a1,c1)
d2<-d1[,-c(8,11,12)]

selmeans<-daply(d2, .(insample), colwise(mean))
selvars<-daply(d2, .(insample), colwise(var))

stdmeandifsel<-(unlist(selmeans[2,]) - unlist(selmeans[1,]))/sqrt(unlist(selvars[2,]))

peach1<-(as.numeric(stdmeandifsel)*100)

pdf("loveplotexampdata.pdf")
#postscript("Fig1Color.eps", paper="special", height=8, width=8, horizontal=FALSE, colormodel="rgb")
#postscript("Fig1BW.eps", paper="special", height=8, width=8, horizontal=FALSE, colormodel="gray")
dotchart(orderedpear1, pch="", labels=nam2, cex=0.75) 
mtext("Standardized Mean Difference (%)", side=1, line=2)
points(orderedpear1, seq(1:length(orderedpear1)), pch=19, col="blue", cex=1.2)
points(peach1, seq(1:length(peach1)), pch=21, col="red", cex=1.2) 
abline(v=0, lty=1) 
abline(v=10, lty=2, lwd=1, col="grey") 
abline(v=-10, lty=2, lwd=1, col="grey")
abline(v=20, lty=2, lwd=1, col="black") 
abline(v=-20, lty=2, lwd=1, col="black")
dev.off()

save(b, file="exampsvydat.Rdata")
save(svysub, file="exampsvysubdat.Rdata")

