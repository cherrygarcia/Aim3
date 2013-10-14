## get data for survey subsample
load("datratelessexclcat5.Rdata")

datsecexclcat5<-list(rep(NA,10))
for(i in 1:10){
	datsecexclcat5[[i]]<-datratelessexclcat5[[i]][(datratelessexclcat5[[i]]$endtime - datratelessexclcat5[[i]]$begtime)/(60*60)<4, ]
	datsecexclcat5[[i]]$racecat1[datsecexclcat5[[i]]$racecat==1]<-0
	datsecexclcat5[[i]]$racecat1[datsecexclcat5[[i]]$racecat==2]<-1
	datsecexclcat5[[i]]$racecat1[datsecexclcat5[[i]]$racecat==3]<-2
	datsecexclcat5[[i]]$racecat1[datsecexclcat5[[i]]$racecat==4]<-3
	datsecexclcat5[[i]]$racecat<-datsecexclcat5[[i]]$racecat1
}

svysub<-datsecexclcat5[[1]]
svysub$insample<-1
save(svysub, file="exampsvysubdatOct2013.Rdata")

## make dataset for survey sample
library(mice)
#read in the data
setwd("/Users/kararudolph/Documents/PhD/NIMH/NCSA/cortisol")
#get date and time of interview for everyone
five<-read.csv("nshs_tot_datetime.csv", header=TRUE)
five$SampleID<-five$SID
five$t1<-as.character(five$StartDateTime)
five$t1t<-as.POSIXlt(strptime(five$t1, "%m/%d/%y %I:%M %p"))
#identify data entry errors
#cortisol sample too early
w<-five$t1t
w$hour2<-ifelse(w$hour<6, w$hour+12, w$hour)
w$min2<-ifelse(w$min<10, as.character(paste(0, w$min, sep="")), w$min)
w$time<-as.character(paste (w$hour2, w$min2, sep = ":"))
five$time<-w$time
# get the interview time in the right format and have the date match the date of those in the survey subsample
five$time3<-paste("2002-01-01", five$time)
five$t1t2<-strptime(five$time3,"%Y-%m-%d %H:%M")
five$first<-as.POSIXlt(five$t1t2)
five$begtime<-as.numeric(five$first)

keep2<-c("SampleID", "first", "t1t", "begtime")
positx.vars<-five[!is.na(five$begtime), keep2]

load("cortimpwholesamp.Rdata")
cortimpwholesamp1<-complete(cortimpwholesamp, action=1, include=FALSE)
cortimpwholesamp2<-complete(cortimpwholesamp, action=2, include=FALSE)
cortimpwholesamp3<-complete(cortimpwholesamp, action=3, include=FALSE)
cortimpwholesamp4<-complete(cortimpwholesamp, action=4, include=FALSE)
cortimpwholesamp5<-complete(cortimpwholesamp, action=5, include=FALSE)
cortimpwholesamp6<-complete(cortimpwholesamp, action=6, include=FALSE)
cortimpwholesamp7<-complete(cortimpwholesamp, action=7, include=FALSE)
cortimpwholesamp8<-complete(cortimpwholesamp, action=8, include=FALSE)
cortimpwholesamp9<-complete(cortimpwholesamp, action=9, include=FALSE)
cortimpwholesamp10<-complete(cortimpwholesamp, action=10, include=FALSE)

cortimpwholesampten<-list(cortimpwholesamp1, cortimpwholesamp2, cortimpwholesamp3, cortimpwholesamp4, cortimpwholesamp5, cortimpwholesamp6, cortimpwholesamp7, cortimpwholesamp8, cortimpwholesamp9, cortimpwholesamp10)

tmp.excl<-list(rep(NA, 10))
tmp2.excl<-list(rep(NA, 10))

## add back in date/time variables
for(i in 1:10){
tmp.excl[[i]]<-merge(cortimpwholesampten[[i]], positx.vars, by="SampleID", all.x=FALSE, all.y=FALSE) }

## exclude those who had measures taken during CAR
for(i in 1:10){
	tmp.excl[[i]]$schyr<-ifelse(with(tmp.excl[[i]], t1t$mo==6 | t1t$mo==7 | t1t$mo==8), 0, 1 )
	tmp.excl[[i]]$late<-ifelse(tmp.excl[[i]]$schyr==0 | tmp.excl[[i]]$weekend==1, 1, 0)
	tmp.excl[[i]]$wm[tmp.excl[[i]]$late==0 & with(tmp.excl[[i]], first$hour<5 | first$hour>9 )]<-1
	tmp.excl[[i]]$wm[tmp.excl[[i]]$late==1 & with(tmp.excl[[i]], first$hour<7 | first$hour>11 )]<-1
}

for(i in 1:10){
tmp2.excl[[i]]<-tmp.excl[[i]][!is.na(tmp.excl[[i]]$wm) & tmp.excl[[i]]$wm==1,]}

## we could also further limit to those kids who would not be excluded by cortisol studies
imp.exclcat5<-list(rep(NA, 10))
for(i in 1:10){
imp.exclcat5[[i]]<-tmp2.excl[[i]][tmp2.excl[[i]]$asthmatrt==0 & tmp2.excl[[i]]$pregnant==0 & tmp2.excl[[i]]$diabetes==0 & tmp2.excl[[i]]$psychmeds==0 & tmp2.excl[[i]]$oc==0 & tmp2.excl[[i]]$smoke==0 & tmp2.excl[[i]]$currentdruguse==0,]
} 
for(i in 1:10){
imp.exclcat5[[i]]$insample<-ifelse(imp.exclcat5[[i]]$SampleID %in% svysub$SampleID, 1, 0)
}

setwd("/Users/kararudolph/Documents/PhD/NIMH/NCSA/cortisol/simulationproject")

svynosub<-imp.exclcat5[[1]][imp.exclcat5[[1]]$insample==0,]

a<-merge(svynosub, svysub,all.x=TRUE, all.y=TRUE)
b<-a
save(b, file="exampsvydatOct2013.Rdata")