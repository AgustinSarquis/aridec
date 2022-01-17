library(aridec)
db=loadEntries()

# Fig 3 a: pubblications per year
citationKey=lapply(db, FUN=function(x){x$citationKey})
year=as.numeric(gsub("[^[:digit:]]", "", citationKey))
library(plyr)
freq.year=count(year)
library(ggplot2)
fig3a=ggplot(data=freq.year, aes(x, freq)) + geom_bar(stat= "identity", fill="blue") +
  xlab("Year of Pubblication") + ylab("Number of pubblications") + ggtitle("(a)") +theme_bw() +
  scale_x_continuous(breaks=c(1975,1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=30))

# Fig 3 b: length of study period
duration=as.data.frame(unlist(lapply(db, FUN=function(x){x$experimentInfo$duration})))
colnames(duration) = c( "days")
durationmedian=median(duration$days)
durationmean=mean(duration$days)
fig3b=ggplot(data=duration, aes(days)) + geom_histogram(binwidth= 200,fill="blue") +
 xlab("Duration of study (days)") + ylab("Number of pubblications") + ggtitle("(b)") +theme_bw() +
 scale_x_continuous(breaks=c(300,600, 900,1200,1800, 3600)) +
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=30)) +
 geom_vline(xintercept=durationmean, linetype= "dashed") + geom_vline(xintercept = durationmedian, linetype="dotted")

# Fig 3 c: number of sampling dates
NAomitseries=function(x){lapply(x, na.omit)}
harvestentry= function(x) {lapply(NAomitseries(x$timeSeries), length)}
harvestdb=function(x){lapply(x, harvestentry)}
harvests=harvestdb(db)
harvestsfinal=unlist(lapply(harvests, FUN=function(x){x[-1]}))
harvestmean=mean(harvestsfinal)
harvestmedian=median(harvestsfinal)
harvestsfinal=as.data.frame(harvestsfinal)
fig3c=ggplot(data=harvestsfinal, aes(harvestsfinal)) + geom_histogram(binwidth= 1,fill="blue") +
  xlab("Number of harvests") + ylab("Number of pubblications") + ggtitle("(c)") +theme_bw() +
  scale_x_continuous(breaks=c(5,10,15,20,23)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=30)) +
  geom_vline(xintercept=harvestmean, linetype= "dashed") + geom_vline(xintercept = harvestmedian, linetype="dotted")

# Fig 3 d: sampling frequency 
# using objects from the duration and sampling dates figures
meanharvest=as.data.frame(unlist(lapply(harvests, FUN=function(x){x$Time-1})))
colnames(meanharvest)="harvests"
harvestfreq=as.data.frame(cbind(duration$days, meanharvest$harvests))
colnames(harvestfreq)=c("time", "harvests")
harvestfreq=mutate(harvestfreq, time=time/30) 
frequency=harvestfreq$harvests/harvestfreq$time
harvestfreq=cbind(harvestfreq, frequency)
fig3d=ggplot(data=harvestfreq, aes(frequency)) + geom_histogram(binwidth= 0.5,fill="blue") +
  xlab("Sampling frequency (harvest . months^-1)") + ylab("Number of studies") + ggtitle("(d)") +theme_bw() +
  scale_x_continuous(breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8 ,9, 10, 11)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=30)) +
  geom_vline(xintercept=mean(harvestfreq$frequency), linetype= "dashed") + geom_vline(xintercept = median(harvestfreq$frequency), linetype="dotted")
fig3d

# Fig 3 e: elevation
elevation=as.data.frame(unlist(lapply(db, FUN=function(x){x$siteInfo$elevation})))
colnames(elevation) = c( "elevation")
fig3e=ggplot(data=elevation, aes(elevation)) + geom_histogram(binwidth= 200,fill="blue") +
  xlab("Study site elevation (m a.s.l.)") + ylab("Number of sites") + ggtitle("(e)") +theme_bw() +
  scale_x_continuous(breaks=c(-400, 1000, 2000, 3000, 4000)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=30)) +
  geom_vline(xintercept=mean(elevation$elevation), linetype= "dashed") + geom_vline(xintercept = median(elevation$elevation), linetype="dotted")

# Fig 3 f: MAT
MAT=as.data.frame(unlist(lapply(db, FUN=function(x){x$siteInfo$MAT})))
colnames(MAT) = c( "MAT")
fig3f=ggplot(data=MAT, aes(MAT)) + geom_histogram(binwidth= 5,fill="blue") +
  xlab("Mean annual temperature (°C)") + ylab("Number of sites") + ggtitle("(f)") +theme_bw() +
  scale_x_continuous(breaks=c(0,5,10,15,20,25,30)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=30)) +
  geom_vline(xintercept=mean(MAT$MAT), linetype= "dashed") + geom_vline(xintercept = median(MAT$MAT), linetype="dotted")

# Fig 3 g: MAP
MAP=as.data.frame(unlist(lapply(db, FUN=function(x){x$siteInfo$MAP})))
colnames(MAP) = c( "MAP")
fig3g=ggplot(data=MAP, aes(MAP)) + geom_histogram(binwidth= 150,fill="blue") +
  xlab("Mean annual precipitation (mm)") + ylab("Number of sites") + ggtitle("(g)") +theme_bw() +
  scale_x_continuous(breaks=c(250, 500, 750, 1000, 1250, 1500, 1750)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=30)) +
  geom_vline(xintercept=mean(MAP$MAP), linetype= "dashed") + geom_vline(xintercept = median(MAP$MAP), linetype="dotted")

# Fig 3 h: land cover
cover=as.data.frame(unlist(lapply(db, FUN=function(x){x$siteInfo$landCover})))
colnames(cover) = c( "cover")
freq.cover=count(cover)
freq.cover[3, 2]=56 # desert + sandland
freq.cover[4,2]=38 # farmland + abbandoned farmland
freq.cover[2,1]="coastal"
freq.cover[2,2]=4 # coastal + mangrove + marsh
freq.cover[5,2]=41 # forest + greenbelt
freq.cover[4,1]="agroecosystem"
freq.cover=freq.cover[-c(1, 7, 8, 9, 10),]
freq.cover=cbind(freq.cover, rel.freq=freq.cover$freq/238)
fig3h=ggplot(data=freq.cover, aes(reorder(cover, -freq), freq)) + geom_bar(stat= "identity", fill="blue") +
  xlab("Land cover type") + ylab("Number of sites") + ggtitle("(h)") +theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle=30))
fig3h

library(gridExtra)
windows()
grid.arrange(fig3a,fig3b, fig3c, fig3d, fig3e, fig3f, fig3g, fig3h,ncol=2)


