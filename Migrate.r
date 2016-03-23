setwd('~/Documents/GitHub/STModel-CompAnalysis/')

load('./out/4-LandMig/bandprop.rdata')
bandprop <- res
load('./out/4-LandMig/timetrics.rdata')
timetrics <- res
# library
library(reshape2)

# Create fig with the proportion along the gradient
bandprop[,c('B','T','R','M')] <-bandprop[,c('B','T','R','M')]/rowSums(bandprop[,c('B','T','R','M')])
bandprop$TM <- rowSums(bandprop[,c('T','M')])

agg_bandprop <- aggregate(TM~y+model+year,data=bandprop,mean)
agg_bandprop$y <- agg_bandprop$y/1000
agg_bandprop <- dcast(y~model+year,data=agg_bandprop,value.var="TM")

par(mar=c(4,4,1,1),oma=c(1,1,1,1))
par(mfrow=c(2,3))
plot(agg_bandprop$y,agg_bandprop$local_2015,col="#FEAC18",lwd=2,type="l",ylim=c(0,1))
plot(agg_bandprop$y,agg_bandprop$local_2045,col="#FEAC18",lwd=2,type="l",ylim=c(0,1))
plot(agg_bandprop$y,agg_bandprop$local_2095,col="#FEAC18",lwd=2,type="l",ylim=c(0,1))
plot(agg_bandprop$y,agg_bandprop$solved_2015,col="dodgerblue4",lwd=2,type="l",ylim=c(0,1))
plot(agg_bandprop$y,agg_bandprop$solved_2045,col="dodgerblue4",lwd=2,type="l",ylim=c(0,1))
plot(agg_bandprop$y,agg_bandprop$solved_2095,col="green4",lwd=2,type="l",ylim=c(0,1))

par(mfrow=c(2,3))
for(i in 2:ncol(bandprop)) hist(bandprop[,i],main=names(bandprop)[i])
# polygon(x=c(local_2015$lat,rev(local_2015$lat)),y=c(local_2015$avg_prop+local_2015$sd_prop,rev(local_2015$avg_prop-local_2015$sd_prop)),border=NA,col=rgb(254,172,1,alpha=40,max=255))
# polygon(x=c(local_2095$lat,rev(local_2095$lat)),y=c(local_2095$avg_prop+local_2095$sd_prop,rev(local_2095$avg_prop-local_2095$sd_prop)),border=NA,col=rgb(255,0,0,alpha=40,max=255))
# polygon(x=c(solved2095$lat,rev(solved2095$lat)),y=c(solved2095$avg_prop+solved2095$sd_prop,rev(solved2095$avg_prop-solved2095$sd_prop)),border=NA,col=rgb(16,78,139,alpha=40,max=255))
# polygon(x=c(equi_2095$lat,rev(equi_2095$lat)),y=c(equi_2095$avg_prop+equi_2095$sd_prop,rev(equi_2095$avg_prop-equi_2095$sd_prop)),border=NA,col=rgb(0,139,0,alpha=40,max=255))
# legend(-10,0.2, c("STM - 2015", "STM - 2095", "STM - 2095 (+10.000)","STM at equilibrium"), col = c("#FEAC18", "red","green4","dodgerblue4"), lty = rep(1,4),lwd=rep(2,4),merge = TRUE,cex=0.9,bty = "n")
dev.off()

### Compute metrics 

# Rate of T+M
library(reshape2)

sub_metric <- timetrics[,c("state" ,"model","pars","landCI", "gcm","rep","year", "Q95" )]
sub_metric <- subset(sub_metric,model!="t10000" & (year=='2000' | year=='2095'))
sub_metric <- aggregate(Q95~model+pars+landCI+gcm+rep+year,data=sub_metric,max)
sub_metric <- dcast(data=sub_metric,model+pars+landCI+gcm+rep~year,value.var="Q95")
names(sub_metric)[(ncol(sub_metric)-1):ncol(sub_metric)] <- c('t0','t1')
sub_metric$delta <- sub_metric$t1 -sub_metric$t0

avg_metric <- aggregate(delta~model,data=sub_metric,mean)
sd_metric <- aggregate(delta~model,data=sub_metric,sd)
metric <- merge(avg_metric,sd_metric,by=c('model'))
names(metric)[2:3] <- c("avg","sd")
metric[,2:3] <- metric[,2:3]/1000


























