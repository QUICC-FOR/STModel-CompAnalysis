# Validation via TSS
rm(list=ls())
# Setwd
setwd("~/Documents/Github/STModel-CompAnalysis")
load('./extra/RandomForest_complete.robj')
load('./extra/transitions_r1.rdata')
load('./extra/scale_info.robj')

source('./fcts/fcts_resAnalytic.r')

## Function
# evaluation statistics (revised from Boulangeat et al. 2012 to handle large numbers)
HK <- function(Pred, Obs)
{
	Misc = table(Pred, Obs)

    if (nrow(Misc)!=ncol(Misc)) stop("wrong misclassification table")
    Misc <- unclass(Misc)
    k  <- ncol(Misc)
    Nobs <- apply(Misc, 2, sum)
    Npred <- apply(Misc, 1, sum)
    N <- sum(Nobs)


   HK <- (sum(diag(Misc))/N - sum(as.numeric(Nobs)*as.numeric(Npred))/N/N ) / ( 1 - sum(as.numeric(Nobs)*as.numeric(Nobs))/N/N )

    return(HK)√ 
}

#### TSS EQUILIBRIUM
# Libraries
library(randomForest)

# Filter plots: 1 measure and <10 degree
data_plt <- subset(stateData, annual_mean_temp <=10)
id_uq <- table(stateData$plot_id)
id_uq <- names(id_uq[which(id_uq==1)])
uq_plt <- subset(data_plt, plot_id %in% id_uq)

# PrepData for RF
rfdata <- data_plt[,c('annual_mean_temp','tot_annual_pp','mean_diurnal_range','ph_2cm','slp','lat','lon')]
pred <- predict(SDM2,new=rfdata,"prob")

# Data for eq classification
obs_eq <- cbind(uq_plt[,c('lat','lon','annual_mean_temp','tot_annual_pp')])
obs_eq$annual_mean_temp <- (obs_eq$annual_mean_temp - vars.means['annual_mean_temp'])/vars.sd['annual_mean_temp']
obs_eq$tot_annual_pp <- (obs_eq$tot_annual_pp - vars.means['tot_annual_pp'])/vars.sd['tot_annual_pp']
names(obs_eq) <- c('lat','lon','tp','pp')
pars <- read.table("./pars/GenSA_rf_all_2_5y_rep1.txt",row.names=1) # pars

# solve and merge
pred_eq <- solve_stm(obs_eq,pars)
pred_eq$R <- 1 - (pred_eq$T+pred_eq$M+pred_eq$B)
eq <- merge(uq_plt[,c('lon','lat','plot_id','year_measured','state')],pred_eq,by=c('lat','lon'),all=TRUE)
eq$state <- as.factor(eq$state)
eq$pred <- names(eq)[sapply(1:nrow(eq),function(x) which.max(eq[x,c('B','T','M')]))+5]
names(eq)[5] <- 'obs'
eq$pred <- factor(eq$pred,levels=levels(eq$obs))

# Compute HK for eq distribution
HK_eq <- HK(eq$pred,eq$obs)

#### TSS DYNAMIC
data_trans <- subset(transitionData, annual_mean_temp <=10)

# rf
rfdata <- data_trans[,c('annual_mean_temp','tot_annual_pp','mean_diurnal_range','ph_2cm','slp','lat','lon')]
pred <- predict(SDM2,new=rfdata,"prob")

# mod
data_trans <- cbind(data_trans[,c('state1','state2','year1','year2','annual_mean_temp','tot_annual_pp')],pred)
data_trans$itime <- data_trans$year2 - data_trans$year1
names(data_trans) <- c('st0','st1','year1','year2','ENV1','ENV2','EB','EM','ER','ET','itime')




































