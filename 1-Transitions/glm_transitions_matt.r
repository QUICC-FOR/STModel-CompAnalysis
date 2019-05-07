# Load library , data
load("transitions_r5.rdata")
library("nnet")
library("data.table")
library("knitr")

#0. Filter data
transitionData <- as.data.table(subset(transitionData,subset=annual_mean_temp<=10))


## todo first:
## use offset() time interval and compare to subsetting by 5-year intervals
## needs to be in the formula: y ~ temp + pp + offset(time)
## offsets don't really add much, so instead use 5-year intervals


#1. Scale variables
#2. Add time interval
# scale_tp<- scale(transitionData$annual_mean_temp)
# transitionData$annual_mean_temp <- as.numeric(scale_tp)
# scale_pp <- scale(transitionData$tot_annual_pp)
# transitionData$annual_pp <- as.numeric(scale_pp)
# transitionData$time <- transitionData$year2 - transitionData$year1
vars <- cbind(poly(transitionData$annual_mean_temp, 3), poly(transitionData$tot_annual_pp, 3))
colnames(vars) <- c(paste0("temp", 1:3), paste0("pp", 1:3))
vars <- apply(vars, 2, function(x) x / sd(x))

modelVars <- cbind(transitionData[,.(plot, year1, year2, 
	interval = year2-year1, state1, state2)], vars)

#3. Remove imposible transitions & filter on time intervals
# transitionData$trobs <- paste0(transitionData$state1,transitionData$state2)
# transitionData <- subset(transitionData,trobs!="BT" & trobs!="TB")
modelVars <- modelVars[!((state1 == 'B' & state2=='T') | (state1 == 'T' & state2 == 'B'))]
modelVars <- modelVars[interval <= 10]

# now create subset data tables with states recoded to get nice parameters
modelVarsT <- modelVars[state1 == 'T' & interval == 5]
modelVarsT[,state2 := factor(state2, levels=c('T', 'M', 'R'))]

modelVarsB <- modelVars[state1 == 'B' & interval == 5]
modelVarsB[,state2 := factor(state2, levels=c('B', 'M', 'R'))]

modelVarsM <- modelVars[state1 == 'M' & interval == 5]
modelVarsM[,state2 := factor(state2, levels=c('M', 'B', 'T', 'R'))]

modelVarsR <- modelVars[state1 == 'R' & interval == 5]
modelVarsR[,state2 := factor(state2, levels=c('R', 'B', 'T', 'M'))]


T_multi_5 <- multinom(state2 ~ temp1 + temp2 + temp3 + pp1 + pp2 + pp3, data = modelVarsT)
B_multi_5 <- multinom(state2 ~ temp1 + temp2 + temp3 + pp1 + pp2 + pp3, data = modelVarsB)
M_multi_5 <- multinom(state2 ~ temp1 + temp2 + temp3 + pp1 + pp2 + pp3, data = modelVarsM)
R_multi_5 <- multinom(state2 ~ temp1 + temp2 + temp3 + pp1 + pp2 + pp3, data = modelVarsR)

T_multi_2o <- multinom(state2 ~ temp1 + temp2 + pp1 + pp2, data = modelVarsT)
B_multi_2o <- multinom(state2 ~ temp1 + temp2 + pp1 + pp2, data = modelVarsB)
M_multi_2o <- multinom(state2 ~ temp1 + temp2 + pp1 + pp2, data = modelVarsM)
R_multi_2o <- multinom(state2 ~ temp1 + temp2 + pp1 + pp2, data = modelVarsR)

# compute delta AIC for dropping the third order term; drop if > 10
AIC(T_multi_5) - AIC(T_multi_2o) 
AIC(B_multi_5) - AIC(B_multi_2o) 
AIC(M_multi_5) - AIC(M_multi_2o) 
AIC(R_multi_5) - AIC(R_multi_2o) 

## Summarize information in table and export in latex
## nested Rbinds to make an ugly table less ugly
##
outTab <- signif(rbind(
	rbind(coefficients(T_multi_5), summary(T_multi_5)$standard.errors)[c(1,3,2,4),],
	rbind(coefficients(B_multi_5), summary(B_multi_5)$standard.errors)[c(1,3,2,4),],
	rbind(coefficients(M_multi_5), summary(M_multi_5)$standard.errors)[c(1,4,2,5,3,6),],
	rbind(coefficients(R_multi_5), summary(R_multi_5)$standard.errors)[c(1,4,2,5,3,6),]), digits=3)
rownames(outTab) <- paste0(c(rep('T', 4), rep('B', 4), rep('M', 6), rep('R', 6)), rownames(outTab), c('.mean', '.se'))
kable(outTab,format="latex")

