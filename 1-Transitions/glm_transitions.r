## Evaluate if transitions are responsive to the climate using multinomiale regression
## By Steve Vissault
## Jeudi 10 février

#Clean wd
rm(list=ls())

# Load library , data
load("~/Documents/Git/stm/STModel-Data-4S/out_files/transitions_r1.rdata")
library("nnet")
library("reshape2")

#0. Filter data
transitionData <- subset(transitionData,subset=annual_mean_temp<=10)

#1. Scale variable
scale_tp<- scale(transitionData$annual_mean_temp)
transitionData$annual_mean_temp <- as.numeric(scale_tp)
scale_pp <-  scale(transitionData$tot_annual_pp)
transitionData$annual_pp <- as.numeric(scale_pp)


#2. Add time interval
transitionData$time <- transitionData$year2 - transitionData$year1

#3. Remove imposible transitions
transitionData$trobs <- paste0(transitionData$state1,transitionData$state2)
transitionData <- subset(transitionData,trobs!="BT" & trobs!="TB")

#4. Evaluate contribution of the coef: build all formula
coefs <-c("annual_mean_temp" , "I(annual_mean_temp^2)" ,"I(annual_mean_temp^3)"  , "annual_pp" , "I(annual_pp^2)" ,"I(annual_pp^3)", "annual_mean_temp:annual_pp" , "time")

contrib_form <- as.character()
tested_coefs <- as.character()

for(i in 1:length(coefs)){
  temp_coefs <- coefs[-i]
  tested_coefs  <- append(tested_coefs ,coefs[i])
  form <- paste0("state2 ~ ",paste(temp_coefs,collapse="+"))
  contrib_form <- append(contrib_form,form)
}
# Add full model
full_formula <- paste0("state2 ~ ",paste(coefs,collapse="+"))
null_formula <- "state2 ~ 1"

# Create DF of formula
formulas <- data.frame(tested=tested_coefs,formula=contrib_form,stringsAsFactors=FALSE)

#5. Full Models
T_multi_full <- multinom(full_formula,data=transitionData,subset= state1=='T')
M_multi_full <- multinom(full_formula,data=transitionData,subset= state1=='M')
R_multi_full <- multinom(full_formula,data=transitionData,subset= state1=='R')
B_multi_full <- multinom(full_formula,data=transitionData,subset= state1=='B')

#6. Null Models
T_multi_null <- multinom(null_formula,data=transitionData,subset= state1=='T')
M_multi_null <- multinom(null_formula,data=transitionData,subset= state1=='M')
R_multi_null <- multinom(null_formula,data=transitionData,subset= state1=='R')
B_multi_null <- multinom(null_formula,data=transitionData,subset= state1=='B')

#7. Compute McFadden pseudo-R2
T_Rs <- 1 - (T_multi_full$deviance / T_multi_null$deviance)
M_Rs <- 1 - (M_multi_full$deviance / M_multi_null$deviance)
R_Rs <- 1 - (R_multi_full$deviance / R_multi_null$deviance)
B_Rs <- 1 - (B_multi_full$deviance / B_multi_null$deviance)

#8. Compute contributions to the full model
states <- c("T","M","R","B")
ls_contrib <-list()


for (state in 1:length(states)){
  contrib <- as.numeric()
  for(i in 1:nrow(formulas)){
    multi <-multinom(formulas$formula[i],data=transitionData,subset= state1==states[state])
    contrib <- append(contrib,(multi$AIC-get(paste0(states[state],'_multi_full'))$AIC))
  }
  ls_contrib[[state]] <- contrib
}

#9. Summarize information in table and export in latex
contrib <- do.call("rbind",ls_contrib)
colnames(contrib) <- formulas$tested
rownames(contrib) <- states
contrib <- data.frame(contrib, AIC=c(T_multi_full$AIC, M_multi_full$AIC, R_multi_full$AIC, B_multi_full$AIC), R2=c(T_Rs,M_Rs,R_Rs,B_Rs))

library("knitr")
kable(contrib,format="latex")

#10. Evaluate probability along gradient
library(RColorBrewer)
library(colorRamps)
pal <-colorRampPalette(rev(brewer.pal(11,"RdYlBu")))

# get range for the specific transition
#rg_tp <- range(transitionData[which(transitionData$state1==state),]$annual_mean_temp)
#rg_pp <- range(transitionData[which(transitionData$state1==state),]$annual_pp)
rg_tp <- range(transitionData$annual_mean_temp)
rg_pp <- range(transitionData$annual_pp)

states <- c("T","M","B","R")

## Figure FOUR PANELS BY STATE
for (from in 1:length(states)){

# STATE FROM
state_from = states[from]

# predict on scaled data
pp <- seq(rg_pp[1], rg_pp[2], length.out = 200)
tp <- seq(rg_tp[1], rg_tp[2], length.out = 200)
prob = as.data.frame(predict(get(paste0(state_from,'_multi_full')), newdata = data.frame(expand.grid(annual_mean_temp = tp, annual_pp = pp),time=5 ), type = "probs"))

# unscaled data
pp <- pp  * attr(scale_pp,"scaled:scale") + attr(scale_pp,"scaled:center")
tp <- tp * attr(scale_tp,"scaled:scale") + attr(scale_tp,"scaled:center")

png(paste0("~/Desktop/From",state_from,".png"),res=300, width = 5, height = 5.5, units = 'in')
# Set graphic with 4 panels
par(mfrow=c(2,2),mar=c(2,2,1,1),oma=c(1,1,2,1))

for(to in 1:length(states)){
  # STATE TO
  state_to <- states[to]
  if(state_to %in% names(prob)){
    image(x=pp, y=tp, z = t(matrix(prob[,which(names(prob)==state_to)], ncol = length(pp), nrow = length(tp))),xlab = "", ylab = "", col = pal(100),zlim=c(0,1))
    contour(x=pp, y=tp, z = t(matrix(prob[,which(names(prob)==state_to)], ncol = length(pp), nrow = length(tp))),cex=0.5,add=TRUE)
  title(main=paste("to",state_to),cex=0.8,font.main = 2)
 }
}
mtext(paste("From",state_from),side=3,outer=TRUE,cex=1.2,font=2)
mtext("Temperature",side=2,outer=TRUE,cex=1)
mtext("Precipitation",side=1,outer=TRUE,cex=1)
dev.off()
}

# Figure INDIVIDUAL PANEL

### Begin function
fig_prob <- function(state_from,state_to){
  # predict on scaled data
  pp <- seq(rg_pp[1], rg_pp[2], length.out = 200)
  tp <- seq(rg_tp[1], rg_tp[2], length.out = 200)
  prob = as.data.frame(predict(get(paste0(state_from,'_multi_full')), newdata = data.frame(expand.grid(annual_mean_temp = tp, annual_pp = pp),time=5 ), type = "probs"))

  # unscaled data
  pp <- pp  * attr(scale_pp,"scaled:scale") + attr(scale_pp,"scaled:center")
  tp <- tp * attr(scale_tp,"scaled:scale") + attr(scale_tp,"scaled:center")

  image(x=pp, y=tp, z = t(matrix(prob[,which(names(prob)==state_to)], ncol = length(pp), nrow = length(tp))),xlab = NA, ylab = NA, col = pal(100),zlim=c(0,1), cex.axis=1)
  contour(x=pp, y=tp, z = t(matrix(prob[,which(names(prob)==state_to)], ncol = length(pp), nrow = length(tp))),add=TRUE, labcex=0.8)
}

# Multinom FIG1
figs_state <- matrix(c('B','R','R','T','B','M','M','T'),ncol=2,nrow=4,byrow=TRUE)
for(i in 1:nrow(figs_state)){
  tiff(paste0("~/Documents/Git/stm/STModel-CompAnalysis/",paste0(figs_state[i,],collapse='-'),".tiff"),res=300, width = 5, height = 5.5, units = 'in',bg =   "transparent",type='cairo')
  fig_prob(figs_state[i,1],figs_state[i,2])
  dev.off()
}

# Multinom SI_1
figs_state <- matrix(c('T','R','R','B','T','M','M','B'),ncol=2,nrow=4,byrow=TRUE)
for(i in 1:nrow(figs_state)){
  tiff(paste0("~/Documents/Git/stm/STModel-CompAnalysis/",paste0(figs_state[i,],collapse='-'),".tiff"),res=300, width = 5, height = 5.5, units = 'in',bg =   "transparent",type='cairo')
  fig_prob(figs_state[i,1],figs_state[i,2])
  dev.off()
}
