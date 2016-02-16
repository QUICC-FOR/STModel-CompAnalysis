
rm(list = ls())

setwd("~/Desktop/glm_explo_stm/")

# Open data
load("~/Documents/GitHub/STModel-Data-4S/out_files/transitions_r5.RData")
# Rename and columns
transitionData$annual_pp = transitionData$tot_annual_pp
transitionData$st0 = transitionData$state1
transitionData$st1 = transitionData$state2

pair.dat <- transitionData

# subset 10 degree
select = unique(pair.dat$plot[which(pair.dat$annual_mean_temp<=10)])
pair.dat_subset10 = pair.dat[pair.dat$plot %in% select,]

#plot(pair.dat_subset10$annual_mean_temp, pair.dat_subset10$tot_annual_pp, cex =.2, xlab = "temperature", ylab ="precipitations")
#title("Calibration plots")


# Create transition column
pair.dat_subset10$transition <- paste(pair.dat_subset10$st0,pair.dat_subset10$st1,sep="")

# Datset without filters
pair_dat0 <- pair.dat_subset10

# Add delta t
pair_dat0$delta_t <-pair_dat0$year2-pair_dat0$year1

# Graph lim
rg_pp <- range(pair.dat$annual_pp)
rg_tp <- range(pair.dat$annual_mean_temp)



###################################################################
#####    Analyses     GLM                                      #######
###################################################################
library(MASS)
library(ROCR)
library(fmsb)
library(colorRamps)

pal <-colorRampPalette(rev(brewer.pal(11,"RdYlBu")))

#####    glm climate                            #######

modelTransition <- function(st0 , st1, pair.dat, name = NULL)
{
  print(st0)
  print("->")
  print(st1)
  
  datst0 = pair.dat[pair.dat$st0%in%st0, ]
  datst0$transition = ifelse(datst0$st1 %in% st1, 1, 0)
  
  n <- dim(datst0)[1]
  
  mod = glm(transition ~ annual_mean_temp + I(scale(annual_mean_temp)^2) + I(scale(annual_mean_temp)^3) + annual_pp + I(scale(annual_pp)^2) + I(scale(annual_pp)^3) + delta_t + annual_mean_temp:annual_pp, family = "binomial", data =datst0)
  stepMod  = stepAIC(mod)
  #print(summary(stepMod))
  
  
  pred = predict(stepMod,new=datst0,"response")
  # overall performance
  R2 = NagelkerkeR2(stepMod)$R2
  #discrimination
  perf = performance(prediction(pred, datst0$transition), "auc")
  AUC = perf@y.values[[1]]
  
  ## selected vars
  ## selected vars
  coeff = summary(stepMod)$coefficients
  vars = rownames(coeff)[-1]
  effect = coeff[-1,1]
  pval = coeff[-1,4]
  print(pval)
  
  if(length(st0)>2) st0 <- paste(st0,collapse = ",")
  if(is.null(name)) name = paste(st0, st1, sep = "->")
  
  return(list(mod = stepMod,n=n, vars = vars,effect = effect, pval = pval , R2 = R2, AUC = AUC, ranges = apply(datst0[unlist(lapply(1:ncol(datst0), function(x)is.numeric(datst0[,x])))], 2, range), name = name))
}

#test
#mod = modelTransition_climate("R", c("T", "M"), pair.dat = pair_dat0)


###################################################################
#####    figures                            #######
###################################################################


fig_glm <- function(mod)
{
  
  temp = seq(as.numeric(mod$ranges[1, "annual_mean_temp"]), as.numeric(mod$ranges[2, "annual_mean_temp"]), length.out = 50)
  pp = seq(as.numeric(mod$ranges[1, "annual_pp"]), as.numeric(mod$ranges[2, "annual_pp"]), length.out = 50)
  prob = predict(mod$mod, newdata = data.frame(expand.grid(annual_mean_temp = temp, annual_pp = pp),delta_t=5 ), type = "response")
  
  image(x=temp, y=pp, z = matrix(prob, ncol = length(pp), nrow = length(temp)),xlab = "Temperature", ylab = "Precipitations", col = pal(100), main = mod$name,zlim=c(0,1))
  contour(x=temp, y=pp, z = matrix(prob, ncol = length(pp), nrow = length(temp)),cex=0.5,add=TRUE)
}

#test
#fig_glm(mod)

###################################################################
#####    figures  completes                          #######
###################################################################


#models
# regeneration
modRT = modelTransition("R", "T", pair.dat = pair_dat0)
modRB = modelTransition("R", "B", pair.dat = pair_dat0)
modRM = modelTransition("R", "M", pair.dat = pair_dat0)
# exclusion
modMT = modelTransition(c("M"), c("T"), pair.dat = pair_dat0)
modMB = modelTransition(c("M"), c("B"), pair.dat = pair_dat0)
# colonisation
modTM = modelTransition(c("T"), c("M"), pair.dat = pair_dat0)
modBM = modelTransition(c("B"), c("M"), pair.dat = pair_dat0)
# disturbance
modR = modelTransition(c("T", "B", "M"), c("R"), pair.dat = pair_dat0)
modMR = modelTransition(c("M"), c("R"), pair.dat = pair_dat0)
modTR = modelTransition("T", c("R"), pair.dat = pair_dat0)
modBR = modelTransition(c("B"), c("R"), pair.dat = pair_dat0)


png("glm_transitions_dom.png", width = 10, height = 10,units = "in",res=300)

layout(matrix(c(1, 2, 3, 4, 5, 6 ), ncol = 2))
par(mar=c(2,2,3,2),oma=c(1,1,0,0),cex.main=1.5)
fig_glm(modRT)
fig_glm(modRB)
fig_glm(modRM)


fig_glm(modTM)
fig_glm(modBM)
fig_glm(modR)


dev.off()

# Prepare table for publication

mod <- list(RT = modRT, RB = modRB, RM = modRM, MT = modMT, MB = modMB, TM = modTM, BM = modBM, R = modR, TR = modTR, BR = modBR, MR = modMR)

pval.star <- function(pval)
{
  star=""
  
  if(pval<0.1) star = "."
  if(pval<0.05) star = "*"
  if(pval<0.01) star= "**"
  if(pval<0.001) star = "***"
  
  return(star)
}

library(xtable)
library(plyr)

out <- list()
df <- data.frame()

for (i in 1:length(mod)){
  
  stars <- character()
  pvals <- mod[[i]]$pval

  for(pos in 1:length(pvals)) stars <- append(stars,pval.star(pvals[pos]))
  
  effects <- paste(round(mod[[i]]$effect,4),stars,sep="")
  names(effects) <- names(mod[[i]]$effect)
  
  if(length(mod[[i]]$name)>1) mod[[i]]$name <- "T,M,B->R"
  
  r <- data.frame(name=mod[[i]]$name,R2=round(mod[[i]]$R2,2),AUC=round(mod[[i]]$AUC,2),n=mod[[i]]$n)
  ef <- data.frame(t(effects))
  df <- cbind(r,ef)
  out[[i]] <- df
}

df <- do.call("rbind.fill",out)
df <- df[,c(1,4:9,2,3)]

# Ajouter le n par transition
