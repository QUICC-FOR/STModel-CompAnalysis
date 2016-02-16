# Purpose:	Create map showing the rate of turnover (states transitons) for each time step.

# Clean workspace
rm(list=ls())

# load library
library('raster')
library('doParallel')

# Archive path
arch <- Sys.getenv("ARCHIVE")

# Get node addresses given by Torque and assigned to your task
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile,header=FALSE) # Node adresses

# Get the number of cores by node
nCores <- detectCores()

# Open the cluster
cl <- makeCluster(rep(c(as.character(hostlist$V1)),nCores),type='SOCK')
# We replicated the addresses nCores times.
registerDoParallel(cl)

################################################
# Some usefull functions
stateToId <- function(state){
	# Function stateToId
	# Purpose: Convert state to id
	state <- as.character(state)
	state[state=="B"] <- 1
	state[state=="T"] <- 2
	state[state=="M"] <- 3
	state[state=="R"] <- 4
	state[is.na(state)] <- 0

	return(as.numeric(state))
}


################################################
# list all folders
dirs <- list.dirs(paste0(arch,"/stm-out/stmodel-local"))
dirs <- grep("rep_",dirs,value=TRUE)

# list of files in the folder
files <- paste0(c('0000',as.character(seq(2000,2095,5))),"_land_rs.robj")

# This one can be parralelized
MtoT <- foreach(dir=1:length(dirs),.packages='raster') %dopar% {
  ls_over <- list()
  for (file in 2:length(files)){
  	# for each time step in the folder
  	load(paste(dirs[dir],files[file],sep="/"))
    rs[rs==0] <- NA
    rst1 <- rs
    load(paste(dirs[dir],files[file-1],sep="/"))
    rs[rs==0] <- NA
    rst0 <- rs

    ls_over[[file-1]] <- overlay(rst1, rst0, fun = function(t1,t0) { ifelse( t1 == stateToId('T') &  t0 == stateToId('M') , 1, 0) })
  }
  sum(stack(ls_over))
}

MtoT <- mean(stack(MtoT))

BtoM <- foreach(dir=1:length(dirs),.packages='raster') %dopar% {
  ls_over <- list()
  for (file in 2:length(files)){
  	# for each time step in the folder
  	load(paste(dirs[dir],files[file],sep="/"))
    rs[rs==0] <- NA
    rst1 <- rs
    load(paste(dirs[dir],files[file-1],sep="/"))
    rs[rs==0] <- NA
    rst0 <- rs

    ls_over[[file-1]] <- overlay(rst1, rst0, fun = function(t1,t0) { ifelse( t1 == stateToId('M') &  t0 == stateToId('B') , 1, 0) })
  }
  sum(stack(ls_over))
}

BtoM <- mean(stack(BtoM))

BTMtoR <- foreach(dir=1:length(dirs),.packages='raster') %dopar% {
  ls_over <- list()
  for (file in 2:length(files)){
  	# for each time step in the folder
  	load(paste(dirs[dir],files[file],sep="/"))
    rs[rs==0] <- NA
    rst1 <- rs
    load(paste(dirs[dir],files[file-1],sep="/"))
    rs[rs==0] <- NA
    rst0 <- rs

    ls_over[[file-1]] <- overlay(rst1, rst0, fun = function(t1,t0) { ifelse( t1 == stateToId('R')  &  (t0 == stateToId('T') | t0 == stateToId('B') | t0 == stateToId('M')), 1, 0) })
  }
  sum(stack(ls_over))
}
BTMtoR <- mean(stack(BTMtoR))

# Save all object
save(BtoM,BTMtoR,MtoT,file="~/STModel-CompAnalysis/out/turnover.rdata")

stopCluster(cl)
