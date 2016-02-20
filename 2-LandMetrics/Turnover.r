# Purpose:	Create map showing the rate of turnover (states transitons) for each time step.

# Clean workspace
rm(list=ls())

# load library
library('doParallel')


# Get node addresses given by Torque and assigned to your task
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile,header=FALSE) # Node adresses

# Get the number of cores by node
nCores <- 24

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
dirs <- list.dirs("~/research/STModel-CompAnalysis/stmodel-local")
dirs <- grep("rep_",dirs,value=TRUE)

# list of files in the folder
files <- paste0(c('0000',as.character(seq(2000,2095,5))),"_land_rs.robj")

# This one can be parralelize
res <- foreach(dir=1:nrow(),.packages=c('raster')) %dopar% {
  	load(paste(dirs[dir],files[1],sep="/"))
    rs[rs==0] <- NA
    rst1 <- rs
    load(paste(dirs[dir],files[21],sep="/"))
    rs[rs==0] <- NA
    rst0 <- rs

    MtoT <- overlay(rst1, rst0, fun = function(t1,t0) { ifelse( t1 == stateToId('T') &  t0 == stateToId('M') , 1, 0) })
		BtoM <- overlay(rst1, rst0, fun = function(t1,t0) { ifelse( t1 == stateToId('M') &  t0 == stateToId('B') , 1, 0) })
		BTMtoR <- overlay(rst1, rst0, fun = function(t1,t0) { ifelse( t1 == stateToId('R')  &  (t0 == stateToId('T') | t0 == stateToId('B') | t0 == stateToId('M')), 1, 0) })

		folder_out <- paste0('~/research/STModel-CompAnalysis/out/',paste(unlist(strsplit(dirs[dir],'/'))[7:9],collapse="-"))
		dir.create(folder_out, showWarnings = FALSE, recursive = TRUE)

		save(MtoT,BtoM,BTMtoR,file=paste(out_folder,'overlay_rs.rdata',sep='/'))
		return(1)
}

stopCluster(cl)
