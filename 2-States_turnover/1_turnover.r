# Purpose:	Create map showing the rate of turnover (states transitons) for each time step.

# Clean workspace
rm(list=ls())

# load library
library('doParallel')
source("~/research/STModel-CompAnalysis/fcts/fcts.r")


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
# list all folders
dirs <- list.dirs("/mnt/parallel_scratch_mp2_wipe_on_august_2016/dgravel/sviss/stmodel-local")
dirs <- dirs[grep('rep_',dirs)]

# This one can be parralelize
res <- foreach(dir=1:length(dirs),.packages=c('raster')) %dopar% {
  	load(paste(dirs[dir],'2015_land_rs.robj',sep="/"))
    rs[rs==0] <- NA
    rst0 <- rs
    load(paste(dirs[dir],'2095_land_rs.robj',sep="/"))
    rs[rs==0] <- NA
    rst1 <- rs

    MtoT <- overlay(rst1, rst0, fun = function(t1,t0) { ifelse( t1 == stateToId('T') &  t0 == stateToId('M') , 1, 0) })
		BtoM <- overlay(rst1, rst0, fun = function(t1,t0) { ifelse( t1 == stateToId('M') &  t0 == stateToId('B') , 1, 0) })
		#TMtoR <- overlay(rst1, rst0, fun = function(t1,t0) { ifelse( t1 == stateToId('R')  &  (t0 == stateToId('T') | t0 == stateToId('B') | t0 == stateToId('M')), 1, 0) })

		folder_out <- paste0('/mnt/parallel_scratch_mp2_wipe_on_august_2016/dgravel/sviss/STModel-CompAnalysis/out/',paste(unlist(strsplit(dirs[dir],'/'))[7:9],collapse="-"))
		dir.create(folder_out, showWarnings = FALSE, recursive = TRUE)

		save(MtoT,BtoM,file=paste(folder_out,'overlay_rs.rdata',sep='/'))
		return(0)
}

stopCluster(cl)
