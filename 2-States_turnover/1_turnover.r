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
dirs <- list.dirs("~/research/STModel-CompAnalysis/stmodel-local",pattern="_rep")

# list of files in the folder
files <- paste0(c('0000',as.character(seq(2000,2095,5))),"_land_rs.robj")

# This one can be parralelize
res <- foreach(dir=1:length(dirs),.packages=c('raster')) %dopar% {
  	load(paste(dirs[dir],files[1],sep="/"))
    rs[rs==0] <- NA
    rst1 <- rs
    load(paste(dirs[dir],files[length(files)],sep="/"))
    rs[rs==0] <- NA
    rst0 <- rs

      MtoT <- overlay(rst1, rst0, fun = function(t1,t0) { ifelse( t1 == stateToId('T') &  t0 == stateToId('M') , 1, 0) })
		BtoM <- overlay(rst1, rst0, fun = function(t1,t0) { ifelse( t1 == stateToId('M') &  t0 == stateToId('B') , 1, 0) })
		BTMtoR <- overlay(rst1, rst0, fun = function(t1,t0) { ifelse( t1 == stateToId('R')  &  (t0 == stateToId('T') | t0 == stateToId('B') | t0 == stateToId('M')), 1, 0) })

		folder_out <- paste0('/mnt/parallel_scratch_mp2_wipe_on_august_2016/dgravel/sviss/STModel-CompAnalysis/out/',paste(unlist(strsplit(dirs[dir],'/'))[7:9],collapse="-"))
		dir.create(folder_out, showWarnings = FALSE, recursive = TRUE)

		save(MtoT,BtoM,BTMtoR,file=paste(out_folder,'overlay_rs.rdata',sep='/'))
		return(1)
}

stopCluster(cl)
