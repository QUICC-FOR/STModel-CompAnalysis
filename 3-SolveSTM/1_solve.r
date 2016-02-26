# Solve STM for each clim projected
rm(list=ls())

# Setwd
setwd("~/research/STModel-CompAnalysis")

# Load library and extra objects
source('./fcts/fcts_resAnalytic.r')
source('./fcts/fcts.r')
load('./extra/rsRef_grid.robj')
pars <- read.table("./pars/GenSA_rf_all_2_5y_rep1.txt",row.names=1)
########## Mammouth Section
# load library
library('doParallel')

# Get node addresses given by Torque and assigned to your task
nodefile <- Sys.getenv("PBS_NODEFILE")
hostlist <- read.table(nodefile,header=FALSE) # Node adresses

# Get the number of cores by node
nCores <- detectCores()

# Open the cluster
cl <- makeCluster(rep(c(as.character(hostlist$V1)),nCores),type='SOCK')
#cl <- makeCluster(1)

# We replicated the addresses nCores times.
registerDoParallel(cl)
########## End Mammouth Section

files <- list.files("~/research/STModel-CompAnalysis/in/in_stm_solved",recursive=TRUE,full.names=TRUE,pattern=".csv")

res <- foreach(file=1:length(files),.packages=c('rootSolve','raster','stringr')) %dopar% {
  clim <- read.csv(files[file])[,c('lon','lat','annual_mean_temp','tot_annual_pp')]
  names(clim)[3:4] <- c('tp','pp')
  stm_solved <- solve_stm(clim,pars)
  stm_solved$state <- names(stm_solved[,c('B','T','M')])[apply(stm_solved[,c('B','T','M')],1,which.max)]
  stm_solved$state <- stateToId(stm_solved$state)
  rs_solved <- stm_solved[,c('lon','lat','state')]

  coordinates(rs_solved) <- ~lon + lat
  gridded(rs_solved) <- TRUE
  projection(rs_solved) <- CRS("+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs")
  rs_solved <- raster(rs_solved)
  rs_solved <- crop(rs_solved,extent(ref_rs))

  # get id GCM
  md <- unlist(strsplit(files[file],"/"))
  gcm_id <- md[8]
  year<- unlist(str_extract_all(md[9],"[0-9]+"))[1]

  # Save object
  dir.create(paste0("/mnt/parallel_scratch_mp2_wipe_on_august_2016/dgravel/sviss/stmodel-solved/GenSA_rf_all_2_5y_rep1/",gcm_id,"/"),showWarnings=FALSE,recursive=TRUE)
  save(rs_solved,file=paste0("/mnt/parallel_scratch_mp2_wipe_on_august_2016/dgravel/sviss/stmodel-solved/GenSA_rf_all_2_5y_rep1/",gcm_id,"/",year,"_land_rs.robj"))
}

stopCluster(cl)
