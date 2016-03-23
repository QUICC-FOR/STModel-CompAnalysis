# Get metrics from all simulations
rm(list=ls())

# Setwd
setwd("~/research/STModel-CompAnalysis")
source('./fcts/fcts.r')
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

# We replicated the addresses nCores times.
registerDoParallel(cl)
########## End Mammouth Section

ext_path <- '/mnt/parallel_scratch_mp2_wipe_on_august_2016/dgravel/sviss/'
#ext_path <- './'
#2. Load the four (STM-Global, local, SDM, Solved) into Brick (Not the same resolution but should have the same extent)!
# Create archi folder
com_files <- list.files(paste0(ext_path,"stmodel-global/"),recursive=TRUE,pattern="2095_land_rs.robj")

res <- foreach(file=1:length(com_files),.packages=c('raster','reshape2'),.combine='rbind') %dopar% {
    #1. Load the four rasters and add the stm solved analyticaly
    # Get metadata
    md <- unlist(strsplit(com_files[file],"/"))
    md <- md[which(md == "GenSA_rf_all_2_5y_rep1"):length(md)]
    names(md) <- c("pars","landCI","gcm","rep","year")
    md['year'] <- gsub("[^0-9]", "", md['year'])
    if(md['year']=="0000") md['year'] <- "2000"
    

    
    # stmodel-global
    load(paste0(ext_path,"stmodel-global/",com_files[file]))
    stm_global <- rs

    # stmodel-local
    load(paste0(ext_path,"stmodel-local/",com_files[file]))
    stm_local <- rs

    # SDM
    load(paste0(ext_path,"stmodel-sdm/",md['gcm'],"/",md['year'],"_rf_sdm.robj"))
    sdm <- rs_class

    # stm-1000
    if(md['year']=="2095" & file.exists(paste0("stmodel-10000/",gsub("_land_rs","_land_eq",com_files[file])))){
      load(paste0(ext_path,"stmodel-10000/",gsub("_land_rs","_land_eq",com_files[file])))
      stm_10000 <- rs
    }

    # stm-solved
    load(paste0(ext_path,"stmodel-solved/",md['pars'],"/",md['gcm'],"/",md['year'],"_land_rs.robj"))
    stm_solved <- rs_solved

    # Clean ws
    rm(rs,rs_class,stack_prob,rs_solved,stack_vars)

    #####################################
    # PREP
    #####################################

    # Resample
    stm_global<-resample(stm_global, stm_local, method="ngb")
    stm_solved<-resample(stm_solved, stm_local, method="ngb")
    sdm<-resample(sdm, stm_local, method="ngb")
    if(exists("stm_10000")) stm_10000<-resample(stm_10000, stm_local, method="ngb")

    # Stack rasters
    stm_st_wgs84 <- stack(stm_local,stm_global,stm_solved,sdm)
    if(exists("stm_10000")) stm_st_wgs84 <- addLayer(stm_st_wgs84,stm_10000)
    names(stm_st_wgs84)[1:4]<-c('local','global','solved','sdm')
    if(exists("stm_10000")) names(stm_st_wgs84)[5]<-'t10000'

    # NA Management
    stm_st_wgs84[stm_st_wgs84[]==0] <- NA

    # Reproject
    new_proj <- "+proj=lcc +lat_1=60 +lat_2=46 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
    stm_st_lcc <- projectRaster(stm_st_wgs84,crs=new_proj,method='ngb')

    #####################################
    # BAND
    #####################################

    # Extract surface covered by each state
    # Crop
    ext_study_area <- c(-74.6,-74,45.6,50.3)
    stm_crop_wgs84 <- crop(stm_st_wgs84,ext_study_area)
    new_proj <- paste0("+proj=lcc +lat_1=",extent(stm_crop_wgs84)[4]," +lat_2=",extent(stm_crop_wgs84)[3]," +lat_0=44 +lon_0=",mean(extent(stm_crop_wgs84)[1:2])," +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
    stm_crop_lcc <- projectRaster(stm_crop_wgs84,crs=new_proj,method='ngb')

    # Lat metrics
    stm_df <- as.data.frame(stm_crop_lcc,xy=TRUE,na.rm=TRUE)
    stm_df$x <- as.factor(stm_df$x)
    stm_df <- melt(stm_df,id.vars=c("y","x"),variable.name="model",value.name="state")
    stm_df$state <- as.factor(idToState(stm_df$state))
    stm_df$y <- stm_df$y - min(stm_df$y)


    band_prop <- dcast(y+model~state,data=stm_df, fun.aggregate = length)
    # Add md
    md_out <- as.data.frame(matrix(rep(md,each=nrow(band_prop)),nrow=nrow(band_prop)))
    names(md_out) <- names(md)
    out <- cbind(md_out,band_prop)

    return(out)
}

save(res,file="~/research/STModel-CompAnalysis/out/4-LandMig/bandprop.rdata")
