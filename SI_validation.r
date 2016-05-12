## Create map with equilibrium distrib

# Clean workspace
rm(list=ls())

# load library
library('doParallel')
source("./fcts/fcts.r")

# Open the cluster
cl <- makeCluster(24)

# We replicated the addresses nCores times.
registerDoParallel(cl)

ls_files <- list.files('./stmodel-local',full.names=TRUE,recursive=TRUE)
ls_files <- ls_files[grep('2015',ls_files)]
ls_files <- ls_files[grep('landCI_1',ls_files)]

# State_t
ls_state <- foreach(i=1:length(ls_files),.packages=c('raster')) %dopar% {

  load(ls_files[i])

  # get only T
  int_state <- 'T'
  id_state <- stateToId(int_state)
  rs[rs == 0] <- NA
  rs[rs != id_state] <- 0
  rs[rs == id_state] <- 1
  
  return(rs)

}

ls_state <- stack(ls_state)
state_t <- sum(ls_state)/length(ls_files)

# State_m
ls_state <- foreach(i=1:length(ls_files),.packages=c('raster')) %dopar% {

  load(ls_files[i])

  # get only M
  int_state <- 'M'
  id_state <- stateToId(int_state)
  rs[rs == 0] <- NA
  rs[rs != id_state] <- 0
  rs[rs == id_state] <- 1
  
  return(rs)

}

ls_state <- stack(ls_state)
state_m <- sum(ls_state)/length(ls_files)

# Sum both raster
state_tm <- state_t + state_m

# Load supp libraries
library('rgdal')
library('RColorBrewer')
library('rgeos')

# Load Shapefiles
load("../STModel-Wrapper/inputs/shp/shp_stm_area.rdata")
load("../STModel-Wrapper/inputs/shp/altBorder.robj")
veg <- readOGR(dsn="extra", layer='zv05073g')
veg <- crop(veg,extent(border))
border <- gSimplify(border,tol=0.005,topologyPreserve=TRUE)
veg <- gSimplify(veg,tol=0.005,topologyPreserve=TRUE)

# Project in LCC
new_proj <- "+proj=lcc +lat_1=60 +lat_2=46 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
state_tm <- projectRaster(state_tm, crs=new_proj)
border <- spTransform(border,CRS(new_proj))
state_tm <- mask(state_tm,border)
veg <- spTransform(veg,CRS(new_proj))
great_lakes <- spTransform(great_lakes,CRS(new_proj))

library(colorspace)   ## hsv colorspace manipulations

## Function for desaturating colors by specified proportion
desat <- function(cols, sat=0.5) {
    X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
    hsv(X[1,], X[2,], X[3,])
}

# Set params for figures
pal <-colorRampPalette(rev(brewer.pal(11,"RdYlBu")))
brk <- 100


png(file="~/Documents/GitHub/ms_stm/ms_4states_steve/figs/pseudo_validation.png",width=8,height=7,units="in",res=300)

layout(matrix(c(1,2), 2,1 ,byrow=TRUE),height=c(1,0.2),TRUE)
par(mar=c(0.5,0.5,0.5,0.5),oma=c(0, 0, 0, 0),xpd=TRUE)
plot(border,border=NA,col="grey90")
image(state_tm,axes=FALSE,xlab="",ylab=,asp=1,breaks=round(seq(0,1,length.out=brk),2),col=desat(pal(brk-1),0.3),add=TRUE)
plot(great_lakes,lwd=0.4,col="white",border="white",add=TRUE)
plot(veg[145],add=TRUE,lwd=1.5,border="tomato")
plot(border,add=TRUE,lwd=1.5,border="white",col=NA)
llgridlines(border,cex=0.7,lwd=0.75,lty=3,col='grey60')
legend("bottomright",'Bioclimatic domain of temperate vegetation in Quebec (MFFP, 2016)',lty=1,col=c("tomato"),
lwd=c(2),inset = c(0, -0.02), xpd = TRUE, bty='n', cex=.75)
labs<-as.character(round(seq(0,1,length.out=5),2))
par(mar=c(4,20,1,1),xpd = TRUE)
image(matrix(1:brk-1),col=desat(pal(brk-1),0.3),axes=FALSE,ann=FALSE)
axis(1,at=seq(0,1,length.out=5),lab=labs,col='grey30',
    col.ticks='grey30',col.axis="grey30", cex.axis=0.7)
box(lwd=0.8,col='grey30')
mtext("Occupancy of temperate stands",1,line=-2,font=2,cex=0.75,col='grey20')

dev.off()