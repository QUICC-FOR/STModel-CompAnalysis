# Purpose:	Create map showing the rate of turnover (states transitons) for each time step.

# Clean workspace
rm(list=ls())

files <- list.files('~/Documents/Git/stm/STModel-CompAnalysis/out/',full.names=TRUE,recursive=TRUE,pattern='overlay_rs')

# Create ref raster
load(files[1])
ref_rs <- MtoT
ref_rs[!is.na(ref_rs[])] <- 0
rs_MtoT <- rs_BtoM <- rs_BTMtoR <- ref_rs

# Sum over files
for(file in 1:length(files)){
  load(files[file])
  rs_MtoT <- rs_MtoT+MtoT
  rs_BtoM <- rs_BtoM+BtoM
  rs_BTMtoR <- rs_BTMtoR+BTMtoR
}

# Get the mean
MtoT <- rs_MtoT/length(files)
BtoM <- rs_BtoM/length(files)
BTMtoR <- rs_BTMtoR/length(files)

# Create stack
st_transitions <- stack(MtoT,BtoM)
names(st_transitions) <- c('MtoT','BtoM')
st_transitions[st_transitions[]<0] <- 0

# Load supp libraries
library('rgdal')
library('RColorBrewer')
library('viridis')

# Create Maps
st_transitions <- projectRaster(st_transitions,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
load("./shp/altBorder.robj")
load("./shp/shp_stm_area.rdata")
st_transitions <- mask(st_transitions,border)

new_proj <- "+proj=lcc +lat_1=60 +lat_2=46 +lat_0=44 +lon_0=-68.5 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
st_transitions <- projectRaster(st_transitions, crs=new_proj)
border <- spTransform(border,CRS(new_proj))
great_lakes <- spTransform(great_lakes,CRS(new_proj))

# Set params for figures
pal <-colorRampPalette(rev(brewer.pal(11,"RdYlBu")))
brk <- 100

#png(file="./future_distrib_states.pdf",width = 8, height =4,units="in",res=300)
#pdf(file="./future_distrib_states.pdf",width = 8, height =4)
tiff(file="./future_distrib_states.tiff",width = 8, height =4,units="in",res=300)
layout(matrix(c(1:2,3,3), 2,2 ,byrow=TRUE),height=c(1,0.4),TRUE)
par(mar=c(0.5,0.5,0.5,0.5),oma=c(0, 2, 4, 0),xpd=NA)

plot(border,border=NA,col="grey80")
image(st_transitions$BtoM,axes=FALSE,xlab="",ylab=,asp=1,breaks=round(seq(0,0.60,length.out=brk),2),col=pal(brk-1),add=TRUE)
plot(great_lakes,lwd=0.4,col="white",border="white",add=TRUE)
plot(border,add=TRUE,lwd=0.6,border="white",col=NA)
llgridlines(border,cex=0.5,lty=3,col='grey50')

plot(border,border=NA,col="grey90")
image(st_transitions$MtoT,axes=FALSE,xlab="",ylab="",asp=1,zlim=c(0,1),breaks=round(seq(0,0.60,length.out=brk),2),col=pal(brk-1),add=TRUE)
plot(great_lakes,lwd=0.4,col="white",border="white",add=TRUE)
plot(border,add=TRUE,lwd=0.6,border="white",col=NA)
llgridlines(border,cex=0.5,lty=3,col='grey50')

labs<-as.character(round(seq(0,0.60,length.out=5),2))
labs[length(labs)] <- paste0(">",labs[length(labs)])

par(mar=c(3,20,2,2))
image(matrix(1:brk-1),col=pal(brk-1),axes=FALSE,ann=FALSE)
axis(1,at=seq(0,1,length.out=5),lab=labs,col='grey40',col.ticks='grey40',col.axis="grey40", cex.axis=0.65)
box(lwd=0.8,col='grey40')

mtext("Proportion of transitions observed",1,line=-0.75,at=-0.30,font=2,cex=0.65,col='grey20')

dev.off()
