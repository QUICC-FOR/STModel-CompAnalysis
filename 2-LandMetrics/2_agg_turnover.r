# Purpose:	Create map showing the rate of turnover (states transitons) for each time step.

# Clean workspace
rm(list=ls())

files <- list.files('/mnt/parallel_scratch_mp2_wipe_on_august_2016/dgravel/sviss/STModel-CompAnalysis/out/',full.names=TRUE,recursive=TRUE,pattern='overlay_rs')

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

# Save
save(MtoT,BtoM,BTMtoR,file='~/research/STModel-CompAnalysis/2-LandMetrics/Transition_rates.rdata')
