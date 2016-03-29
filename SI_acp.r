## PCA on composition of the plot and climatic variables

# Classify as State
library('reshape2')

tSpecies = c('19481-BET-ALL','19408-QUE-RUB','28728-ACE-RUB','28731-ACE-SAC',
  '32931-FRA-AME','32945-FRA-NIG','19462-FAG-GRA','19511-OST-VIR','21536-TIL-AME',
  '24764-PRU-SER')
bSpecies = c('183302-PIC-MAR','183295-PIC-GLA','18034-PIC-RUB','183412-LAR-LAR',
  '183319-PIN-BAN','18032-ABI-BAL')
enVar = c("annual_mean_temp", "tot_annual_pp", "mean_diurnal_range","ph_2cm", "slp")

r = 1

# read the data
treeDat = read.csv('./extra/treeData.csv')
plotDat = read.csv('./extra/plotInfoData.csv')
climDat = read.csv('./extra/climData.csv')
treeDat = merge(treeDat, plotDat, all.x=TRUE)
treeDat$type = rep('U', nrow(treeDat))
treeDat$type[which(treeDat$id_spe %in% tSpecies)] = 'T'
treeDat$type[which(treeDat$id_spe %in% bSpecies)] = 'B'

# get rid of all plots that NEVER have at least one T or B species
# this prevents them from being classified as R when they never contain even one
# species of interest
trTab = table(treeDat$plot_id, treeDat$type)
filterNames = as.numeric(rownames(trTab[rowSums(trTab[,1:2]) == 0,]))
treeDat.filtered = treeDat[!(treeDat$plot_id %in% filterNames),]

# reshape the data into plot-year samples by state
sampleDat = dcast(treeDat.filtered, plot_id + year_measured + lat + lon + id_spe ~ type, fill = 0, 
		value.var = "basal_area", fun.aggregate = sum)
sampleDat$sumBA = sampleDat$B + sampleDat$T + sampleDat$U
sampleDat$state = rep('U', nrow(sampleDat))
sampleDat$state[sampleDat$B > 0 & sampleDat$T == 0] = 'B'
sampleDat$state[sampleDat$B == 0 & sampleDat$T > 0] = 'T'
sampleDat$state[sampleDat$B > 0 & sampleDat$T > 0] = 'M'
sampleDat$state[sampleDat$sumBA < r] = 'R'

# Drop species row which is not in the list of species studied
sampleDat <- sampleDat[sampleDat$id_spe %in% c(tSpecies,bSpecies),]
sampleDat$id_spe <- droplevels(sampleDat$id_spe)

sampleDat = dcast(sampleDat, plot_id + year_measured + lat + lon ~ id_spe, value.var="sumBA" , fun.aggregate = sum)

stateData = merge(sampleDat, climDat, all = 'T', by=c("plot_id", "year_measured"))
stateData = stateData[!is.na(stateData$lat), ]
stateData = stateData[!is.na(stateData$mean_diurnal_range), ]


stateData$uq_id <- as.numeric(paste0(stateData$plot_id,stateData$year_measured))
plot_w_env <- stateData[,c("uq_id",enVar)]
plot_w_env <- unique(plot_w_env)
rownames(plot_w_env) <- plot_w_env$uq_id
plot_w_env <- plot_w_env[,-1]
plot_w_com <- stateData[,c("uq_id",tSpecies,bSpecies)]
plot_w_com <- unique(plot_w_com)
rownames(plot_w_com) <- plot_w_com$uq_id
plot_w_com <- plot_w_com[,-1]

library("ade4")
var.pca = dudi.pca(plot_w_env[complete.cases(plot_w_env),], scale=TRUE,scann=FALSE, nf = ncol(plot_w_env))
barplot(var.pca$eig)
kip <- 100 * var.pca$eig/sum(var.pca$eig)
inertia.dudi(var.pca, row=F,col=T)
cumsum(kip)

s.corcircle(var.pca$co, xax = 1, yax = 2)
##
###pdf("../figures/PCA_selectedVariables.pdf")
s.arrow(var.pca$co, clab = .6, xlim = c(-2,2), sub = "axe 1: 43% ; axe 2: 27 %")
#s.arrow(var.pca$co[selectedVars,], clab = .6, xlim = c(-2,2),yax=3, sub = "axe 1: 43% ; axe 3: 9 %")
#dev.off()

#varCor2 = cor(dat_subset10[, selectedVars])
#varCor2

library(devtools)
install_github("vqv/ggbiplot")

library('ggbiplot')
data(USArrests)
data(state)
m <- princomp(USArrests)
df <- fortify(m, scale = 1)

g <- ggplot(df, aes(x = PC1, y = PC2)) +
  geom_text(aes(label = state.abb[match(rownames(df),state.name)])) +
  geom_axis(data = attr(df, "basis"), aes(label = .name))
print(g)

