## PCA on composition of the plot and climatic variables

# Classify as State
library('reshape2')

tSpecies = c('19481-BET-ALL','19408-QUE-RUB','28728-ACE-RUB','28731-ACE-SAC',
  '32931-FRA-AME','32945-FRA-NIG','19462-FAG-GRA','19511-OST-VIR','21536-TIL-AME',
  '24764-PRU-SER')
bSpecies = c('183302-PIC-MAR','183295-PIC-GLA','18034-PIC-RUB','183412-LAR-LAR',
  '183319-PIN-BAN','18032-ABI-BAL')
enVar = c("annual_mean_temp", "tot_annual_pp", "mean_diurnal_range","ph_2cm", "slp")

# read the data
treeDat = read.csv('./extra/treeData.csv')
plotDat = read.csv('./extra/plotInfoData.csv')
climDat = read.csv('./extra/climData.csv')

# rm unnecessary columns in climData
climDat <- climDat[,c('plot_id','year_measured',enVar)]

# merge ba with climate based on plot_id
treeDat <- merge(treeDat,climDat,by=c('plot_id','year_measured'),all.x=TRUE)

# apply filter on temperature 
treeDat <- treeDat[treeDat$annual_mean_temp<10,]
treeDat$id_spe <- droplevels(treeDat$id_spe)
treeDat$id <- paste(treeDat$plot_id,treeDat$year_measured,sep='_') 

# applied filter on BA
sum_ba <- aggregate(basal_area ~ id,data=treeDat,FUN=sum)
qu95 <- quantile(sum_ba$basal_area,probs=c(.95))
flt_ba <- sum_ba[sum_ba$basal_area<=qu95,'id']
treeDat <- subset(treeDat, !treeDat$id %in% flt_ba)

# Set descriptors (species)
treeDat <- treeDat[,c('id','id_spe','basal_area')]
mat <- dcast(id ~ id_spe, data=treeDat,fill=0,fun.aggregate=sum)
rownames(mat) <- mat$id
mat <- mat[,-1]


library('vegan')
dist <- vegdist(mat, method="bray")

library("ade4")
var.pca = dudi.pca(mat,scale=TRUE,scannf=TRUE, nf = ncol(mat))
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

pr.pca <- prcomp(mat) 


data(wine)
wine.pca <- prcomp(wine, scale. = TRUE)
ggbiplot(pr.pca) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

fviz_pca_var(var.pca, col.var="contrib") +
scale_color_gradient2(low="white", mid="blue", 
      high="red", midpoint=20) + theme_minimal()

library(devtools)
install_github("vqv/ggbiplot")

library('ggbiplot')
data(USArrests)
data(state)
m <- princomp(USArrests)
df <- fortify(var.pca, scale = 1)

g <- ggplot(df, aes(x = PC1, y = PC2)) +
  geom_text(aes(label = state.abb[match(rownames(df),state.name)])) +
  geom_axis(data = attr(df, "basis"), aes(label = .name))
print(g)

