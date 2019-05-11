 library("data.table")

models <- data.table(readRDS("5-Validation/modComparisons_plotsMeasOnce.rds"))
models$STM_R <- 1 - rowSums(models[,9:11])
models$stmState <- c('B', 'T', 'M', 'R')[apply(models[,c(9:11, 20)], 1, which.max)]
models$rf2State <- c('B', 'M', 'R', 'T')[apply(models[,12:15], 1, which.max)]
models$rf7State <- c('B', 'M', 'R', 'T')[apply(models[,16:19], 1, which.max)]

confMatrices <- list(stm = rbind(table(models$stmState, models$state), R = rep(0,4)),
	rf2 = table(models$rf2State, models$state),
	rf7 = table(models$rf7State, models$state))
confMatrices[[1]] <- confMatrices[[1]][colnames(confMatrices[[1]]), ]

## accuracy
acc <- lapply(confMatrices, function(x) sum(diag(x))/ sum(x))
accByType <- lapply(confMatrices, function(y) 
	sapply(1:4, function(x) (y[x,x] + sum(y[-x, -x]))/sum(y)))

## sensitivity
sens <- lapply(confMatrices, function(y) diag(y) / colSums(y))
## specificity
spec <- lapply(confMatrices, function(y) 
	sapply(1:4, function(x) sum(y[-x, -x])/sum(y[,-x])))
tss <- mapply(function(x,y) apply(rbind(x,y), 2, function(z) sum(z) - 1), sens, spec)


# do the same for the calib data
modelsC <- data.table(readRDS("5-Validation/modComparisons.rds"))
modelsC$STM_R <- 1 - rowSums(modelsC[,9:11])
modelsC$stmState <- c('B', 'T', 'M', 'R')[apply(modelsC[,c(9:11, 20)], 1, which.max)]
modelsC$rf2State <- c('B', 'M', 'R', 'T')[apply(modelsC[,12:15], 1, which.max)]
modelsC$rf7State <- c('B', 'M', 'R', 'T')[apply(modelsC[,16:19], 1, which.max)]

confMatricesC <- list(stm = rbind(table(modelsC$stmState, modelsC$state), R = rep(0,4)),
	rf2 = table(modelsC$rf2State, modelsC$state),
	rf7 = table(modelsC$rf7State, modelsC$state))
confMatricesC[[1]] <- confMatricesC[[1]][colnames(confMatricesC[[1]]), ]

## accuracy
acc_c <- lapply(confMatricesC, function(x) sum(diag(x))/ sum(x))
accByType_c <- lapply(confMatricesC, function(y) 
	sapply(1:4, function(x) (y[x,x] + sum(y[-x, -x]))/sum(y)))

## sensitivity
sens_c <- lapply(confMatricesC, function(y) diag(y) / colSums(y))
## specificity
spec_c <- lapply(confMatricesC, function(y) 
	sapply(1:4, function(x) sum(y[-x, -x])/sum(y[,-x])))
tss_c <- mapply(function(x,y) apply(rbind(x,y), 2, function(z) sum(z) - 1), sens_c, spec_c)
