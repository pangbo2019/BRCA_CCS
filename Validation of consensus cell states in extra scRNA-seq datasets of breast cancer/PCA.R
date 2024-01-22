#' Principal components analysis
#' @param patient.mean.mat gene average expression matrix for patients, rows are genes and columns are patients
#' @param subtype a vector of patient subtypes

library(stats)
pca <- prcomp(t(patient.mean.mat),scal = T)
pca_sum <- summary(pca)
pc <- pca_sum$importance[2,]*100
pca12 <- as.data.frame(pca_sum$x[,c(1,2)])
pca12 <- cbind.data.frame(sample = colnames(patient.mean.mat), pca12,
                          subtype = subtype)

