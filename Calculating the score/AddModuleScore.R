# AddModuleScore
#' @param object maligant.seurat: seurat object of malignant cell
#' @param features GeneSets: list of genesets
#' @param name names(GeneSets)

library(Seurat)
score <- AddModuleScore(object = maligant.seurat, features = GeneSets, name = names(GeneSets))


