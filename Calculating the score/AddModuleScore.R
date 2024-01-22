# AddModuleScore
#' @param maligant.seurat seurat object of malignant cell
#' @param GeneSet list of genesets

library(Seurat)
score <- AddModuleScore(object = maligant.seurat, features = GeneSet, name = names(GeneSet))


