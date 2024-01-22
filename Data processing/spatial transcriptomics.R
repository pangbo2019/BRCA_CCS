# ST

# Data preprocessing
#' @param spe raw seurat object
#' @param n_spots genes expressed at least in how many spots
#' @param n_feature spots expressed at least in how many genes
#' @param percent_mito the largest proportion of mitochondrial gene sets

library(Seurat)
ST_preprogress <- function(spe, n_spots = 2, n_feature = 200, percent_mito = 20){
  mt.genes <- grep(pattern = "^MT-", x = rownames(spe), value = TRUE)
  spe$percent.mito <- (Matrix::colSums(spe[["Spatial"]]$counts[mt.genes, ])/Matrix::colSums(spe[["Spatial"]]$counts))*100
  genes_to_keep <- setdiff(names(which(Matrix::rowSums(spe[["Spatial"]]$counts != 0)>n_spots)),mt.genes)
  spe <- subset(spe,
                features = genes_to_keep, 
                subset = nFeature_Spatial > n_feature & percent.mito < percent_mito
  )
  return(spe)
}

# SCTransform
#' @param spe seurat object after running ST_preprogress function
#' @param assay name of assay
#' @param variable.features.n number of HVGs

ST_normalize <- function(spe, assay = "Spatial", variable.features.n = 3000){
  spe <- SCTransform(spe, assay = assay,
                     variable.features.n = variable.features.n,
                     vars.to.regress = NULL,
                     verbose = FALSE)
}

