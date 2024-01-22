# 1. singleR 
#' after clustering
#' @param processed.seurat seurat object

library(Seurat)
library(celldex)
library(SingleR)
hpca.se <- celldex::HumanPrimaryCellAtlasData()
processed.seurat_SingleR <- GetAssayData(processed.seurat, slot="data") 
hesc <- SingleR(test = processed.seurat_SingleR, ref = hpca.se, labels = hpca.se$label.main)
processed.seurat$celltype <- hesc$labels


