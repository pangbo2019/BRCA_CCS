# 1. Differential analysis
#' @param maligant.seurat seurat object of malignant cell

library(Seurat)
levels(maligant.seurat)
markers <- FindAllMarkers(maligant.seurat, group.by = 'seurat_clusters', only.pos = TRUE, min.pct = 0.1, logfc.threshold = 1)
markers <- markers[markers$p_val_adj<0.05,] 

##DoHeatmap
maligant.seurat <- ScaleData(object = maligant.seurat, features = rownames(maligant.seurat))

pdf("differentially.expressed.genes.pdf",height = 10,width = 14)
DoHeatmap(maligant.seurat, features = markers$gene)
dev.off()


