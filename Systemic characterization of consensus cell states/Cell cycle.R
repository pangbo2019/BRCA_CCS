# 2. cell cycle
#' @param maligant.seurat seurat object of malignant cell

library(Seurat)
g2m_genes = cc.genes$g2m.genes
g2m_genes = CaseMatch(search = g2m_genes, match = rownames(maligant.seurat))
s_genes = cc.genes$s.genes
s_genes = CaseMatch(search = s_genes, match = rownames(maligant.seurat))
maligant.seurat <- CellCycleScoring(object = maligant.seurat,  g2m.features = g2m_genes, s.features = s_genes)


