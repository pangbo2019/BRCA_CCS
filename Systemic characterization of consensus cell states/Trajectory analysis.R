# Trajectory analysis
#' @param obj seurat object of malignant cell
#' @param color  colors to plot
#' @param pc.method pca or harmony
#' @param ex.umap Whether to import seurat umap
#' @param group.by cell groups
#' @param group_label_size the size of the grouping label
#' @param label_cell_groups whether cell grouping is shown
#' @param label_leaves whether to display a leaf node
#' @param label_branch_points whether the branch point is shown

run_monocle3 <- function(obj, pc.method = "pca", ex.umap = T, group.by = 'seurat_clusters', group_label_size = 5,
                         label_cell_groups = F, label_leaves = T, label_branch_points = F, color) {
  library(Seurat)  
  library(monocle3)
  library(ggplot2)
  expression_matrix <- GetAssayData(object = obj, assay = "RNA", slot = "counts")
  cell_metadata <- obj@meta.data
  all(colnames(expression_matrix) == rownames(cell_metadata))
  gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
  rownames(gene_annotation) <- rownames(expression_matrix)
  
  cds <- new_cell_data_set(expression_matrix,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  
  cds@int_colData$reducedDims$PCA <- obj@reductions[[pc.method]]@cell.embeddings 
  
  cds <- reduce_dimension(cds, reduction_method = 'UMAP', preprocess_method = "PCA")
  
  if(ex.umap){
    cds.embed <- cds@int_colData$reducedDims$UMAP
    int.embed <- Embeddings(obj, reduction = "umap")
    int.embed <- int.embed[rownames(cds.embed),]
    cds@int_colData$reducedDims$UMAP <- int.embed
  }
  pdf('int.umap.pdf',width = 6,height = 5)
  p1 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = group.by, show_trajectory_graph = F) + ggtitle('int.umap')+
    scale_color_manual(values = color)
  print(p1)
  dev.off()
  
  cds <- cluster_cells(cds, reduction_method='UMAP')
  print(length(unique(partitions(cds))))
  
  pdf('partition.umap.pdf',width = 6,height = 5)
  p1 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = F, group_label_size = group_label_size)
  print(p1)
  dev.off()
  
  cds <- learn_graph(cds)
  
  pdf('trajectory.umap.pdf',width = 6,height = 4)
  p1 <- plot_cells(cds, color_cells_by = group.by, show_trajectory_graph = T, label_cell_groups = label_cell_groups, 
                   label_groups_by_cluster = F, label_leaves = label_leaves,
                   label_branch_points = label_branch_points, group_label_size = group_label_size)+
    scale_color_manual(values = color)
  print(p1)
  dev.off()
  save(cds, file = "cds.RData")
  return(cds)
}



