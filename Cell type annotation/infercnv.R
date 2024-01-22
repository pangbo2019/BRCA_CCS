# 2. infercnv
#' @param patient.seurat subclasses of processed.seurat

library(infercnv)
# create normal reference  
case_ctrl <- data.frame(cell = colnames(patient.seurat), celltype = patient.seurat$celltype)
case_ctrl$celltype <- ifelse(case_ctrl$celltype %in% "Epithelial_cells","malignant","non-malignant")
write.table(case_ctrl, "cellAnnotations.txt", quote = F, sep = "\t", col.names=F, row.names=F)

# run InferCNV
infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = patient.seurat@assays$RNA@counts, 
                                     annotations_file = "cellAnnotations.txt",
                                     delim = "\t",
                                     gene_order_file = "geneOrderingFile.txt",
                                     ref_group_names = "non-malignant")
infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1, 
                              out_dir = "output_dir", 
                              cluster_by_groups = F,
                              denoise = T, 
                              HMM = F,
                              analysis_mode = "samples",
                              tumor_subcluster_partition_method = "leiden",
                              hclust_method = "ward.D2",
                              k_obs_groups = 1)

# generate cluster metric plots for epithelial cluster
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

expr <- infercnv_obj@expr.data
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <- normal_loc$`non-malignant`
test_loc <- infercnv_obj@observation_grouped_cell_indices
test_loc <- test_loc$malignant

anno.df <- data.frame(
  CB = c(colnames(expr)[normal_loc],colnames(expr)[test_loc]),
  class = c(rep("normal",length(normal_loc)),rep("test",length(test_loc)))
)
head(anno.df)

gn <- rownames(expr)
geneFile <- read.table("geneOrderingFile.txt",header = F,sep = "\t",stringsAsFactors = F)
rownames(geneFile) <-  geneFile$V1
sub_geneFile <-  geneFile[intersect(gn,geneFile$V1),]
expr <- expr[intersect(gn,geneFile$V1),]
sub_geneFile$V2 <- paste("chr",sub_geneFile$V2,sep = "")

# cluster cells
set.seed(123)
kmeans.result <- kmeans(t(expr), 6)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$CB <- rownames(kmeans_df)
kmeans_df <- kmeans_df%>%inner_join(anno.df,by="CB") 
kmeans_df_s <- arrange(kmeans_df,kmeans_class) 
rownames(kmeans_df_s) <- kmeans_df_s$CB
kmeans_df_s$CB <- NULL
kmeans_df_s$kmeans_class <- as.factor(kmeans_df_s$kmeans_class) 

# Define heatmap annotations and color matching
top_anno <- HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "NA",col="NA"), labels = 1:22,labels_gp = gpar(cex = 1.5)))
color_v <- RColorBrewer::brewer.pal(8, "Dark2")[1:6] 
names(color_v) <- as.character(1:6)
left_anno <- rowAnnotation(df = kmeans_df_s,col=list(class=c("test"="red","normal" = "blue"),kmeans_class=color_v))

# plot
pdf("infercnv.pdf",width = 15,height = 10)
Heatmap(t(expr)[rownames(kmeans_df_s),], 
        col = colorRamp2(c(0.85,1,1.15), c("#377EB8","#F0F0F0","#E41A1C")), 
        cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = F,
        column_split = factor(sub_geneFile$V2, paste("chr",1:22,sep = "")), 
        column_gap = unit(2, "mm"),
        
        heatmap_legend_param = list(title = "Modified expression",direction = "vertical",title_position = "leftcenter-rot",at=c(0.4,1,1.6),legend_height = unit(3, "cm")),
        
        top_annotation = top_anno,left_annotation = left_anno,
        row_title = NULL,column_title = NULL)
dev.off()


