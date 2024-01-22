# scRNA-seq 

# raw.mat: raw count matrix

# Empty drops 

library(Seurat)
library(DropletUtils)
set.seed(100)
out <- emptyDrops(raw.mat)

# Low quality cells

raw.seurat <- CreateSeuratObject(counts = raw.mat, project = "brca", min.features = 200, min.cells =5, names.field = 1,names.delim = "-")
raw.seurat[["percent.mt"]] <- PercentageFeatureSet(raw.seurat, pattern = "^MT-")
VlnPlot(raw.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)
raw.seurat <- subset(raw.seurat, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 20)

# Cell clustering

raw.seurat <- NormalizeData(raw.seurat, normalization.method =  "LogNormalize", scale.factor = 10000)
GetAssay(raw.seurat,assay = "RNA")
raw.seurat <- FindVariableFeatures(raw.seurat, selection.method = "vst", nfeatures = 2000) 
raw.seurat <- ScaleData(raw.seurat) 
raw.seurat <- RunPCA(object = raw.seurat, pc.genes = VariableFeatures(raw.seurat)) 

# Choosing the best pc

ElbowPlot(raw.seurat, ndims=30, reduction="pca")
pct <- raw.seurat[["pca"]]@stdev / sum(raw.seurat[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)

# Cell clustering

raw.seurat <- FindNeighbors(raw.seurat, dims =1:pcs)
raw.seurat <- FindClusters(raw.seurat, resolution = 0.4)
raw.seurat <- RunUMAP(raw.seurat, dims = 1:pcs)

# Batch effect correction harmony

library(harmony)
raw.seurat <- RunHarmony(raw.seurat, group.by.vars = "patients", plot_convergence = T)

# Cell clustering

raw.seurat <- FindNeighbors(raw.seurat, dims =1:pcs, reduction = "harmony")
raw.seurat <- FindClusters(raw.seurat, resolution = 0.4)
raw.seurat <- RunUMAP(raw.seurat, dims = 1:pcs, reduction = "harmony")

# Doublet detection and removal
# filtering in combination with cell type markers
#' @param object seurat object
#' @param group cell groups

run_fDblClusters <- function(object, group = "seurat_clusters"){
  library(scDblFinder)
  db.test <- findDoubletClusters(GetAssayData(object, slot="counts", assay="RNA"),
                                 clusters = object@meta.data[[group]])
  library(scater)
  chosen.doublet <- rownames(db.test)[isOutlier(db.test$num.de, type="lower", log=TRUE)]
  return(chosen.doublet)
}


