# Deconvolution analyses of spatial data
#' @param sc.object seurat object of scRNA-seq with annotation of patient labels 
#' @param spe seurat object of ST
#' @param CCS CCS classification labels of cells
#' @param slicename name of ST slice
#' @param resolution original resolution or high resolution 
#' @param NumGridsThe initial number of newly grided spatial locations.

ST_CARD_deconvolution_calculation <- function(sc.object, spe, CCS, slicename, 
                                              resolution=c('original','refined'), 
                                              NumGrids = 2000){
  library(devtools)
  library(CARD)
  library(SingleCellExperiment)
  library(TOAST)
  library(MuSiC)
  
  sc_count <- GetAssayData(sc.object, assay="RNA", slot='counts')
  
  cellID <- colnames(sc.object)
  CCS <- CCS
  sampleInfo <- sc.object@meta.data$patient
  sc_meta <- data.frame(cellID = cellID,CCS = CCS,sampleInfo = sampleInfo)
  rownames(sc_meta) <- cellID
  
  x <- spe@images[[slicename]]@coordinates[["row"]]
  y <- spe@images[[slicename]]@coordinates[["col"]]
  spatial_location <- data.frame(x = x,y = y)
  rownames(spatial_location) <- colnames(spe@assays$SCT$counts) 
  
  spatial_count <- spe@assays$SCT$counts
  
  # 1.Create an CARD object
  CARD_obj = createCARDObject(
    sc_count = sc_count,
    sc_meta = sc_meta,
    spatial_count = spatial_count,
    spatial_location = spatial_location,
    ct.varname = "CCS",
    ct.select = unique(sc_meta$CCS),
    sample.varname = "sampleInfo",
    minCountGene = 100,
    minCountSpot = 5) 
  # 2.Deconvolution using CARD
  CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)  
  if(resolution=='refined'){
    # Refined spatial map
    # Imputation on the newly grided spatial locations
    CARD_obj = CARD.imputation(CARD_obj,NumGrids,ineibor = 10,exclude = NULL)  
  }
  return(CARD_obj)
}


