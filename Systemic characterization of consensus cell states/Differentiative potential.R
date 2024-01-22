# 4. Evaluating cell stemness

#' CytoTRACE
#' @param maligant.seurat seurat object of malignant cell

library(CytoTRACE)
CytoTRACE_results <- CytoTRACE(as.matrix(maligant.seurat@assays$RNA@counts))

#' SCENT
#' @param sc_exp count matrix, rows Entrez_id, columns cells

library(SCENT)
run_SCENT <- function(sc_exp, ver = 2){
  lscChu0.m <- log2(sc_exp+1)
  if(ver==1){
    integ.l <- SCENT::DoIntegPPI(exp.m = lscChu0.m, ppiA.m = SCENT::net13Jun12.m)
    SR <- SCENT::CompSRana(integ.l, local = FALSE, mc.cores = 1)
    return(SR)
  }
  else if(ver==2){
    SR <- SCENT::CompCCAT(exp = lscChu0.m, ppiA = SCENT::net13Jun12.m)
    return(SR)
  }
  else{
    print("wrong version parameter")
  }
}



