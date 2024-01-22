# 7. Multimodal intersection analysis: MIA
#' @param spe seurat object of ST
#' @param sc.markers.list Signatures of CCSs
#' @param st.markers.list Signatures of ST pathological tissues

MIA <- function(spe, sc.markers.list, st.markers.list){
  
  library(Seurat)
  library(stats)
  
  M <- length(sc.markers.list)
  N <- length(st.markers.list) 
  MIA.results <- matrix(0,nrow = M, ncol = N)
  rownames(MIA.results) <- names(sc.markers.list)
  colnames(MIA.results) <- names(st.markers.list)
  # Gene universe
  gene.universe <- length(rownames(spe))
  
  # Loop over ST clusters
  for (i in 1:N) {
    # Then loop over SC clusters
    for (j in 1:M) {
      genes1 <- st.markers.list[[i]]
      genes2 <- sc.markers.list[[j]]
      
      # Hypergeometric    
      A <- length(intersect(genes1,genes2))
      B <- length(genes1)
      C <- length(genes2)
      enr <- -log10(phyper(A, B, gene.universe-B, C, lower.tail = FALSE))
      dep <- -log10(1-phyper(A, B, gene.universe-B, C, lower.tail = FALSE))
      if (enr < dep) {
        MIA.results[j,i] = -dep
      } else {
        MIA.results[j,i] = enr
      } } }
  
  MIA.results <- t(MIA.results)
  # Some results were -Inf...check why this is the case...
  MIA.results[is.infinite(MIA.results)] <- 0
  return(MIA.results)
} 



