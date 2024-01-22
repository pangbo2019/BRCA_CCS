# 3. copykat 
#' @param patient.seurat subclasses of processed.seurat

library(copykat)
copykat <- copykat(rawmat=patient.seurat@assays$RNA@counts, 
                   id.type="S",   
                   ngene.chr=5, 
                   win.size=25, 
                   KS.cut=0.1, 
                   sam.name="test", 
                   distance="euclidean", 
                   n.cores=1,
                   genome = "hg20",
                   plot.genes="TRUE")

pred.test <- data.frame(copykat$prediction)


