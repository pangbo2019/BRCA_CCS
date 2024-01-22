#' Cell-cell interaction 
#' @param all.seurat seurat object with cell type lables

library(Seurat)
library(CellChat)

# Create cellchat object
cellchat <- createCellChat(all.seurat, meta = all.seurat@meta.data, group.by = "celltype")
# Set the reference database
CellChatDB.use <- subsetDB(CellChatDB = CellChatDB.human, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use 

# Ligand-receptor analysis
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

# Calculate the communication probability
# Infer the cellchat network
cellchat <- computeCommunProb(cellchat, raw.use = F,population.size = T)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Return a dataframe
df.net <- subsetCommunication(cellchat)


