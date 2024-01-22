# Go enrichment
# gene: each differentially expressed geneset of hcs

library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

data(geneList, package = "DOSE") 
GO <- bitr(gene,
           fromType = "SYMBOL", 
           toType = c("ENSEMBL", "ENTREZID"), 
           OrgDb = org.Hs.eg.db)

enrich <- enrichGO(gene = GO$ENTREZID, 
                   universe = names(geneList), 
                   OrgDb = org.Hs.eg.db, 
                   keyType = "ENTREZID",
                   ont = "BP",
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)



