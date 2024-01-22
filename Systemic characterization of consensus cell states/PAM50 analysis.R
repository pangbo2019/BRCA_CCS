# 5. PAM50 analysis
#' @param PAM50.logtpm logTPM matrix, rows cells, columns 50 subtype-genes

library(genefu)
library(org.Hs.eg.db)
s2g <- toTable(org.Hs.egSYMBOL)
g <- s2g[match(colnames(PAM50.logtpm), s2g$symbol),1];

dannot <- data.frame(probe = colnames(PAM50.logtpm),
                     "Gene.Symbol" = colnames(PAM50.logtpm), 
                     "EntrezGene.ID" = g)
data(pam50.robust)
data(pam50)
# run genefu
S.tpm <- molecular.subtyping(sbt.model ="pam50", data=PAM50.logtpm, annot=dannot, do.mapping=F)


