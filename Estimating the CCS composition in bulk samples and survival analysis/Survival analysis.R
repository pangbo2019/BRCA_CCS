# Survival analysis
# survival_dat: a dataframe with time, status and PATIENT_ID
# CIBERSORT.Res: a matrix, rows are CCSs, colums are samples

p <- pheatmap(CIBERSORT.Res, color= colorRampPalette(c("Darkblue", "white","red"))(100),
              border_color = "grey", fontsize = 10, fontsize_row = 10, fontsize_col = 10,
              show_rownames = T, show_colnames = F)

col_cluster <- cutree(p$tree_col, k = 4)
all(survival_dat$PATIENT_ID == names(col_cluster))
survival_dat$group <- col_cluster

library(survival)
coxph(Surv(time, status) ~ group,data = survival_dat)

library(survminer)
f1 <- survfit(Surv(time, status)~ group, data = survival_dat)
pdf("cibersort.survival.pdf",width = 8,height = 8)
ggsurvplot(f1, 
           data = survival_dat,
           surv.median.line = "hv",
           legend.title = "Group",
           legend.labs = c("1", "2","3","4"),
           pval = TRUE,            
           conf.int = F,
           risk.table = T)
dev.off()


