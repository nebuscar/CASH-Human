setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(irGSEA)
library(AUCell)
library(dplyr)
library(ggplot2)


# 2. params
dir_for_Myeloid_obj <- "./data/obj/Myeloid_obj.rds"
Subset_color_panel <- c(
  # Myeloid cells
  "c01_Monocyte_CD14CD16" = "#5A7B8F",
  "c02_Monocyte_CD14" = "#3D806F",
  "c03_Monocyte_CD16" = "#BC3C29",
  "c04_Monocyte_CD14" = "#9FAFA3",
  "c05_Macrophage_" = "#F37C95",
  "c06_Macrophage_" = "#A08634",
  "c07_Macrophage_MARCO" = "#7D4E57",
  "c08_cDC1" = "#608541",
  "c09_cDC2" = "#EE4C97",
  "c10_cDC2" = "#958056",
  "c11_cDC2" = "#6F99AD",
  "c12_pDC" = "#0072B5",
  "c13_Actived_DC" = "#CFC59A")
## load data
Myeloid_obj <- readRDS(dir_for_Myeloid_obj)
meta <- Myeloid_obj@meta.data

########## 3. AUCell-Signature_score ##########
## geneSets
geneSets <- list(
  # M1
  M1_signature=c("CCL5", "CCR7", "CD40", "CD86", "CXCL9", 
                 "CXCL10", "CXCL11", "IDO1", "IL1A", "IL1B", 
                 "IL6", "IRF1", "IRF5", "KYNU"),
  # M2
  M2_signature=c("CCL4", "CCL13", "CCL18", "CCL20", "CCL22", 
                 "CD276", "CLEC7A", "CTSA", "CTSB", "CTSC", 
                 "CTSD", "FN1", "IL4R", "IRF4","LYVE1", 
                 "MSR1", "TGFB1", "TGFB2", "TGFB3", "TNFSF8", 
                 "TNFSF12", "VEGFA", "VEGFB", "VEGFC")
  )
geneSets[["M1_signature"]] <- intersect(geneSets[["M1_signature"]], row.names(Myeloid_obj))
geneSets[["M2_signature"]] <- intersect(geneSets[["M2_signature"]], row.names(Myeloid_obj))
length(geneSets[["M1_signature"]])
length(geneSets[["M2_signature"]])

##AUC_score
exprMatrix <- Myeloid_obj@assays$RNA$counts
exprMatrix <- as(exprMatrix, "dgCMatrix")

cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
aucs <- getAUC(cells_AUC)
aucs_t <- data.frame(t(aucs))

Myeloid_obj@meta.data$m1_signature <- aucs_t[,1]
Myeloid_obj@meta.data$m2_signature <- aucs_t[,2]

tmp <- Myeloid_obj@meta.data[, c("sub_celltype", "m1_signature", "m2_signature")]
colnames(tmp) <- c("cluster", "M1 phenotype", "M2 phenotype")
tmp <- reshape2::melt(tmp, id.var = "cluster", variable.name = "pathway", value.name = "score")
plot_df <- tmp[tmp$cluster %in% c("c05_Macrophage_", "c06_Macrophage_", "c07_Macrophage_MARCO"),]
plot_df$cluster <- factor(plot_df$cluster, levels = c("c06_Macrophage_", "c05_Macrophage_", "c07_Macrophage_MARCO"))
plot_df$cluster_short <- stringr::str_split(plot_df$cluster, "_", simplify = T)[,1]
plot_df$cluster_short <- factor(plot_df$cluster_short, levels = c("c06", "c05", "c07"))
table(plot_df$cluster_short)

## plot
p <- ggplot(plot_df, aes(x = cluster_short, y = score)) + 
  geom_boxplot(aes(color = cluster_short), outlier.shape = NA, lwd = 1) +
  ## calculate sig
  ggsignif::geom_signif(
    comparisons = list(
      c("c05", "c06"), 
      c("c05", "c07"),
      c("c06", "c07")),
    test = function(x, y) wilcox.test(x, y, alternative = "two.sided"),
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
    step_increase = 0.1,
    textsize = 4,
    size = 1) +
  facet_wrap(~ pathway, nrow = 1) +
  scale_color_manual(values = c("c05" = "#F37C95",
                                "c06" = "#A08634",
                                "c07" = "#7D4E57")) +
  cowplot::theme_cowplot() +
  labs(x = "", y = "Signature expression", title = "") +
  theme(
    ## axis
    axis.text.x = element_text(size = 14, angle = 0),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.line = element_line(linewidth = 1),
    
    ## facet
    strip.text.x = element_text(size = 14, face = "bold"),
    strip.background = element_rect(color = NA, fill = NA, linewidth = 1),
    
    ## legend
    legend.justification = "center",
    legend.position = "none"
  )
print(p)
ggsave("./fig/Fig02.Myeloid/Fig2D.AUCell_signature_score_Myeloid.pdf", plot = p,  width = 8, height = 8, dpi = 300)
ggsave("./fig/Fig02.Myeloid/Fig2D.AUCell_signature_score_Myeloid.png", plot = p,  width = 8, height = 8, dpi = 300)

