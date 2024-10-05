setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(irGSEA)
library(AUCell)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)


# 2. params
dir_for_T_obj <- "./data/obj/T_obj.rds"
dir_for_meta_tcr <- "./data/tcr/meta_tcr.rds"
Subset_color_panel <- c(
  # T cells
  "c14_CD8_MAIT" = "#5A7B8F",
  "c15_CD8_MAIT" = "#EE4C97",
  "c16_CD8_MAIT" = "#3D806F",
  "c17_CD8_" = "#A08634",
  "c18_CD8_" = "#F37C95",
  "c19_CD8_" = "#608541",
  "c20_CD8_" = "#7D4E57",
  "c21_CD8_" = "#BC3C29",
  "c22_CD8_" = "#958056",
  "c23_CD4_Naïve" = "#9FAFA3",
  "c24_CD4_" = "#6F99AD",
  "c25_CD4_" = "#0072B5",
  "c26_CD4_Treg" = "#CFC59A",
  "c27_gdT" = "#E18727")

## load data
T_obj <- readRDS(dir_for_T_obj)
meta_tcr <- readRDS(dir_for_meta_tcr)

## match meta
meta_match <- meta_tcr[match(rownames(T_obj@meta.data), rownames(meta_tcr)),]
T_obj@meta.data$cloneSize <- meta_match$cloneSize
## add cloneGroup meta
T_obj@meta.data <- T_obj@meta.data %>%
  mutate(cloneGroup = case_when(
    cloneSize %in% c("Hyperexpanded (0.1 < X <= 1)", "Large (0.01 < X <= 0.1)", "Medium (0.001 < X <= 0.01)") ~ "highClone",
    cloneSize %in% c("Small (1e-04 < X <= 0.001)", "Rare (0 < X <= 1e-04)", "None ( < X <= 0)") ~ "lowClone"
  ))
T_obj <- subset(T_obj, cloneGroup %in% c("highClone", "lowClone"))
## subset CD8T
meta_sub <- T_obj@meta.data %>%
  filter((str_detect(sub_celltype, "CD8")))
CD8T_obj <- T_obj[, colnames(T_obj) %in% rownames(meta_sub)]
rm(T_obj, meta_tcr, meta_match, meta_sub)

########## 3. AUCell-Signature_score ##########
## geneSets
geneSets <- list(
  exhaustion = c("CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT", "TIM3", "TOX"),
  effector = c("PRF1", "IFNG", "GNLY", "NKG7", "GZMB", "GZMA", "CST7", "TNFSF10"),
  resident = c("RUNX3", "NR4A1", "CD69", "CXCR6", "NR4A3"))

markers <- unique(unlist(geneSets))
geneSets[["exhaustion"]] <- intersect(geneSets[["exhaustion"]], row.names(CD8T_obj))
geneSets[["effector"]] <- intersect(geneSets[["effector"]], row.names(CD8T_obj))
geneSets[["resident"]] <- intersect(geneSets[["resident"]], row.names(CD8T_obj))
length(geneSets[["exhaustion"]])
length(geneSets[["effector"]])
length(geneSets[["resident"]])

##AUC_score
cells_rankings <- AUCell_buildRankings(CD8T_obj@assays$RNA$data) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
# cells_AUC2 <- AUCell_run(CD8T_obj@assays$RNA$data, geneSets)
# score <- as.data.frame(t(getAUC(cells_AUC)))
CD8T_obj@assays$AUCell <- cells_AUC@assays@data$AUC

scScore <- data.frame(t(CD8T_obj[["AUCell"]]), 
                      cloneGroup = CD8T_obj$cloneGroup)
scScore <- reshape2::melt(scScore, id.vars = "cloneGroup", variable.name = "AUCell", value.name = "Score")

plot_df <- scScore
plot_df$cloneGroup <- factor(plot_df$cloneGroup, levels = c("highClone", "lowClone"))
table(plot_df$cloneGroup)

## plot
p <- ggplot(plot_df, aes(x = cloneGroup, y = Score)) + 
  geom_boxplot(aes(color = cloneGroup), outlier.shape = NA, lwd = 1) +
  ## calculate sig
  ggsignif::geom_signif(
    comparisons = list(
      c("highClone", "lowClone")),
    test = function(x, y) wilcox.test(x, y, alternative = "two.sided"),
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
    step_increase = 0.1,
    textsize = 4,
    size = 1) +
  facet_wrap(~ AUCell, nrow = 1) +
  scale_color_manual(values = c("highClone" = "#F37C95",
                                "lowClone" = "#A08634")) +
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
ggsave("./fig/Fig04.T/Fig4F.AUCell_signature_score_T.pdf", plot = p,  width = 8, height = 8, dpi = 300)
ggsave("./fig/Fig04.T/Fig4F.AUCell_signature_score_T.png", plot = p,  width = 8, height = 8, dpi = 300)

