setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(dplyr)
library(stringr)
library(AUCell)
library(ggplot2)
library(ggpubr)


# 2.params
dir_for_T_obj <- "./data/obj/T_obj.rds"
dir_for_meta_tcr <- "./data/tcr/meta_tcr.rds"

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


# 3. geneSets
geneSets <- list(
  exhaustion = c("CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT", "TIM3", "TOX"),
  effector = c("PRF1", "IFNG", "GNLY", "NKG7", "GZMB", "GZMA", "CST7", "TNFSF10"),
  resident = c("RUNX3", "NR4A1", "CD69", "CXCR6", "NR4A3"))
markers <- unique(unlist(geneSets))

# 4. AUCell
cells_rankings <- AUCell_buildRankings(CD8T_obj@assays$RNA$data) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
# cells_AUC2 <- AUCell_run(CD8T_obj@assays$RNA$data, geneSets)
# score <- as.data.frame(t(getAUC(cells_AUC)))
CD8T_obj@assays$AUCell <- cells_AUC@assays@data$AUC


##########
scScore <- data.frame(t(CD8T_obj[["AUCell"]]), 
                      cloneGroup = CD8T_obj$cloneGroup)
scScore <- reshape2::melt(scScore, id.vars = "cloneGroup", variable.name = "AUCell", value.name = "Score")
scScore %>%
  ggplot(aes(x = cloneGroup, y = as.numeric(Score), color = cloneGroup)) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.1)+
  facet_wrap(~ AUCell, scales = "free_y") +
  stat_compare_means(aes(group = cloneGroup), method = "wilcox.test", label = "p.signif", label.y = 1) +
  theme(axis.text.x = element_text(angle = 0, hjust=1, size = 8),
        axis.text.y = element_text(angle = 0, hjust=1, size = 8),
        # axis.title = element_blank(),
        legend.position="right",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


## calculate mean value
mean_scores <- scScore %>%
  group_by(cloneGroup, AUCell) %>%
  summarize(mean_score = mean(as.numeric(Score)))
## calculate p value
wilcox_tests <- scScore %>%
  group_by(AUCell) %>%
  summarize(p_value = wilcox.test(Score ~ cloneGroup)$p.value)
plots <- list()
for (aucell in unique(scScore$AUCell)) {
  data_subset <- scScore %>% filter(AUCell == aucell)
  mean_subset <- mean_scores %>% filter(AUCell == aucell)
  p_value <- wilcox_tests %>% filter(AUCell == aucell) %>% pull(p_value)
  
  p <- ggplot(data_subset, aes(x = cloneGroup, y = as.numeric(Score), color = cloneGroup)) + 
    geom_violin(trim = F) +
    geom_boxplot(width = 0.1) +
    geom_text(data = mean_subset, aes(x = cloneGroup, y = mean_score, label = sprintf("%.4f", mean_score)),
              vjust = -0.5, size = 4, color = "black") +
    annotate("text", x = 1.5, y = max(data_subset$Score) + 0.05, 
             label = sprintf("Wilcoxon, p=%.4g", p_value), size = 4) +
    labs(title = aucell, x = "Clone Group", y = "AUCellScore") +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 1, size = 10),
      # axis.text.x = element_blank(),
      axis.text.y = element_text(angle = 0, hjust = 1, size = 10),
      axis.title = element_blank(),
      axis.title.y = element_text(size = 14),
      legend.position = "none",
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14)) 
  
  plots[[aucell]] <- p
}
combined_plot <- patchwork::wrap_plots(plots, ncol = 3)
print(combined_plot)
ggsave("./Vlnplot_signature_expr_CD8T_cloneGroup.pdf", combined_plot, width = 12, height = 6)



