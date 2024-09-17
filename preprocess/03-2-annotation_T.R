setwd("D:/Projects/CASH-Human/")
# 1. Load packages
library(Seurat)
library(ggplot2)
library(dplyr)
source("./bin/function.R")


# 2. params
dir_for_T_obj <- "./data/tmp/01.preprocess/sce_subset/T_obj_sub2.rds"
## load data
T_obj <- readRDS(dir_for_T_obj)

# 3. markers
markers_list <- list(
  mixture = c("ALB", "C1QA", "S100A9", "S100A8", "CALCRL", "FCN1"), 
  clone = c("KLRF1", "FCER1G"), 
  T_cells = c("CD3D", "CD4", "CD8A"), 
  gdT = c("TRDV2", "TRGV9"), 
  Naive = c("SELL", "CCR7", "LEF1"), 
  Cytotoxicity = c("GZMB"), 
  Exhaustion = c("LAG3", "TIGIT", "PDCD1", "TOX"), 
  Effect = c("CX3CR1", "GNLY"), 
  EffectMemory = c("GZMK", "CXCR3"), 
  MAIT = c("SLC4A10"), 
  TT = c(
    #"CD3", 
    "NCAM1"),
  fat = c("MS4A1", "CPA3"
          #, "TPSAB2"
          )
)
markers <- unique(unlist(markers_list))


# 4. UMAP_Pct_Dotplot
pdf("./fig/tmp/01.preprocess/T_UMAP_Pct_Dotplot.pdf", width = 14, height = 14)
resolution_list <- paste0("RNA_snn_res.", seq(0.1, 1, 0.1))
for (res in resolution_list) {
  # res <- "RNA_snn_res.0.7"
  ## UMAP
  p1 <- DimPlot(T_obj, label = T, group.by = res, label.size = 6) + NoLegend() + ggtitle(paste0("UMAP-T_", res))
  
  ## Dotplot
  p2 <-  yy_Dotplot(seuratObj = T_obj,
                    genes = markers_list,
                    group.by = res) + ggtitle(paste0("Dotplot-T_", res))
  
  ## Pct bar plot
  meta <- T_obj@meta.data
  tmp <- meta[, c(res, "orig.ident", "group")]
  colnames(tmp) <- c("cluster", "orig.ident", "group")
  
  num_cluster <- tmp %>%
    group_by(cluster, orig.ident) %>%
    summarise(n_clsuter = n(), .groups = 'drop')
  
  num_total <- tmp %>%
    group_by(orig.ident) %>%
    summarise(n_total = n(), .groups = 'drop')
  
  plot_df <- num_cluster %>%
    left_join(num_total, by = "orig.ident") %>%
    mutate(pct = (n_clsuter / n_total) * 100)
  plot_df$group <- stringr::str_split(plot_df$orig.ident, "_", simplify = T)[, 1]
  
  p3 <- ggplot(plot_df, aes(x = group, y = pct)) +
    geom_boxplot(aes(color = group), outlier.shape = NA, lwd = 1) +
    geom_jitter(aes(color = group), size = 1, width = 0.2) +
    facet_wrap(~ cluster, nrow = 1) +
    ggsignif::geom_signif(
      comparisons = list(c("CASH", "Normal")),
      test = wilcox.test,
      map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
      step_increase = -0.1,
      textsize = 4,
      size = 1
    )+
    cowplot::theme_cowplot() +
    scale_color_manual(values = c("CASH" = "#5bc0eb", "Normal" = "#9bc53d")) +
    labs(x = "group", y = "pct of cluster(%)", title = "proportion of each subpopulation cell in different samples", fill = "group") +
    theme(axis.text = element_text(size = 12, angle = 0, hjust = 1),
          axis.title = element_text(size = 14),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 1),
          legend.position =  "bottom",
          legend.justification = "center", 
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 12),
          plot.title = element_text(size = 14, hjust = 0.5))
  combined_plot <- (p1+p3) / p2
  print(combined_plot)
}
dev.off()

# 5. Annotation
Idents(T_obj) <- "RNA_snn_res.0.7"
T_obj <- RenameIdents(T_obj,
                      "0" = "c14_CD8_MAIT",
                      "1" = "c17_CD8_",
                      "2" = "c18_CD8_",
                      "3" = "c15_CD8_MAIT",
                      "4" = "c25_CD4_",
                      "5" = "c20_CD8_",
                      "6" = "c21_CD8_",
                      "7" = "c22_CD8_",
                      "8" = "c24_CD4_",
                      "9" = "c16_CD8_MAIT",
                      "10" = "c19_CD8_",
                      "11" = "c27_gdT",
                      "12" = "c23_CD4_Naïve",
                      "13" = "c26_CD4_Treg")
T_obj$sub_celltype <- Idents(T_obj)
levels(T_obj) <- unique(sort(levels(T_obj)))
T_obj$sub_celltype <- factor(T_obj$sub_celltype, levels = unique(sort(levels(T_obj))))

## check
DimPlot(T_obj, group.by = c("RNA_snn_res.0.7", "sub_celltype"), label = T) + NoLegend()
saveRDS(T_obj, "./data/obj/T_obj.rds")
