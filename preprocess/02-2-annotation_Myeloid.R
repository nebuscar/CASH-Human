setwd("D:/Projects/CASH-Human/")
# 1. Load packages
library(Seurat)
library(ggplot2)
library(dplyr)
source("./bin/function.R")


# 2. params
dir_for_Myeloid_obj <- "./data/tmp/01.preprocess/sce_subset/Myeloid_obj_sub2.rds"
## load data
Myeloid_obj <- readRDS(dir_for_Myeloid_obj)


# 3. markers
markers_list <- list(
  c("CD68", "C1QA", "CD163", "SPP1", "LYVE1", "FOLR2", "MARCO"),
  c("FCN1", "VCAN", "CD14", "FCGR3A", "S100A9", "S100A8", "MMP19"),
  c("CLEC4F"),
  c("FLT3", 
    #"CD123", 
    "CD1E", "LAMP3", "IDO1", "IDO2"),
  c("CLEC9A"),
  c("CD1C"),
  c("CCR7"),
  c("LILRA4", "CLEC4C", "JCHAIN"),
  c("MS4A1", "CPA3"
    #, "TPSAB2"
    ),
  c("NKG7", "CD3D"),
  c("MKI67")
)
markers <- unique(unlist(markers_list))


########## 4. UMAP_Pct_Dotplot ########## 
pdf("./fig/tmp/01.preprocess/Myeloid_UMAP_Pct_Dotplot.pdf", width = 14, height = 14)
resolution_list <- paste0("RNA_snn_res.", seq(0.1, 1, 0.1))
for (res in resolution_list) {
  ## UMAP
  p1 <- DimPlot(Myeloid_obj, label = T, group.by = res, label.size = 6) + NoLegend() + ggtitle(paste0("UMAP-Myeloid_", res))
  
  ## Dotplot
  p2 <-  yy_Dotplot(seuratObj = Myeloid_obj,
                    genes = markers_list,
                    group.by = res) + ggtitle(paste0("Dotplot-Myeloid_", res))
  
  ## Pct bar plot
  meta <- Myeloid_obj@meta.data
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

########## 5. Annotation ########## 
Idents(Myeloid_obj) <- "RNA_snn_res.0.7"
Myeloid_obj <- RenameIdents(Myeloid_obj,
                            "0" = "c01_Monocyte_CD14CD16",
                            "1" = "c09_cDC2",
                            "2" = "c02_Monocyte_CD14",
                            "3" = "c06_Macrophage_",
                            "4" = "c05_Macrophage_",
                            "5" = "c08_cDC1",
                            "6" = "c07_Macrophage_MARCO",
                            "7" = "c03_Monocyte_CD16",
                            "8" = "c10_cDC2",
                            "9" = "c04_Monocyte_CD14",
                            "10" = "c11_cDC2",
                            "11" = "c12_pDC",
                            "12" = "c13_Actived_DC")
Myeloid_obj$sub_celltype <- Idents(Myeloid_obj)
levels(Myeloid_obj) <- unique(sort(levels(Myeloid_obj)))
Myeloid_obj$sub_celltype <- factor(Myeloid_obj$sub_celltype, levels = unique(sort(levels(Myeloid_obj))))

## check
DimPlot(Myeloid_obj, group.by = "sub_celltype", label = T) + NoLegend()
saveRDS(Myeloid_obj, "./data/obj/Myeloid_obj.rds")
