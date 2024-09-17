setwd("D:/Projects/CASH-Human/")
# 1. Load packages
library(Seurat)
library(ggplot2)
library(dplyr)
source("./bin/function.R")


# 2. params
dir_for_B_obj <- "./data/tmp/01.preprocess/sce_subset/B_obj_sub.rds"
## load data
B_obj <- readRDS(dir_for_B_obj)


# 3. markers
markers_list = list(
  c("FCER2", "TCL1A", "IL4R", "CD72", "BACH2", "IGHD", "IGHM"),
  c("NR4A1", 'NR4A2', "CREM", "CD83"),
  c("ISG15", "IFI44L", "IFI6", "IFIT3"),
  
  c("CD27", "TNFRSF13B","TXNIP", "GPR183"),
  c("HSPA1A", "HSPA1B","DNAJB1", "EGR1"),
  c("FCRL4","FCRL5", "ITGAX", "TBX21", "CR2"),
  c("CCR1", "CXCR3", "PDCD1", "HCK", "FCRL3", "FGR"),
  c("NME1", "APEX1", "POLD2", "POLE3","MYC"),
  
  c("BCL6", "RGS13", "AICDA","IL21R"),
  c("MKI67","STMN1","HMGB2","TOP2A"),
  
  c("CD38","MZB1","PRDM1","IRF4","XBP1"),
  c("MS4A1","LTB","HLA-DRA", "HLA-DRB1","HLA-DPA1","HLA-DQA1"),
  c("IGHG1", "IGHG2", "IGHG3", "IGHG4","IGHA1", "IGHA2"),
  c("IL10","IL12A","EBI3","TGFB1")
  )
markers <- unique(unlist(markers_list))

# 4. UMAP_Pct_Dotplot
pdf("./fig/tmp/01.preprocess/B_UMAP_Pct_Dotplot.pdf", width = 14, height = 14)
resolution_list <- paste0("RNA_snn_res.", seq(0.1, 1, 0.1))
for (res in resolution_list) {
  # res <- "RNA_snn_res.0.3"
  ## UMAP
  p1 <- DimPlot(B_obj, label = T, group.by = res, label.size = 6) + NoLegend() + ggtitle(paste0("UMAP-B_", res))
  
  ## Dotplot
  p2 <-  yy_Dotplot(seuratObj = B_obj,
                  genes = markers_list,
                  group.by = res) + ggtitle(paste0("Dotplot-B_", res))
  
  ## Pct bar plot
  meta <- B_obj@meta.data
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
Idents(B_obj) <- "RNA_snn_res.0.3"
B_obj <- RenameIdents(B_obj, 
                      "0" = "c35_PC",
                      "1" = "c36_Bn_TCL1A",
                      "2" = "c37_classical-Bm_GRP183",
                      "3" = "c38_classical-Bm_TXNIP",
                      "4" = "c39_ABC_FGR",
                      "5" = "c40_Bgc_MKI67"
                      )
B_obj$sub_celltype <- Idents(B_obj)
levels(B_obj) <- unique(sort(levels(B_obj)))
B_obj$sub_celltype <- factor(B_obj$sub_celltype, levels = unique(sort(levels(B_obj))))

## check
DimPlot(B_obj, group.by = c("RNA_snn_res.0.7", "sub_celltype"), label = T) + NoLegend()
saveRDS(B_obj, "./data/obj/B_obj.rds")
