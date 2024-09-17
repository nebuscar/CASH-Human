setwd("D:/Projects/CASH-Human/")
# 1. Load packages
library(Seurat)
library(ggplot2)
library(dplyr)
source("./bin/function.R")


# 2. params
dir_for_NK_obj <- "./data/tmp/01.preprocess/sce_subset/NK_obj_clustered.rds"
## load data
NK_obj <- readRDS(dir_for_NK_obj)


# 3. markers
markers_list <- list(
  c("NCAM1", "FCGR3A",
    "IL32", 
    "IFI44L","CX3CR1",
    "BEND7", "ZNF90", "NACA2",
    "CCL3", "NFKBIA",
    "SPC25", "MKI67",
    "HSPA1A", "DNAJB1", 
    "CREM", "NR4A3",
    "GCSAM", "KLRC2",
    "CCL4", "CD160")
  )
markers <- unique(unlist(markers_list))


# 4. UMAP_Pct_Dotplot
pdf("./fig/tmp/01.preprocess/NK_UMAP_Pct_Dotplot.pdf", width = 14, height = 14)
resolution_list <- paste0("RNA_snn_res.", seq(0.1, 1, 0.1))
for (res in resolution_list) {
  # res <- "RNA_snn_res.0.5"
  ## UMAP
  p1 <- DimPlot(NK_obj, label = T, group.by = res, label.size = 6) + NoLegend() + ggtitle(paste0("UMAP-NK_", res))
  
  ## Dotplot
  p2 <-  yy_Dotplot(seuratObj = NK_obj,
                    genes = markers_list,
                    group.by = res) + ggtitle(paste0("Dotplot-NK_", res))
  
  ## Pct bar plot
  meta <- NK_obj@meta.data
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
Idents(NK_obj) <- "RNA_snn_res.0.4"
NK_obj <- RenameIdents(NK_obj,
                       "0" = "c28_NK_CD160",
                       "1" = "c32_NK_CREM",
                       "2" = "c33_NK_IL32",
                       "3" = "c29_NK_NFKBIA",
                       "4" = "c34_NK_CX3CR1",
                       "5" = "c30_NK_DNAJB1",
                       "6" = "c31_NK_CCL4")
NK_obj$sub_celltype <- Idents(NK_obj)
levels(NK_obj) <- unique(sort(levels(NK_obj)))
NK_obj$sub_celltype <- factor(NK_obj$sub_celltype, levels = unique(sort(levels(NK_obj))))

## check
DimPlot(NK_obj, group.by = c("RNA_snn_res.0.4", "sub_celltype"), label = T) + NoLegend()
saveRDS(NK_obj, "./data/obj/NK_obj.rds")
