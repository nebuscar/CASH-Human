setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(ggplot2)
library(stringr)

# 2. params
dir_for_obj <- "./data/tmp/01.preprocess/sce_sub_clustered.rds"
source("./bin/function.R")
# dir_for_markers <- "./data/markers.xlsx"

## load data
all_obj <- readRDS(dir_for_obj)


# 3. markers
markers_list <- list(
  NK = c("FCGR3A", "NCAM1", "NKG7"),
  CD8CD4 = c("CD3D", "CD8A", "CD4"),
  Treg = c("FOXP3", "IL2RA"),
  Naive = c("SELL", "CCR7", "LEF1"),
  Cytotoxicity = c("GZMB", "PRF1"),
  Exhaustion = c("LAG3", "TIGIT", "PDCD1", "TOX"),
  Effect = c("CX3CR1", "GNLY"),
  EffectMemory = c("GZMK", "CXCR3"),
  ResidentMemory = c("CD6"),
  MAIT = c("SLC4A10"),
  Macrophage = c("C1QA", "CD163", "SPP1", "LYVE1", "FOLR2", "MARCO"),
  Monocyte = c("FCN1", "VCAN", "CD14", "FCGR3A"),
  cDC1 = c("CLEC9A"),
  cDC2 = c("CD1C"),
  ActivedDC = c("LAMP3", "CCR7"),
  pDC = c("LILRA4", "CLEC4C", "JCHAIN"),
  Kupffer = c("CLEC4F", "LYVE1"),
  Hepatocytes = c("ALB", "HAMP", "ARG1", "APOC3", "FABP1", "APOA1"),
  HSCs = c("ACTA2", "COL1A1", "COL1A2", "COL3A1", "FGF7", "MME"),
  Cholangiocytes = c("SOX9", "EPCAM", "KRT19"),
  NK_DotPlot = c("AREG", "KLRB1", "NCR1"),
  LSECs = c("CALCRL", "PECAM1", "VWF"),
  Macrophages = c("CD68", "C1QA", "C1QB"),
  Monocytes_DotPlot = c("CD68", "S100A9", "S100A8", "MMP19"),
  B_cells = c("CD19", "MS4A1", "CD79A"),
  Plasma_cells = c("IGHG1", "MZB1", "SDC1", "CD79A"),
  DCs = c("FLT3", 
          #"CD123", 
          "CD1E", "CD1C", "LAMP3", "IDO1", "IDO2"),
  Neutrophil = c("MMP9")
)
# markers <- readxl::read_excel(dir_for_markers)
markers_list <- unique(unlist(markers_list))


# 4. Add group
all_obj$group <- str_split(all_obj$orig.ident, "_", simplify = T)[, 1]
all_obj$patient <- str_split(all_obj$orig.ident, "_", simplify = T)[, 2]
table(all_obj$group, all_obj$patient)


# 5. UMAP_Pct_Dotplot
pdf("./fig/tmp/01.preprocess/All_UMAP_Pct_Dotplot.pdf", width = 14, height = 14)
resolution_list <- paste0("RNA_snn_res.", seq(0.1, 1, 0.1))
for (res in resolution_list) {
  # res <- "RNA_snn_res.0.5"
  ## UMAP
  p1 <- DimPlot(all_obj, label = T, group.by = res, label.size = 6) + NoLegend() + ggtitle(paste0("UMAP-All_", res))
  
  ## Dotplot
  p2 <-  yy_Dotplot(seuratObj = all_obj,
                    genes = markers_list,
                    group.by = res) + ggtitle(paste0("Dotplot-All_", res))
  
  ## Pct bar plot
  meta <- all_obj@meta.data
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

# 6. annotation
Idents(all_obj) <- "RNA_snn_res.0.5"
all_obj <- RenameIdents(all_obj,
                        "0"="T cells",
                        "1"="T cells",
                        "2"="Macrophage",
                        "3"="NK",
                        "4"="T cells",
                        "5"="NK",
                        "6"="T cells",
                        "7"="Monocytes",
                        "8"="T cells",
                        "9"="T cells",
                        "10"="T cells",
                        "11"="LSECs",
                        "12"="B cells",
                        "13"="Plasma cells",
                        "14"="cDC1",
                        "15"="Cycling",
                        "16"="Hepatocytes",
                        "17"="cDC2",
                        "18"="pDC",
                        "19"="HSCs")

all_obj@meta.data$celltype <- Idents(all_obj)
Idents(all_obj) <- "celltype"
saveRDS(all_obj, "./data/tmp/01.preprocess/sce_sub_annotated.rds")


