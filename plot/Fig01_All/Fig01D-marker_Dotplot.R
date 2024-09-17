setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(RColorBrewer)
library(ggplot2)


# 2. params
dir_for_obj <- "./data/obj/All_obj.rds"
## load data
all_obj <- readRDS(dir_for_obj)
all_obj@meta.data[c("UMAP1", "UMAP2")] <- Embeddings(all_obj, reduction = "umap.cca")
Idents(all_obj) <- "celltype"
levels(all_obj) <- c("Monocytes",  "Macrophage", "cDC2", "cDC1", "pDC", "T cells", "NK", "B cells", "Plasma cells")
all_obj$celltype <- factor(all_obj$celltype, levels = levels(all_obj))
meta <- all_obj@meta.data

########## Fig01D.Dotplot ##########
## markers
markers_list <- list(
  Monocytes = c("FCN1", "VCAN", "S100A9", "S100A8"),
  Macrophages = c("CD68", "C1QA", "CD163", "MS4A7"),
  cDCs = c("CLEC9A", "IDO1", "FLT3", "LILRA4"),
  T_cells = c("CD3D", "GZMK", "IL7R", "CD8A"),
  NK_cells = c("NKG7", "NCAM1", "FCGR3A", "FGFBP2"),
  B_cells = c("CD79A", "MS4A1", "BANK1"),
  Plasma_cells = c("JCHAIN", "IGHG1", "TNFRSF17")
  # Mast_cells = c("MS4A2", 
  #                #"TPSAB2", 
  #                "CPA3"),
  # Cycling = c("MKI67", "RRM2", "TOP2A")
)


## plot
dir.create("./fig/Fig01.All", showWarnings = F)

p <- DotPlot(all_obj, features = markers_list, dot.scale = 4) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_radius(range = c(0, 4)) +
  scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression")) +
  theme(
    ## axis text
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
    axis.text.y = element_text(size = 12),
    text = element_text(size = 8),
    
    ## legend
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    
    ## panel
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    strip.text = element_text(size = 14),
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  labs(x = "", y = "")
print(p)
ggsave("./fig/Fig01.All/Fig1D.Dotplot_major.pdf", plot = p, width = 16, height = 8, dpi = 300)


