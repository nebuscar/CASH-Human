setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(RColorBrewer)
library(ggplot2)


# 2. params
dir_for_T_obj <- "./data/obj/T_obj.rds"
## load data
T_obj <- readRDS(dir_for_T_obj)
T_obj@meta.data[c("UMAP1", "UMAP2")] <- Embeddings(T_obj, reduction = "umap.cca")
Idents(T_obj) <- "sub_celltype"
levels(T_obj) <- unique(sort(levels(T_obj)))
T_obj$sub_celltype <- factor(T_obj$sub_celltype, levels = levels(T_obj))
meta <- T_obj@meta.data

########## 4. Dotplot ########## 
## markers
markers_list <- list(
  T_Cells = c("CD3D", "CD4", "CD8A"),
  MAIT = c("SLC4A10", "KLRB1"),
  Effect_Memory = c("GZMK", "CXCR3"),
  Exhausted = c("TIGIT", "PDCD1", "TOX"),
  EMRA = c("CX3CR1", "GNLY", "GZMB"),
  gdT = c("TRDV2", "TRGV9"),
  Naïve = c("SELL", "CCR7", "LEF1"),
  Treg = c("FOXP3", "IL2RA")
)
## plot
dir.create("./fig/Fig04.T/", showWarnings = F)

p <- DotPlot(T_obj, features = markers_list, dot.scale = 4) +
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
ggsave("./fig/Fig04.T/Fig4B.Dotplot_T.pdf", plot = p, width = 16, height = 8, dpi = 300)
ggsave("./fig/Fig04.T/Fig4B.Dotplot_T.png", plot = p, width = 16, height = 8, dpi = 300)