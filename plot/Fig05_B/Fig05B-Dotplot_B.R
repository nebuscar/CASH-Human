setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(RColorBrewer)
library(ggplot2)


# 2. params
dir_for_B_obj <- "./data/obj/B_obj.rds"
## load data
B_obj <- readRDS(dir_for_B_obj)
B_obj@meta.data[c("UMAP1", "UMAP2")] <- Embeddings(B_obj, reduction = "umap.cca")
Idents(B_obj) <- "sub_celltype"
levels(B_obj) <- unique(sort(levels(B_obj)))
B_obj$sub_celltype <- factor(B_obj$sub_celltype, levels = levels(B_obj))
meta <- B_obj@meta.data

########## 4. Dotplot ########## 
## markers
markers_list <- list(
  Plasma_cells = c("CD38", "MZB1", "PRDM1", "IRF4", "XBP1", "IGHG1", "IGHA1"),
  Naive_B_cells = c("FCER2", "TCL1A", "IL4R", "CD72", "BACH2", "IGHD", "IGHM"),
  Memory_B_cells = c("MS4A1", "LTB", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DQA1", "TNFRSF13B", "TXNIP", "GPR183"),
  atypical_B_cells = c("FCRL5", "ITGAX", "TBX21", "CR2", "HCK", "FCRL3", "FGR"),
  Germinal_center_B_cells = c("MKI67", "STMN1", "HMGB2", "TOP2A")
)

## plot
dir.create("./fig/Fig05.B/", showWarnings = F)

p <- DotPlot(B_obj, features = markers_list, dot.scale = 4) +
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
    strip.text = element_text(size = 10),
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed", linewidth = 0.5),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5)
  ) +
  labs(x = "", y = "")
print(p)
ggsave("./fig/Fig05.B/Fig2B.Dotplot_B.pdf", plot = p, width = 10, height = 10, dpi = 300)
ggsave("./fig/Fig05.B/Fig2B.Dotplot_B.png", plot = p, width = 10, height = 10, dpi = 300)

