setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(RColorBrewer)
library(ggplot2)


# 2. params
dir_for_NK_obj <- "./data/obj/NK_obj.rds"
## load data
NK_obj <- readRDS(dir_for_NK_obj)
NK_obj@meta.data[c("UMAP1", "UMAP2")] <- Embeddings(NK_obj, reduction = "umap.cca")
Idents(NK_obj) <- "sub_celltype"
levels(NK_obj) <- unique(sort(levels(NK_obj)))
NK_obj$sub_celltype <- factor(NK_obj$sub_celltype, levels = levels(NK_obj))
meta <- NK_obj@meta.data

########## 4. Dotplot ########## 
## markers
markers_list <- list(
  "CD56[bright]CD16[-]NK" = c("NCAM1", "CD160", "NFKBIA", "CCL3", "CCL4", "NR4A3", "DNAJB1", "HSPA1A"),
  "CD56[dim]CD16[+]NK" = c("FCGR3A", "CREM", "IL32", "KLRC2", "CX3CR1")
)

## plot
dir.create("./fig/Fig03.NK/", showWarnings = F)

p <- DotPlot(NK_obj, features = markers_list, dot.scale = 4) +
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
ggsave("./fig/Fig03.NK/Fig3B.Dotplot_NK.pdf", plot = p, width = 16, height = 8, dpi = 300)
ggsave("./fig/Fig03.NK/Fig3B.Dotplot_NK.png", plot = p, width = 16, height = 8, dpi = 300)