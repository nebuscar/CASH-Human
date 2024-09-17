setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(RColorBrewer)
library(ggplot2)


# 2. params
dir_for_obj <- "./data/obj/Myeloid_obj.rds"
## load data
Myeloid_obj <- readRDS(dir_for_obj)
Myeloid_obj@meta.data[c("UMAP1", "UMAP2")] <- Embeddings(Myeloid_obj, reduction = "umap.cca")
Idents(Myeloid_obj) <- "sub_celltype"
levels(Myeloid_obj) <- unique(sort(levels(Myeloid_obj)))
Myeloid_obj$sub_celltype <- factor(Myeloid_obj$sub_celltype, levels = levels(Myeloid_obj))
meta <- Myeloid_obj@meta.data

########## 4. Dotplot ########## 
## markers
markers_list <- list(
  Monocytes = c("FCN1","VCAN",'S100A9','S100A8',"CD14","FCGR3A"),
  Macrophages = c("CD68","C1QA","CD163","MARCO","FOLR2","LYVE1"),
  DCs = c("FLT3",'IDO1',"CLEC9A",'CD1C',"LILRA4","CLEC4C","JCHAIN","LAMP3","CCR7")
)

## plot
dir.create("./fig/Fig02.Myeloid/", showWarnings = F)

p <- DotPlot(Myeloid_obj, features = markers_list, dot.scale = 4) +
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
ggsave("./fig/Fig02.Myeloid/Fig2B.Dotplot_Myeloid.pdf", plot = p, width = 8, height = 8, dpi = 300)
ggsave("./fig/Fig02.Myeloid/Fig2B.Dotplot_Myeloid.png", plot = p, width = 8, height = 8, dpi = 300)

