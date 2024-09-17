setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(RColorBrewer)
library(ggplot2)


# 2. params
dir_for_obj <- "./data/obj/All_obj.rds"
celltype_color_panel <- c(
  "T cells" = "#98c9dd",       
  "NK" = "#a6d38e",           
  "Monocytes" = "#37a849",     
  "Macrophage" = "#207cb5", 
  "B cells" = "#f69595",      
  "Plasma cells" = "#eb2a2a", 
  "pDC" = "#9467bd",         
  "cDC1" = "#fcba71",         
  "cDC2" = "#f78200" 
)
## load data
all_obj <- readRDS(dir_for_obj)
all_obj@meta.data[c("UMAP1", "UMAP2")] <- Embeddings(all_obj, reduction = "umap.cca")
Idents(all_obj) <- "celltype"
levels(all_obj) <- c("Monocytes",  "Macrophage", "cDC2", "cDC1", "pDC", "T cells", "NK", "B cells", "Plasma cells")
all_obj$celltype <- factor(all_obj$celltype, levels = levels(all_obj))
meta <- all_obj@meta.data

########## Fig01C.Dmplot ########## 
dir.create("./fig/Fig01.All", showWarnings = F)

## UMAP
p <- ggplot(meta, aes(x = UMAP1, y = UMAP2, color = celltype)) +
  geom_point(size = 0.4, shape = 16, stroke = 0) +
  theme_void() +
  scale_color_manual(values = celltype_color_panel, name = '') +
  theme(
    ## border
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    
    ## legend
    aspect.ratio = 1, 
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0, 'cm'),
    legend.key.height = unit(0,"cm"),
    legend.box.spacing = unit(0, 'cm')) +
  guides(color = guide_legend(
    # ncol = 2,
    override.aes = list(size = 5, alpha = 1)))
print(p)
ggsave("./fig/Fig01.All/Fig1C.Dimplot_major.pdf", plot = p, 
       width = 8, height = 8, dpi = 300)
ggsave("./fig/Fig01.All/Fig1C.Dimplot_major.png", plot = p, 
       width = 8, height = 8, dpi = 300)
