setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(dplyr)


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
    # panel.border = element_rect(color = "black", fill = NA, size = 1),
    
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
       width = 5, height = 5, dpi = 300)
ggsave("./fig/Fig01.All/Fig1C.Dimplot_major.png", plot = p, 
       width = 5, height = 5, dpi = 300)

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
p <- DotPlot(all_obj, features = markers_list, dot.scale = 4) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_radius(range = c(0, 4)) +
  scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression")) +
  labs(x = "", y = "") +
  theme(
    ## axis text
    axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 6),
    text = element_text(size = 0),
    
    ## facet
    strip.text = element_text(size = 8),
    strip.background = element_blank(),
    
    ## panel
    plot.margin = unit(c(0, 1, 0, 0), "char"),
    panel.spacing = unit(x = 0.2, units = "lines"),
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed", linewidth = 0.3),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    
    ## legend
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.15, "inch"),
    legend.position = "bottom")
print(p)
ggsave("./fig/Fig01.All/Fig1D.Dotplot_major.pdf", plot = p, width = 8, height = 4, dpi = 300)
ggsave("./fig/Fig01.All/Fig1D.Dotplot_major.png", plot = p, width = 8, height = 4, dpi = 300)

########## Fig01E.Pct_bar_plot ##########
# 3. plot
tmp <- meta[, c("orig.ident", "celltype")]

num_celltype <- tmp %>% 
  group_by(celltype, orig.ident) %>% 
  summarise(n_cluster = n(), .groups = "drop")
num_total <- tmp %>% 
  group_by(orig.ident)%>% 
  summarise(n_total = n(), .groups = "drop")
plot_df <- num_celltype %>%
  left_join(num_total, by = "orig.ident") %>%
  mutate(pct = (n_cluster / n_total) * 100)
plot_df$orig.ident <- factor(plot_df$orig.ident, levels = c("CASH_Patient1", "CASH_Patient2", "CASH_Patient3", "CASH_Patient4", "CASH_Patient5",
                                                            "Normal_Patient1", "Normal_Patient2", "Normal_Patient3", "Normal_Patient4", "Normal_Patient5",
                                                            "Normal_Patient6", "Normal_Patient7", "Normal_Patient8", "Normal_Patient9", "Normal_Patient10"))

p <- ggplot(plot_df, aes(x = orig.ident, y = pct, fill =  celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  cowplot::theme_cowplot() +
  scale_fill_manual(values = celltype_color_panel) + 
  labs(x = "", y = "Proportion", title = "") +
  theme(
    ## axis
    axis.text.x = element_text(size = 14, angle = 90),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    
    ## legend
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.position = "right",
    
    ## panel
    axis.line.x = element_line(linewidth = 1),
    axis.line.y = element_line(linewidth = 1)) +
  coord_flip()
print(p)
ggsave("./fig/Fig01.All/Fig1E.Pct_bar_plot.pdf", plot = p, width = 8, height = 8, dpi = 300)
ggsave("./fig/Fig01.All/Fig1E.Pct_bar_plot.png", plot = p, width = 8, height = 8, dpi = 300)



