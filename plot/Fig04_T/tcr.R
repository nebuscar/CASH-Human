setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(dplyr)
library(ggplot2)


# 2. params
dir_for_T_obj <- "./data/obj/T_obj.rds"
dir_for_meta_tcr <- "./data/tcr/meta_tcr.rds"

Clonalsize_color_panel <- c(
  "Hyperexpanded (0.1 < X <= 1)" = "#78290f",
  "Large (0.01 < X <= 0.1)" = "#78290f",
  "Medium (0.001 < X <= 0.01)" = "#ff7d10",
  "Small (1e-04 < X <= 0.001)" = "#ffecd1",
  "Rare (0 < X <= 1e-04)" = "#FFF9DA",
  "None ( < X <= 0)" = "#BFBFBF")

Subset_color_panel <- c(
  # T cells
  "c14_CD8_MAIT" = "#5A7B8F",
  "c15_CD8_MAIT" = "#EE4C97",
  "c16_CD8_MAIT" = "#3D806F",
  "c17_CD8_" = "#A08634",
  "c18_CD8_" = "#F37C95",
  "c19_CD8_" = "#608541",
  "c20_CD8_" = "#7D4E57",
  "c21_CD8_" = "#BC3C29",
  "c22_CD8_" = "#958056",
  "c23_CD4_Naïve" = "#9FAFA3",
  "c24_CD4_" = "#6F99AD",
  "c25_CD4_" = "#0072B5",
  "c26_CD4_Treg" = "#CFC59A",
  "c27_gdT" = "#E18727"
  )

## load data
T_obj <- readRDS(dir_for_T_obj)
Idents(T_obj) <- "sub_celltype"
levels(T_obj) <- unique(sort(levels(T_obj)))
T_obj$sub_celltype <- factor(T_obj$sub_celltype, levels = levels(T_obj))
meta_tcr <- readRDS(dir_for_meta_tcr)

## match meta
meta_match <- meta_tcr[match(rownames(T_obj@meta.data), rownames(meta_tcr)),]
T_obj@meta.data$cloneSize <- meta_match$cloneSize


########## 3.Count of TCR clones ##########
tmp <- T_obj@meta.data[, c("group", "orig.ident", "cloneSize")]
tmp <- na.omit(tmp)
tmp <- table(tmp) %>% 
  as.data.frame() %>%
  group_by(group, orig.ident) %>%
  summarise(Count = sum(Freq))
colnames(tmp)[2] <- "Donor"

ggplot(tmp, aes(x = Donor, y = Count, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("CASH" = "#0f7b9f", "Normal" = "#d83215"), name = "CloneSize") +
  cowplot::theme_cowplot() +
  xlab("") + ylab("Count") +
  theme(
    ## axis
    axis.text.x = element_text(size = 7, angle = 0),
    text = element_text(size = 7),
    axis.line = element_line(size = 0.3),
    axis.ticks = element_line(size = 0.3),
    axis.text.y = element_text(size = 7, angle = 0),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    
    ## panel
    panel.border = element_blank(),
    
    ## panel
    legend.position = "none") + 
  guides(fill = guide_legend(ncol = 3)) +
  coord_flip()

########## UMAP-clonSize##########
colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)
T_obj@meta.data[c("UMAP1", "UMAP2")] <- Embeddings(T_obj, reduction = "umap.cca")
tmp <- T_obj@meta.data
p <- ggplot(tmp, aes(x = UMAP1, y = UMAP2, color = cloneSize)) +
  ## panel
  geom_point(size = 0.8, shape = 16, stroke = 0) +
  theme_void() +
  ## split
  facet_wrap(~ group) +
  ## color
  scale_color_manual(values = rev(colorblind_vector[c(1,3,4,5,7)]), name = '') +
  theme(
    aspect.ratio = 1,
    ## facet
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.background = element_rect(color = NA, fill = NA, size = 1),
    
    ## legend
    legend.text = element_text(size = 12),
    legend.spacing.y = unit(0, 'cm'),
    legend.key.height = unit(0,"cm"),
    legend.box.spacing = unit(0, 'cm'),
    legend.position = "right") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 5, alpha = 1)))
## Dimplot
p <- DimPlot(T_obj, group.by = "cloneSize", split.by = "group",
              pt.size = 0.8, order = T, shuffle  = T, raster=FALSE, reduction = "umap.cca") +
  scale_color_manual(values=rev(colorblind_vector[c(1,3,4,5,7)]))
print(p)
ggsave("./fig/Fig04.T/Fig4C.Dimplot_TCR.pdf", plot = p, 
       width = 10, height = 6, dpi = 300)
ggsave("./fig/Fig04.T/Fig4C.Dimplot_TCR.png", plot = p, 
       width = 10, height = 6, dpi = 300)

########## 4.Proportion ##########
tmp <- T_obj@meta.data
tmp$cloneSize <- as.character(tmp$cloneSize)
# tmp$cloneSize[is.na(tmp$cloneSize)] <- "None ( < X <= 0)"

tmp <- tmp[, c("sub_celltype", "group", "cloneSize")]
tmp$cluster <- unique(sort(tmp$cluster))
tmp <- table(tmp) %>% 
  as.data.frame() 
colnames(tmp)[1] <- "cluster" 

p <- ggplot(tmp, aes(x = cluster, y = Freq, fill = as.character(cloneSize))) +
  ## panel
  geom_bar(stat = "identity", position = "fill") +
  cowplot::theme_cowplot() +
  ## color
  scale_fill_manual(values = Clonalsize_color_panel, name = "CloneSize") +
  ## split
  facet_grid(~group, scales = "free") +
  labs(x = "", y = "Proportion") +
  theme(
    ## axis
    axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.line = element_line(size = 0.5),
    axis.ticks = element_line(size = 0.5),
    
    ## panel
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    
    ## facet
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.background = element_rect(color = NA, fill = NA, size = 1),
    
    ## legend
    legend.position = "bottom") + 
  guides(fill = guide_legend(ncol = 3))
print(p)
ggsave("./fig/Fig04.T/Fig4D.TCR_propotion.pdf", width = 10, height = 8, dpi = 300)
ggsave("./fig/Fig04.T/Fig4D.TCR_propotion.png", width = 10, height = 8, dpi = 300)
