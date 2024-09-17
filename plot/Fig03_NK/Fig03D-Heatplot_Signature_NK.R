setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)


# 2. params
fdir <- ""
dir_for_NK_obj <- "./data/obj/NK_obj.rds"
dir_for_cellRankings <- "./data/AUCell/AUC_NK.rds"
## load data
NK_obj <- readRDS(dir_for_NK_obj)


########## Heatplot_Signature ##########
## markers
markers_list <- list(
  inflammatory = c('IL18', 'IL15', 'IL7', 'IL6', 'IL1B', 'CXCL9', 'CXCL10', 'CCL5', 'CCL4', 'CCL3', 'CCL2'),
  stresss = c("ZFP36L1", "ZFP36", "ZFAND2A", "UBC", "SOCS3", "SLC2A3", "RGS2", "NFKBIZ", "NFKBIA", "JUNB", 
              "JUN", "IER2", "HSPH1", "HSPB1", "HSPA6", "HSPA1B", "HSPA1A", "HSP90B1", "HSP90AB1", "HSP90AA1", 
              "HIF1A", "FOSB", "FOS", "EGR1", "DUSP1", "DNAJB1", "CALU", "BAG3"),
  cytotoxicity = c("CTSW", "PRF1", "GNLY", "GZMK", "GZMM", "GZMH", "GZMB", "GZMA"))
markers <- unique(unlist(markers_list))

plot_df = FetchData(NK_obj,vars = c(markers,"sub_celltype"),slot = "data")
plot_df = plot_df %>% group_by(sub_celltype) %>% summarise_all(function(x) mean(x = expm1(x = x)))
plot_df = plot_df %>% tibble::column_to_rownames(var = "sub_celltype")
plot_df = scale(plot_df)
range(plot_df)

new_order <- c("c28_NK_CD160", "c29_NK_NFKBIA", "c30_NK_DNAJB1", "c31_NK_CCL4", "c32_NK_CREM", "c33_NK_IL32", "c34_NK_CX3CR1")
plot_df <- plot_df[new_order, ]

col_fun <- circlize::colorRamp2(c(-2.5, 0, 2.5), c("#0F7B9F", "white", "#D83215"))

## plot
p <- ComplexHeatmap::Heatmap(
  plot_df,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  column_names_side = "bottom",
  column_names_rot = 90,
  name = "Z-score",
  heatmap_legend_param = list(
    legend_direction = "vertical",
    title_position = "topcenter",
    legend_width = unit(0.05, "inch"),
    title_gp = gpar(fontsize = 7),
    labels_gp = gpar(fontsize = 6)),
  column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp = grid::gpar(fontsize = 6),
  width = ncol(plot_df) * unit(0.1, "inch"),
  height = nrow(plot_df) * unit(0.1, "inch")
  )
pdf("./fig/Fig03.NK/Fig3D.Heatplot_Signature_NK.pdf", width = 6, height = 6)
print(p)
dev.off()

png("./fig/Fig03.NK/Fig3D.Heatplot_Signature_NK.png", width = 6, height = 6, units = "in", res = 300)
print(p)
dev.off()
