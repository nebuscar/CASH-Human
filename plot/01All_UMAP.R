setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(dplyr)
library(ggplot2)

# 2. params
dir_for_obj <- "./data/obj/"

## color_panel
Major_color_panel = c(
  "Myeloid" = "#E18727FF",
  "T" =  "#0072B5FF",
  "NK" = "#7876B1FF",
  "B" = "#20854EFF"
)

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

Subset_color_panel <- c(
  # Myeloid cells
  "c01_ Monocyte_CD14CD16" = "#5A7B8F",
  "c02_ Monocyte_CD14" = "#3D806F",
  "c03_ Monocyte_CD16" = "#BC3C29",
  "c04_ Monocyte_CD14" = "#9FAFA3",
  "c05_ Macrophage_" = "#F37C95",
  "c06_ Macrophage_" = "#A08634",
  "c07_ Macrophage_ MARCO" = "#7D4E57",
  "c08_ cDC1" = "#608541",
  "c09_cDC2" = "#EE4C97",
  "c10_cDC2" = "#958056",
  "c11_cDC2" = "#6F99AD",
  "c12_pDC" = "#0072B5",
  "c13_Actived_DC" = "#CFC59A",
  
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
  "c27_gdT" = "#E18727",
  
  # NK cells
  "c28_NK_CD160" = "#E18727FF", 
  "c29_NK_NFKBIA" = "#7876B1FF",
  "c30_NK_DNAJB1" = "red",
  "c31_NK_CCL4" = "#7D4E57",
  "c32_NK_CREM" = "#0072B5FF", 
  "c33_NK_IL32" = "#20854EFF",
  "c34_NK_CX3CR1" = "black", 


  
  # B cells
  "c35_PC" = "#E18727FF", 
  "c36_Bn_TCL1A" = "#0072B5FF", 
  "c37_classical-Bm_GRP183" = "#20854EFF",
  "c38_classical-Bm_TXNIP" = "#7876B1FF",
  "c39_ABC_FGR" = "black", 
  "c40_Bgc_MKI67" = "red"
  )


########## 3. Dimplot-major_celltype ########## 
## 3-1. loda data
all_obj <- readRDS(file.path(dir_for_obj, "All_obj.rds"))
all_obj@meta.data[c("UMAP1", "UMAP2")] <- Embeddings(all_obj, reduction = "umap.cca")
Idents(all_obj) <- "celltype"
levels(all_obj) <- c("Monocytes",  "Macrophage", "cDC2", "cDC1", "pDC", "T cells", "NK", "B cells", "Plasma cells")
all_obj$celltype <- factor(all_obj$celltype, levels = levels(all_obj))
meta <- all_obj@meta.data
## 3-2. plot

## UMAP
p <- ggplot(meta, aes(x = UMAP1, y = UMAP2, color = celltype)) +
  geom_point(size = 0.4, shape = 16, stroke = 0) +
  theme_void() +
  scale_color_manual(values = celltype_color_panel, name = '') +
  theme(
    ## legend
    aspect.ratio = 1, 
    legend.position = "right",
    legend.text = element_text(size = 14),
    legend.spacing.y = unit(0, 'cm'),
    legend.key.height = unit(0,"cm"),
    legend.box.spacing = unit(0, 'cm')) +
  guides(
    color = guide_legend(
      # ncol = 2,
      override.aes = list(size = 5, alpha = 1)))

ggsave("./fig/01Dimplot/all_obj_Dimplot.pdf", plot = p, 
       width = 8, height = 8, dpi = 300)


########## 4. Dimplot-sub_celltype ########## 
file_list <- list.files(dir_for_obj)
file_list <- file_list[!file_list %in% c("CASH-Normal-harmony.rds", "all_obj.rds")]

for (i_file in file_list) {
  print(i_file)
  cell <- stringr::str_split(i_file, "_", simplify = T)[,1]
  sce_sub <- readRDS(file.path(dir_for_obj, i_file))
  sce_sub@meta.data[c("UMAP1", "UMAP2")] <- Embeddings(sce_sub, reduction = "umap.cca")
  
  levels(sce_sub) <- unique(sort(levels(sce_sub)))
  sce_sub$sub_celltype <- factor(sce_sub$sub_celltype, levels = unique(sort(levels(sce_sub))))
  
  ## plot
  p <- ggplot(sce_sub@meta.data, aes(x = UMAP1, y = UMAP2, color = sub_celltype)) +
    geom_point(size = 1, shape = 16, stroke = 0) +
    theme_void() +
    scale_color_manual(values = Subset_color_panel, name = '') +
    theme(
      ## legend
      aspect.ratio = 1, 
      legend.position = "right",
      legend.text = element_text(size = 14),
      legend.spacing.y = unit(0, 'cm'),
      legend.key.height = unit(0,"cm"),
      legend.box.spacing = unit(0, 'cm')) +
    guides(
      color = guide_legend(
        # ncol = 2,
        override.aes = list(size = 5, alpha = 1)))
  print(p)
  ggsave(paste0("./fig/01Dimplot/", cell, "_Dimplot.pdf"), plot = p, width = 8, height = 8, dpi = 300)
}

