setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(RColorBrewer)
library(ggplot2)


# 2. params
dir_for_NK_obj <- "./data/obj/NK_obj.rds"
Subset_color_panel <- c(
  # Myeloid cells
  "c01_Monocyte_CD14CD16" = "#5A7B8F",
  "c02_Monocyte_CD14" = "#3D806F",
  "c03_Monocyte_CD16" = "#BC3C29",
  "c04_Monocyte_CD14" = "#9FAFA3",
  "c05_Macrophage_" = "#F37C95",
  "c06_Macrophage_" = "#A08634",
  "c07_Macrophage_MARCO" = "#7D4E57",
  "c08_cDC1" = "#608541",
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

## load data
NK_obj <- readRDS(dir_for_NK_obj)
NK_obj@meta.data[c("UMAP1", "UMAP2")] <- Embeddings(NK_obj, reduction = "umap.cca")
Idents(NK_obj) <- "sub_celltype"
levels(NK_obj) <- unique(sort(levels(NK_obj)))
NK_obj$sub_celltype <- factor(NK_obj$sub_celltype, levels = levels(NK_obj))
meta <- NK_obj@meta.data

########## 3. UMAP ########## 
dir.create("./fig/Fig03.NK/", showWarnings = F)

## UMAP
p <- ggplot(meta, aes(x = UMAP1, y = UMAP2, color = sub_celltype)) +
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
  guides(color = guide_legend(
    # ncol = 2,
    override.aes = list(size = 5, alpha = 1)))
print(p)
ggsave("./fig/Fig03.NK/Fig3A.Dimplot_NK.pdf", plot = p, 
       width = 8, height = 8, dpi = 300)
ggsave("./fig/Fig03.NK/Fig3A.Dimplot_NK.png", plot = p, 
       width = 8, height = 8, dpi = 300)