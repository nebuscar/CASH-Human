setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(ggplot2)

# 2. params
dir_for_data <- "./data"
dir_for_raw_obj <- "./data/obj/CASH-Normal-harmony.rds"

## load data
sce <- readRDS(dir_for_raw_obj)

# 3. QC
Idents(sce) <- sce@meta.data$orig.ident
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, 
        raster = F, pt.size = 0)
table(sce@meta.data$orig.ident)

# 4. down sample
samples <- unique(sort(sce$orig.ident))
samples <- grep("CASH|Normal", samples, value = T)

sce_sub <- subset(sce, subset = nFeature_RNA > 800 & nFeature_RNA < 6000 & percent.mito < 5 & orig.ident %in% samples[-3])
table(sce_sub@meta.data$orig.ident)

# 5. save
dir.create(file.path(dir_for_data, "tmp", "01.preprocess"))
saveRDS(sce_sub, "./data/tmp/01.preprocess/sce_sub_qc.rds")