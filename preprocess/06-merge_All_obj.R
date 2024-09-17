setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
source("./bin/sc-pipeline.R")

# 2. params
dir_for_data <- "./data/obj/"
dir_for_B_obj <- "./data/obj/B_obj.rds"
dir_for_Myeloid_obj <- "./data/obj/Myeloid_obj.rds"
dir_for_NK_obj <- "./data/obj/NK_obj.rds"
dir_for_T_obj <- "./data/obj/T_obj.rds"
dir_for_All_obj_raw <- "./data/tmp/01.preprocess/sce_sub2_annotated.rds"
## load data
B_obj <- readRDS(dir_for_B_obj)
Myeloid_obj <- readRDS(dir_for_Myeloid_obj)
Nk_obj <- readRDS(dir_for_NK_obj)
T_obj <- readRDS(dir_for_T_obj)

# 3. merge All
All_obj_merge <- merge(B_obj, y = c(Myeloid_obj, Nk_obj, T_obj), add.cell.ids = NULL, meta.data = F, merge.dr = F)
rm(B_obj, Myeloid_obj, Nk_obj, T_obj)
meta1 <- All_obj_merge@meta.data

# 4. match
All_obj_raw <- readRDS(dir_for_All_obj_raw)
meta2 <- All_obj_raw@meta.data
meta <- meta2[match(rownames(meta1), rownames(meta2)),]
meta$sub_celltype <- meta1$sub_celltype
All_obj_match <- All_obj_raw[, match(rownames(meta), colnames(All_obj_raw))]
All_obj_match@meta.data <- meta
rm(All_obj_merge, All_obj_raw, meta1, meta2, meta)

# 5. create sce
meta <- All_obj_match@meta.data
counts <- All_obj_match@assays$RNA$counts
sce_sub <- CreateSeuratObject(counts, meta.data =  meta)

# 6. Re Integration and clustering
sce_sub <- preprocess(sce_sub)
sce_sub <- clustering(sce_sub)

levels(sce_sub) <- factor(levels = unique(sort(levels(sce_sub))))

# 7. add meta
sce_sub$sub_celltype_short <- stringr::str_split(sce_sub$sub_celltype, "_", simplify = T)[,1]

## Check UMAP
DimPlot(sce_sub, group.by = "sub_celltype", label = T)

saveRDS(sce_sub, "./data/obj/all_obj.rds")
