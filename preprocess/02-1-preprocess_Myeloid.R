setwd("D:/Projects/CASH-Human/")
# 1. Load function
library(Seurat)
source("./bin/sc-pipeline.R")


# 2. params
dir_for_data <- "./data/tmp/01.preprocess/"
dir_for_Myeloid_obj <- "./data/tmp/01.preprocess/sce_subset/Myeloid_obj.rds"
dir_for_Myeloid_obj_clustered <- "./data/tmp/01.preprocess/sce_subset/Myeloid_obj_clustered.rds"
dir_for_Myeloid_obj_sub <- "./data/tmp/01.preprocess/sce_subset/Myeloid_obj_sub.rds"
## load data
Myeloid_obj <- load_data(dir_for_Myeloid_obj)


# 3. Integration and clustering
Myeloid_obj <- preprocess(Myeloid_obj)
Myeloid_obj <- clustering(Myeloid_obj)
saveRDS(Myeloid_obj, "./data/tmp/01.preprocess/sce_subset/Myeloid_obj_clustered.rds")


# 4. Subset1
## load data
Myeloid_obj_clustered <- load_data(dir_for_Myeloid_obj_clustered)
Idents(Myeloid_obj_clustered) <- "RNA_snn_res.0.8"
idents_to_keep <- levels(Myeloid_obj_clustered)[!levels(Myeloid_obj_clustered) %in% c("7", "8", "10")]
## create suerat obj_sub
Myeloid_obj_clustered <- subset(Myeloid_obj_clustered, idents = c("7", "8", "10"), invert = TRUE)
# Myeloid_obj_clustered <- subset(Myeloid_obj_clustered, idents = idents_to_keep)
meta <- Myeloid_obj_clustered@meta.data
counts <- Myeloid_obj_clustered@assays$RNA$counts
Myeloid_obj_clustered <- CreateSeuratObject(counts, meta.data =  meta)
## Re Integration and clustering
Myeloid_obj_clustered <- preprocess(Myeloid_obj_clustered)
Myeloid_obj_clustered <- clustering(Myeloid_obj_clustered)
saveRDS(Myeloid_obj_clustered, "./data/tmp/01.preprocess/sce_subset/Myeloid_obj_sub.rds")


# 5. Subset2
## load data
Myeloid_obj_sub <- load_data(dir_for_Myeloid_obj_sub)
Idents(Myeloid_obj_sub) <- "RNA_snn_res.0.8"
idents_to_keep <- levels(Myeloid_obj_sub)[!levels(Myeloid_obj_sub) %in% "8"]
## create suerat obj_sub2
Myeloid_obj_sub <- subset(Myeloid_obj_sub, idents = idents_to_keep)
meta <- Myeloid_obj_sub@meta.data
counts <- Myeloid_obj_sub@assays$RNA$counts
Myeloid_obj_sub <- CreateSeuratObject(counts, meta.data =  meta)
## Re Integration and clustering
Myeloid_obj_sub <- preprocess(Myeloid_obj_sub)
Myeloid_obj_sub <- clustering(Myeloid_obj_sub)
saveRDS(Myeloid_obj_sub, "./data/tmp/01.preprocess/sce_subset/Myeloid_obj_sub2.rds")
