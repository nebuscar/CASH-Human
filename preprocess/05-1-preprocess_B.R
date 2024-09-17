setwd("D:/Projects/CASH-Human/")
# 1. Load function
library(Seurat)
source("./bin/sc-pipeline.R")


# 2. params
dir_for_data <- "./data/tmp/01.preprocess/"
dir_for_B_obj <- "./data/tmp/01.preprocess/sce_subset/B_obj.rds"
dir_for_B_obj_clustered <- "./data/tmp/01.preprocess/sce_subset/B_obj_clustered.rds"
## load data
B_obj <- readRDS(dir_for_B_obj)


# 3. Integration and clustering
B_obj <- preprocess(B_obj)
B_obj <- clustering(B_obj)
saveRDS(B_obj, "./data/tmp/01.preprocess/sce_subset/B_obj_clustered.rds")


# 4. subset
## load data
B_obj_clustered <- readRDS(dir_for_B_obj_clustered)
Idents(B_obj_clustered) <- "RNA_snn_res.0.1"
idents_to_keep <- levels(B_obj_clustered)[!levels(B_obj_clustered) %in% c("2", "4")]
## create suerat obj
B_obj_clustered <- subset(B_obj_clustered, idents = c("2", "4"), invert = TRUE)
# B_obj_clustered <- subset(B_obj_clustered, idents = idents_to_keep)
meta <- B_obj_clustered@meta.data
counts <- B_obj_clustered@assays$RNA$counts
B_obj_clustered <- CreateSeuratObject(counts, meta.data =  meta)
## Re Integration and clustering
B_obj_clustered <- preprocess(B_obj_clustered)
B_obj_clustered <- clustering(B_obj_clustered)
saveRDS(B_obj_clustered, "./data/tmp/01.preprocess/sce_subset/B_obj_sub.rds")