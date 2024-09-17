setwd("D:/Projects/CASH-Human/")
# 1. Load function
library(Seurat)
source("./bin/sc-pipeline.R")


# 2. params
dir_for_data <- "./data/tmp/01.preprocess/"
dir_for_T_obj <- "./data/tmp/01.preprocess/sce_subset/T_obj.rds"
dir_for_T_obj_clustered <- "./data/tmp/01.preprocess/sce_subset/T_obj_clustered.rds"
dir_for_T_obj_sub <- "./data/tmp/01.preprocess/sce_subset/T_obj_sub.rds"
## load data
T_obj <- load_data(dir_for_T_obj)


# 3. Integration and clustering
T_obj <- preprocess(T_obj)
T_obj <- clustering(T_obj)
saveRDS(T_obj, "./data/tmp/01.preprocess/sce_subset/T_obj_clustered.rds")

# 4. Subset1
## load data
T_obj_clustered <- load_data(dir_for_T_obj_clustered)
Idents(T_obj_clustered) <- "RNA_snn_res.0.6"
idents_to_keep <- levels(T_obj_clustered)[!levels(T_obj_clustered) %in% c("3", "12", "13")]
## create suerat obj_sub
T_obj_clustered <- subset(T_obj_clustered, idents = idents_to_keep)
meta <- T_obj_clustered@meta.data
counts <- T_obj_clustered@assays$RNA$counts
T_obj_clustered <- CreateSeuratObject(counts, meta.data =  meta)
## Re Integration and clustering
T_obj_clustered <- preprocess(T_obj_clustered)
T_obj_clustered <- clustering(T_obj_clustered)
saveRDS(T_obj_clustered, "./data/tmp/01.preprocess/sce_subset/T_obj_sub.rds")

# 5. Subset2
## load data
T_obj_sub <- load_data(dir_for_T_obj_sub)
Idents(T_obj_sub) <- "RNA_snn_res.0.6"
idents_to_keep <- levels(T_obj_sub)[!levels(T_obj_sub) %in% "15"]
## create suerat obj_sub
T_obj_sub <- subset(T_obj_sub, idents = idents_to_keep)
meta <- T_obj_sub@meta.data
counts <- T_obj_sub@assays$RNA$counts
T_obj_sub <- CreateSeuratObject(counts, meta.data =  meta)
## Re Integration and clustering
T_obj_sub <- preprocess(T_obj_sub)
T_obj_sub <- clustering(T_obj_sub)
saveRDS(T_obj_sub, "./data/tmp/01.preprocess/sce_subset/T_obj_sub2.rds")