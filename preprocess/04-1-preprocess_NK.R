setwd("D:/Projects/CASH-Human/")
# 1. Load function
library(Seurat)
source("./bin/sc-pipeline.R")


# 2. params
dir_for_data <- "./data/tmp/01.preprocess/"
dir_for_NK_obj <- "./data/tmp/01.preprocess/sce_subset/NK_obj.rds"
## load data
NK_obj <- load_data(dir_for_NK_obj)


# 3. Integration and clustering
NK_obj <- preprocess(NK_obj)
NK_obj <- clustering(NK_obj)
saveRDS(NK_obj, "./data/tmp/01.preprocess/sce_subset/NK_obj_clustered.rds")