setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(dplyr)

# 2. params
dir_for_data <- "./data/tmp/01.preprocess/"
dir_for_all_obj <- "./data/tmp/01.preprocess/sce_sub_annotated.rds"
dir_for_all_obj_sub <- "./data/tmp/01.preprocess/sce_sub2_annotated.rds"
## load data
all_obj <- readRDS(dir_for_all_obj)

# 3. Add Annotation_major meta
all_obj@meta.data <- all_obj@meta.data %>%
  mutate(major_celltype = case_when(
    celltype %in% "T cells" ~ "T", 
    celltype %in% c("Macrophage", "Monocytes", "cDC1", "cDC2", "pDC", "Cycling") ~ "Myeloid", 
    celltype %in% c("B cells", "Plasma cells") ~ "B", 
    celltype %in% c("LSECs", "Hepatocytes", "HSCs") ~ "Parenchymal", 
    celltype %in% "NK" ~ "NK" 
  ))
Idents(all_obj) <- "major_celltype"
idents_to_keep <- c("T", "Myeloid", "B", "NK")
sce <- subset(all_obj, idents = idents_to_keep)
saveRDS(all_obj, "./data/tmp/01.preprocess/sce_sub2_annotated.rds")

# 4. Subset
dir.create(file.path(dir_for_data, "sce_subset"))
all_obj <- readRDS(dir_for_all_obj_sub)
for (cell in levels(all_obj)) {
  print(cell)
  if (cell %in% "Parenchymal") {
    next
  }
  obj_sub <- subset(sce, subset = major_celltype %in% cell)
  
  # create sce
  meta <- obj_sub@meta.data
  meta <- meta[, c(1:7, which(colnames(meta) %in% c("group", "patient", "major_celltype", "celltype")))]
  counts <- obj_sub@assays$RNA$counts
  obj_sub <- CreateSeuratObject(counts, meta.data =  meta)
  saveRDS(obj_sub, file.path(dir_for_data, "sce_subset", paste0(cell, "_obj.rds")))
}



