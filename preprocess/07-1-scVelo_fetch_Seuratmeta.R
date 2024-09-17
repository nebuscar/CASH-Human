setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(stringr)

# 2. params
dir_for_obj <- "./data/obj/"
dir_for_obj_meta <- "./data/obj_meta/"

file.list <- list.files(dir_for_obj)
file.list <- file.list[file.list %in% c("B_obj.rds", 
                                        "Myeloid_obj.rds", 
                                        "NK_obj.rds" , 
                                        "T_obj.rds")]

# 3. fetch meta
if (dir.exists("./data/obj_meta")) {
  print("directory has been created")
}else{
  dir.create("./data/obj_meta")
  dir_for_obj_meta <- "./data/obj_meta/"
  print(paste0("dir_for_obj_meta: ", dir_for_obj_meta))
}

for (i_file in file.list) {
  print(i_file)
  i_cell <- str_split(i_file, "_", simplify = T)[,1]
  obj <- readRDS(file.path(dir_for_obj, i_file))
  obj <- AddMetaData(obj, Embeddings(obj, reduction = "umap.cca"), col.name = c("UMAP1", "UMAP2"))
  meta <- obj@meta.data
  meta$barcode <- rownames(meta)
  tmp <- meta[, c("barcode", 
                  "orig.ident", "group", "patient", "sub_celltype", "UMAP1", "UMAP2")]
  # tmp <- rownames_to_column(tmp, var = "barcode")
  tmp <- subset(tmp, group %in% "CASH")
  # write.csv(tmp, file.path(dir_for_obj_meta, paste0(i_cell, "_obj_meta.csv")),row.names = F)
  }
