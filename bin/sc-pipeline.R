# 1. Load packages
load_required_packages <- function() {
  library(Seurat)
}

# 2. Load data
load_data <- function(file_path) {
  if (file.exists(file_path)) {
    return(readRDS(file_path))
  } else {
    stop(paste("File not found:", file_path))
  }
}

# 3. preprocess
preprocess <- function(sce) {
  sce[["RNA"]] <- split(sce[["RNA"]], f = sce$orig.ident)
  sce <- NormalizeData(sce, verbose = F)
  sce <- FindVariableFeatures(sce, verbose = F)
  sce <- ScaleData(sce, verbose = F)
  sce <- RunPCA(sce, verbose = F)
  
  # Intergration
  t <- table(sce@meta.data[["orig.ident"]])
  # k.weight <- 30
  if (min(t) < 100) {
    k.weight = 30
  } else  {
    k.weight = 100
  }
  print(paste0("k.weight = ", k.weight))
  sce <- IntegrateLayers(sce,
                         method = CCAIntegration,
                         orig.reduction = "pca",
                         new.reduction = "integrated.cca",
                         verbose = FALSE,
                         k.weight = k.weight)
  sce[["RNA"]] <- JoinLayers(sce@assays$RNA)
  sce <- FindNeighbors(sce, reduction = "integrated.cca", dims = 1:30)
  
  return(sce)
}

# 4. clustering
clustering <- function(sce, resolution_seq = seq(0.1, 1, 0.1)) {
  for (res in resolution_seq) {
    print(paste("Processing resolution:", res))
    sce <- FindClusters(sce, resolution = res, verbose = FALSE)
    sce <- RunUMAP(sce, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca", verbose = FALSE)
  }
  return(sce)
}

# 5. Save
save_sce <- function(sce, file_path) {
  saveRDS(sce, file_path)
  print(paste("File saved:", file_path))
}
