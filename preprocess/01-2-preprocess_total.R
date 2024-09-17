setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)

# 2. params
dir_for_obj <- "./data/tmp/01.preprocess/sce_sub_qc.rds"
dir_for_gtf <- "./data/gtf.rds"

## load data
sce_sub <- readRDS(dir_for_obj)
meta <- sce_sub@meta.data
counts <- sce_sub@assays$RNA@counts

# 3. convert Ensembl to symbol
gtf <- readRDS(dir_for_gtf)
gene_mapping <- gtf[match(rownames(counts), gtf$gene_id),]
gene_mapping <- na.omit(gene_mapping)
gene_mapping_unique <- gene_mapping[!duplicated(gene_mapping$gene_name), ]
counts <- counts[gene_mapping_unique$gene_id, ]
rownames(counts) <- gene_mapping_unique$gene_name
## create seruat obj
sce_sub <- CreateSeuratObject(counts, meta.data = meta[, 1:7])
head(rownames(sce_sub), 10)

# 4. preprocess
Idents(sce_sub) <- "orig.ident"
sce_sub[["RNA"]] <- split(sce_sub[["RNA"]], f = sce_sub$orig.ident)
sce_sub <- NormalizeData(sce_sub, verbose = F)
sce_sub <- FindVariableFeatures(sce_sub, verbose = F)
sce_sub <- ScaleData(sce_sub, verbose = F)
sce_sub <- RunPCA(sce_sub, verbose = F)
sce_sub <- IntegrateLayers(sce_sub,
                             method = CCAIntegration,
                             orig.reduction = "pca",
                             new.reduction = "integrated.cca",
                             verbose = F)
sce_sub[["RNA"]] <- JoinLayers(sce_sub[["RNA"]])
sce_sub <- FindNeighbors(sce_sub, reduction = "integrated.cca", dims = 1:30)

# 5. clustering
resolution_list <- seq(0.1, 1, by = 0.1)
for(res in resolution_list){
  print(res)
  sce_sub <- FindClusters(sce_sub, resolution = res, verbose = T)
  sce_sub <- RunUMAP(sce_sub, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
}


# 5. save
saveRDS(sce_sub, "./data/tmp/01.preprocess/sce_sub_clustered.rds")