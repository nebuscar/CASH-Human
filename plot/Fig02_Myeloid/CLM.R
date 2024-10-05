setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(AUCell)

# 2. params
dir_for_expr <- "./data/expr_matrix/GSE136114_matrix_files.txt"
dir_for_gtf <- "./data/gtf.rds"

## load data
expr <- read.csv(dir_for_expr, sep = "\t", row.names = 1)
rownames(expr) <- sub("\\..*", "", rownames(expr))

# 3. convert Ensembl to symbol
gtf <- readRDS(dir_for_gtf)
gene_mapping <- gtf[match(rownames(expr), gtf$gene_id),]
gene_mapping <- na.omit(gene_mapping)
gene_mapping_unique <- gene_mapping[!duplicated(gene_mapping$gene_name), ]
expr <- expr[gene_mapping_unique$gene_id, ]
rownames(expr) <- gene_mapping_unique$gene_name

geneSets <- c("TIMD4", "ITLN1", "MARCO", "CETP", "IFI27",
             "NDST3", "ITGAD", "FABP3", "BCAM", "CD5L", 
             "VCAM1", "C2", "EXT1", "SDC3", "HK3", 
             "ALDH1A1", "PELATON", "SLC7A8", "RND3", "MT1F")

cells_rankings <- AUCell_buildRankings(as.matrix(expr))
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
# cells_AUC2 <- AUCell_run(CD8T_obj@assays$RNA$data, geneSets)
# score <- as.data.frame(t(getAUC(cells_AUC)))
score <- getAUC(cells_AUC)
score <- cells_AUC@assays@data$AUC
