setwd("D:/Projects/CASH-Human/")
# 1. library
library(limma)
library(sva)
library(AUCell)

# 2. params
dir_for_counts <- "./data/expr_matrix/counts/"
dir_for_gtf <- "./data/gtf.rds"
## load data
file.list <- list.files(dir_for_counts)
counts.list <- list()
for (i_file in file.list) {
  print(i_file)
  counts.list[[i_file]] <- read.csv(file.path(dir_for_counts, i_file), sep = "\t", row.names = 1)
  ## convert Ensembl to symbol
  gtf <- readRDS(dir_for_gtf)
  gene_mapping <- gtf[match(rownames(counts.list[[i_file]]), gtf$gene_id),]
  gene_mapping <- na.omit(gene_mapping)
  gene_mapping_unique <- gene_mapping[!duplicated(gene_mapping$gene_name), ]
  counts.list[[i_file]] <- counts.list[[i_file]][gene_mapping_unique$gene_id, ]
  rownames(counts.list[[i_file]]) <- gene_mapping_unique$gene_name
  rm(gtf, gene_mapping, gene_mapping_unique)
  
  counts.list[[i_file]] <- normalizeBetweenArrays(counts.list[[i_file]])
}

# 3. preprocess
gene <- intersect(rownames(counts.list[["HCC_PD1_1.txt"]]), rownames(counts.list[["HCC_PD1_2.txt"]]))
expr_merged <- cbind(counts.list[["HCC_PD1_1.txt"]][gene, ], counts.list[["HCC_PD1_2.txt"]][gene, ])
boxplot(expr_merged, outline = F, notch = T, las = 2)

batch <- c(rep(1, ncol(counts.list[["HCC_PD1_1.txt"]])), 
           rep(2, ncol(counts.list[["HCC_PD1_2.txt"]])))
expr_merged <- ComBat(dat = expr_merged, batch = batch)
boxplot(expr_merged, outline = F, notch = T, las = 2)

# 4. group by AUCell_score
geneSets <- c("TIMD4", "ITLN1", "MARCO", "CETP", "IFI27",
              "NDST3", "ITGAD", "FABP3", "BCAM", "CD5L", 
              "VCAM1", "C2", "EXT1", "SDC3", "HK3", 
              "ALDH1A1", "PELATON", "SLC7A8", "RND3", "MT1F")
cells_rankings <- AUCell_buildRankings(as.matrix(expr_merged), plotStats = TRUE)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
auc_scores <- getAUC(cells_AUC)

median_score <- median(auc_scores, na.rm = TRUE)
group <- factor(ifelse(auc_scores > median_score, "high", "low"), 
                   levels = c("low", "high"))

#
expr_merged <- data.frame(t(expr_merged), scoreGroup = group)
expr_merged <- expr_merged[, c(ncol(expr_merged), 1:(ncol(expr_merged) - 1))]
