setwd("D:/Projects/CASH-Human/")
# 1. library
library(dplyr)
library(tidyr)
library(CellChat)
library(Seurat)
library(ggplot2)

# 2. params
dir_for_cellchat_Myeloid_T <- "./data/cellchat/"
dir_for_cellchat_Myeloid_T <- "./data/cellchat/"
dir_for_cellchat_Myeloid_B <- "./data/cellchat/cellchat_Myeloid_B.rds"
dir_for_cellchat_Myeloid_T <- "./data/cellchat/cellchat_Myeloid_T.rds"

## load data
file_list <- list.files(dir_for_data)

########## Cellchat-Myeloid Vs B ##########
cellchat <- readRDS(dir_for_cellchat_Myeloid_B)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# Bubble plot:Myeloid to Tcell
all_cells <- unique(as.character(sort(cellchat@idents)))
targets_use <- all_cells[!all_cells %in% c("Myeloid_c3", "Myeloid_c6")]

p <- netVisual_bubble(cellchat, sources.use = c("Myeloid_c3", "Myeloid_c6"), 
                      targets.use = targets_use, 
                      remove.isolate = F) + 
  ggtitle("Cellchat_Myeloid_to_B")
pdf("../fig/Cellchat_Myeloid_to_B.pdf", width = 10, height = 10)
print(p)
dev.off()
