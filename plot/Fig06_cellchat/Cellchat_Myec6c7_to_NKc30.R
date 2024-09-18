setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(CellChat)
library(ggplot2)

# 2. params
dir_for_obj <- "./data/obj/all_obj.rds"
dir_for_cellchat <- "./data/cellchat/cellchat_Myeloidc3c6_to_NKc5.rds"
dir_for_cellchiat_CASH <- "./data/cellchat/cellchat_Myeloidc3c6_to_NKc5_CASH.rds"

## load data
all_obj <- readRDS(dir_for_obj)
meta <- all_obj@meta.data


# 3. preprocess
table(meta$sub_celltype)
cell.use <-  rownames(meta)[grepl("c06_ Macrophage_|c07_ Macrophage_ MARCO|c30_NK_DNAJB1", meta$sub_celltype)]
sce.use <- all_obj[, cell.use]
# sce.use <- subset(sce.use, group %in% "CASH")

# 4. cellchat
## setup
options(future.globals.maxSize = 1024 * 1024 * 1024 * 54)  # 54 GB
CellchatDB_use <- CellChatDB.human

## processing
run_cellchat <- function(expr_data,
                         label,
                         do_parrallel = TRUE,
                         workers = 10) {
  meta <- data.frame(
    label = label,
    stringsAsFactors = FALSE
  )
  require(CellChat)
  
  rownames(meta) <- colnames(expr_data)
  meta$label <- droplevels(meta$label, exclude = setdiff(levels(meta$label), unique(meta$label)))
  
  cellchat <- createCellChat(object = expr_data, meta = meta, group.by = "label")
  
  # Set the ligand-receptor interaction database
  CellChatDB <- CellChatDB.human
  # use all CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  
  if (do_parrallel) {
    future::plan("multicore", workers = workers) # do parallel
    
    options(future.globals.maxSize = 5 * 1024^3)
  }
  
  print(paste(Sys.time(), "identifyOverExpressedGenes"))
  cellchat <- identifyOverExpressedGenes(cellchat)
  print(paste(Sys.time(), "identifyOverExpressedInteractions"))
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
  # print(paste(Sys.time(), "projectData"))
  # cellchat <- smoothData(cellchat, adj = PPI.human)
  
  # cell-cell communication calculation
  # Compute the communication probability and infer cellular communication network
  print(paste(Sys.time(), "computeCommunProb"))
  cellchat <- computeCommunProb(cellchat, raw.use = T)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  print(paste(Sys.time(), "filterCommunication"))
  cellchat <- filterCommunication(cellchat, min.cells = 5)
  
  # Infer the cell-cell communication at a signaling pathway level
  print(paste(Sys.time(), "computeCommunProbPathway"))
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Calculate the aggregated cell-cell communication network
  print(paste(Sys.time(), "aggregateNet"))
  cellchat <- aggregateNet(cellchat)
  
  # Compute the network centrality scores
  print(paste(Sys.time(), "netAnalysis_computeCentrality"))
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  return(cellchat)
}
data <- list(
  expr_data = sce.use@assays$RNA$data,
  label = sce.use$sub_celltype
  )
cellchat <- run_cellchat(data$expr_data,
                         data$label,
                         do_parrallel = F,
                         workers = 10
                         )

saveRDS(cellchat, "./data/cellchat/cellchat_Myeloidc3c6_to_NKc5.rds")
# saveRDS(cellchat, "./data/cellchat/cellchat_Myeloidc3c6_to_NKc5_CASH.rds")

########## 
cellchat <- readRDS(dir_for_cellchat)
table(cellchat@idents)


# Bubble plot:Cluster1 to Tcell
all_cells <- unique(as.character(sort(cellchat@idents)))
targets_use <- all_cells[!all_cells %in% c("c06_ Macrophage_", "c07_ Macrophage_ MARCO")]

p <- netVisual_bubble(cellchat, 
                      sources.use = c("c06_ Macrophage_", "c07_ Macrophage_ MARCO"),
                      targets.use = targets_use, 
                      remove.isolate = F, 
                      sort.by.target = T) + 
  ggtitle("Cellchat_Myeloidc3c6_to_NKc5")
pdf("./fig/Fig06.cellchat/Cellchat_Myeloidc3c6_to_NKc5.pdf", width = 8, height = 8)
print(p)
dev.off()
