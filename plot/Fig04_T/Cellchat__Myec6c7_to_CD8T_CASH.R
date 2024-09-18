setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(CellChat)
library(ggplot2)

# 2. params
dir_for_all_obj <- "./data/obj/All_obj.rds"

dir_for_T_obj <- "./data/obj/T_obj.rds"
dir_for_meta_tcr <- "./data/tcr/meta_tcr.rds"

## load data
### obj
all_obj <- readRDS(dir_for_all_obj)
meta <- all_obj@meta.data
### meta_tcr
meta_tcr <- readRDS(dir_for_meta_tcr)


# 3. preprocess
table(meta$sub_celltype)
meta_sub <- meta %>%
  filter((major_celltype == "T" & str_detect(sub_celltype, "CD8")) | 
           sub_celltype %in% c("c06_ Macrophage_", "c07_ Macrophage_ MARCO"))
## match meta
meta_match <- meta_tcr[match(rownames(meta_sub), rownames(meta_tcr)),]
meta_sub$cloneSize <- meta_match$cloneSize

meta_sub <- meta_sub %>%
  mutate(label = case_when(
    major_celltype == "T" & cloneSize %in% c("Hyperexpanded (0.1 < X <= 1)", "Large (0.01 < X <= 0.1)", "Medium (0.001 < X <= 0.01)") ~ "CD8T_highClone",
    major_celltype == "T" & cloneSize %in% c("Small (1e-04 < X <= 0.001)", "Rare (0 < X <= 1e-04)", "None ( < X <= 0)") ~ "CD8T_lowClone",
    major_celltype != "T" ~ as.character(sub_celltype)
    ))
meta_sub <- meta_sub %>%
  filter(!is.na(label))
meta_sub$sub_celltype <- factor(meta_sub$sub_celltype)
meta_sub$label <- factor(meta_sub$label)
cell.use <-  rownames(meta_sub)
sce.use <- all_obj[, cell.use]
sce.use@meta.data <- meta_sub
sce.use <- subset(sce.use, group %in% "CASH")
sce.use$label <- factor(sce.use$label)

# 4. cellchat
## setup
options(future.globals.maxSize = 1024 * 1024 * 1024 * 54)
CellchatDB_use <- CellChatDB.human

## processing
run_cellchat <- function(expr_data, 
                         label, 
                         do_parrallel = TRUE, 
                         workers = 10) {
  meta <- data.frame(
    label = label,
    stringsAsFactors = FALSE)
  require(CellChat)
  
  rownames(meta) <- colnames(expr_data)
  # meta$labels <- droplevels(meta$labels, exclude = setdiff(levels(meta$labels), unique(meta$labels)))
  
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
  label = sce.use$label)
cellchat <- run_cellchat(data$expr_data,
                         data$label,
                         do_parrallel = F)

# saveRDS(cellchat, "./data/cellchat/cellchat_Myeloidc6c7_to_allT.rds")
saveRDS(cellchat, "./data/cellchat/cellchat_Myeloid_c6c7_to_CD8T_CASH.rds")

########## 
dir_for_cellchat_CASH <- "./data/cellchat/cellchat_Myeloid_c6c7_to_CD8t_CASH.rds"
cellchat <- readRDS(dir_for_cellchat_CASH)
table(cellchat@idents)


# Bubble plot:Cluster1 to Tcell
all_cells <- unique(as.character(sort(cellchat@idents)))
targets_use <- all_cells[!all_cells %in% c("c06_ Macrophage_", "c07_ Macrophage_ MARCO")]

p <- netVisual_bubble(cellchat, 
                      sources.use = c("c06_ Macrophage_", "c07_ Macrophage_ MARCO"),
                      targets.use = targets_use, 
                      remove.isolate = F, 
                      sort.by.target = T) + 
  labs(title = "Cellchat_Myeloid_c6c7_to_CD8T_CASH") +
  theme(
    ## text
    axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14),
    
    ## legend
    legend.text = element_text(size = 10),
    legend.position = "right") 
print(p)
ggsave("./fig/Fig04.T/Cellchat_Myeloid_c6c7_to_CD8T_CASH.pdf", plot = p, width = 10, height = 10, dpi = 300)
ggsave("./fig/Fig04.T/Cellchat_Myeloid_c6c7_to_CD8T_CASH.png", plot = p, width = 10, height = 10, dpi = 300)

