setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(irGSEA)
library(dplyr)
library(ggplot2)
library(stringr)


# 2. params
# dir_for_tcr <- "../CASH/tmp/sce_Tcell_TCR.rds"
# sce_Tcell_TCR <- readRDS(dir_for_tcr)
# meta_tcr <- sce_Tcell_TCR@meta.data
# saveRDS(meta_tcr, "./data/tcr/meta_tcr.rds")

dir_for_T_obj <- "./data/obj/T_obj.rds"
dir_for_meta_tcr <- "./data/tcr/meta_tcr.rds"

## loda data
T_obj <- readRDS(dir_for_T_obj)
Idents(T_obj) <- "sub_celltype"
meta_tcr <- readRDS(dir_for_meta_tcr)
## match meta
meta_match <- meta_tcr[match(rownames(T_obj@meta.data), rownames(meta_tcr)),]
T_obj@meta.data$cloneSize <- meta_match$cloneSize

##########
meta_match$sub_celltype <- T_obj@meta.data$sub_celltype
# meta_match$sequence <- stringr::str_split(meta_match$CTgene, "\\.", simplify = T)[,1]
tmp <- tidyr::separate(meta_match, CTgene, into = c("sequence", "rest"), sep = "\\.", extra = "merge", fill = "right")
t <- subset(tmp, sub_celltype %in% c("c17_CD8_", "c18_CD8_"))
t$sub_celltype <- factor(t$sub_celltype)

df <- data.frame(table(t$sequence, t$sub_celltype))
df_trav <- subset(df, grepl("TRAV", Var1))

num_total <- df_trav %>%
  summarise(n_total = sum(Freq))

df_trav <- df_trav %>%
  mutate("pct(%)" = Freq / num_total$n_total * 100)
write.csv(df_trav, "./pct_tcr.csv")
##########
meta_match$sub_celltype <- T_obj@meta.data$sub_celltype
# meta_match$sequence <- stringr::str_split(meta_match$CTgene, "\\.", simplify = T)[,1]
tmp <- tidyr::separate(meta_match, CTgene, into = c("sequence", "rest"), sep = "\\.", extra = "merge", fill = "right")
t <- subset(tmp, sub_celltype %in% c("c14_CD8_MAIT", "c15_CD8_MAIT", "c16_CD8_MAIT"))
t$sub_celltype <- factor(t$sub_celltype)

df <- data.frame(table(t$sequence, t$sub_celltype))
df_trav <- subset(df, grepl("TRAV", Var1))

num_total <- df_trav %>%
  summarise(n_total = sum(Freq))

df_trav <- df_trav %>%
  mutate("pct(%)" = round(Freq / num_total$n_total * 100, 4))
df_trav <- df_trav %>%
  mutate("pct(%)" = Freq / num_total$n_total * 100)
write.csv(df_trav, "./pct_tcrAV_MAIT.csv", row.names = F)
##########

T_obj@meta.data <- T_obj@meta.data %>%
  mutate(cloneGroup = case_when(
    cloneSize %in% c("Hyperexpanded (0.1 < X <= 1)", "Large (0.01 < X <= 0.1)", "Medium (0.001 < X <= 0.01)") ~ "highClone",
    cloneSize %in% c("Small (1e-04 < X <= 0.001)", "Rare (0 < X <= 1e-04)", "None ( < X <= 0)") ~ "lowClone"
  ))
tmp_obj <- subset(T_obj, cloneGroup %in% c("highClone", "lowClone"))
tmp_obj <- subset(tmp_obj, group %in% "CASH")


cluster_to_use <- c("c17_CD8_", "c18_CD8_")
DEGs <- list()
for (i_cluster in cluster_to_use) {
  print(i_cluster)
  tmp_to_use <- subset(tmp_obj, sub_celltype %in% cluster_to_use)
  tmp_to_use$sub_celltype <- factor(tmp_to_use$sub_celltype)
  print(table(tmp_to_use$sub_celltype,  tmp_to_use$cloneGroup))
  ########## 3. DEGs-Pathway_AUCell ##########
  tmp_to_use <- irGSEA.score(object = tmp_to_use, 
                          assay = "RNA", 
                          slot = "data", seeds = 42, ncores = 10,
                          custom = F, geneset = NULL, msigdb = T, 
                          species = "Homo sapiens", category = "C2", subcategory = NULL, geneid = "symbol",
                          method = "AUCell")
  DEGs <- irGSEA.integrate(object = tmp_to_use,
                           group.by = "cloneGroup",
                           method = "AUCell")[["AUCell"]]
  write.csv(DEGs, "./data/DEGs/DEGs_pathway_c17c18T_cloneGroup-CASH.csv")
} 
for (i_cluster in cluster_to_use) {
  print(i_cluster)
  write.csv(DEGs[[i_cluster]], paste0("./data/DEGs/DEGs_pathway_", i_cluster, "_cloneGroup-CASH.csv"))
  saveRDS(DEGs[[i_cluster]], paste0("./data/DEGs/DEGs_pathway_", i_cluster, "_cloneGroup-CASH.rds"))
}

# Barplot
file_list <- c("./data/DEGs/DEGs_pathway_c17_CD8__cloneGroup-CASH.rds", 
                          "./data/DEGs/DEGs_pathway_c18_CD8__cloneGroup-CASH.rds")

for (i_cluster in file_list) {
  print(stringr::str_split(i_cluster, "_", simplify = T)[, 3])
  i_cluster <- stringr::str_split(i_cluster, "_", simplify = T)[, 3]
  DEGs <-  readRDS(i_cluster)
  top <- DEGs %>%
    filter(DEGs$cluster %in% c("highClone", "lowClone")) %>%
    filter(direction == "up") %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    arrange(desc(avg_diff)) %>% 
    slice_max(order_by = avg_diff, n = 10) %>%
    ungroup()
  
  plot_df <- top %>%
    mutate(log10p = -log10(abs(p_val_adj)))
  plot_df <- plot_df %>% 
    mutate(log10p = ifelse(cluster == "lowClone", -log10p, log10p))
  
  #ggplot
  p <- ggplot(plot_df, aes(reorder(Name, log10p), log10p, fill = cluster))+
    geom_col() +
    theme_bw() +
    scale_fill_manual(name = "", values = c("lowClone" = "#1084A4", "highClone" = "#8D4873")) + 
    ## label
    labs(x = '', y = bquote(-Log[10] * italic(P.adj)), 
         title = paste0("TOP Pathways between highClone and lowClone group across c17c18 T cells-CASH")) +
    geom_text(data = plot_df[which(plot_df$log10p > 0),],aes(x = Name, y = -0.5, label = Name),
              hjust = 1, size = 2) +
    geom_text(data = plot_df[which(plot_df$log10p < 0),],aes(x = Name, y = 0.5, label = Name),
              hjust = 0, size = 2) +
    theme(
      ## axis
      text = element_text(size = 7),
      axis.text.x = element_text(color = "black", size = 7, angle = 0),
      axis.text.y = element_blank(),
      
      axis.line.x = element_line(color = "black", linewidth = 0.3),
      axis.ticks.x = element_line(linewidth = 0.3),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      
      
      ## panel
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      panel.border = element_blank(),
      
      ## legend
      legend.title = element_blank(),
      legend.text = element_text(size = 7),
      legend.justification = "center",
      legend.position = "bottom")+
    coord_flip() +
    geom_segment(aes(y = 0, yend = 0, x = 0)) +
    scale_x_discrete(expand = expansion(mult = c(0,0))) +
    ylim(-50, 50)
  print(p)
  }
ggsave("./fig/Fig04.T/Fig4E.Top10_DEGs_pathway_T.pdf", plot = p, width = 5, height = 4, dpi = 300)
ggsave("./fig/Fig04.T/Fig4E.Top10_DEGs_pathway_T.png", plot = p, width = 5, height = 4, dpi = 300)
