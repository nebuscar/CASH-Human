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
# meta_tcr <- sce_Tcell_TCR@meta.data[, c("clonalProportion", "clonalFrequency", "cloneSize")]
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

T_obj@meta.data <- T_obj@meta.data %>%
  mutate(cloneGroup = case_when(
    cloneSize %in% c("Hyperexpanded (0.1 < X <= 1)", "Large (0.01 < X <= 0.1)", "Medium (0.001 < X <= 0.01)") ~ "highClone",
    cloneSize %in% c("Small (1e-04 < X <= 0.001)", "Rare (0 < X <= 1e-04)", "None ( < X <= 0)") ~ "lowClone"
    ))
tmp_obj <- subset(T_obj, cloneGroup %in% c("highClone", "lowClone"))


########## 3. DEGs-Pathway_AUCell ##########
tmp_obj <- irGSEA.score(object = tmp_obj, 
                      assay = "RNA", 
                      slot = "data", seeds = 42, ncores = 10,
                      custom = F, geneset = NULL, msigdb = T, 
                      species = "Homo sapiens", category = "H", subcategory = NULL, geneid = "symbol",
                      method = "AUCell")
DEGs <- irGSEA.integrate(object = tmp_obj,
                         group.by = "cloneGroup",
                         method = "AUCell")[["AUCell"]]
write.csv(DEGs, "./data/DEGs/DEGs_pathway_T_cloneGroup.csv")

saveRDS(DEGs, "./data/DEGs/DEGs_pathway_T_cloneGroup.rds")
# Barplot
dir_for_DEGs_pathway <- "./data/DEGs/DEGs_pathway_T.rds"
DEGs <-  readRDS(dir_for_DEGs_pathway)
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
  mutate(log10p = ifelse(cluster == "highClone", -log10p, log10p))

#ggplot
p <- ggplot(plot_df, aes(reorder(Name, log10p), log10p, fill = cluster))+
  geom_col() +
  theme_bw() +
  scale_fill_manual(name = "", values = c("highClone" = "#1084A4", "lowClone" = "#8D4873")) + 
  ## label
  labs(x = '', y = bquote(-Log[10] * italic(P.adj)), 
       title = "TOP Pathways between highClone and lowClone group across T cells") +
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
  ylim(-300, 300)
print(p)
ggsave("./fig/Fig04.T/Fig4E.Top10_DEGs_pathway_T.pdf", plot = p, width = 5, height = 4, dpi = 300)
ggsave("./fig/Fig04.T/Fig4E.Top10_DEGs_pathway_T.png", plot = p, width = 5, height = 4, dpi = 300)
