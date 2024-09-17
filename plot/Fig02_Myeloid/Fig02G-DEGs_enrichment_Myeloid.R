setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(irGSEA)
library(dplyr)
library(ggplot2)
library(stringr)

# 2. params
dir_for_Myeloid_obj <- "./data/obj/Myeloid_obj.rds"

## load data
Myeloid_obj <- readRDS(dir_for_Myeloid_obj)
Idents(Myeloid_obj) <- "RNA_snn_res.0.7"
tmp_obj = subset(Myeloid_obj, RNA_snn_res.0.7 %in% c("3", "6"))


########## 3. DEGs-Pathway_AUCell ##########
tmp_obj <- irGSEA.score(object = tmp_obj, 
                        assay = "RNA", 
                        slot = "data", seeds = 42, ncores = 10,
                        custom = F, geneset = NULL, msigdb = T, 
                        species = "Homo sapiens", category = "C2", subcategory = NULL, geneid = "symbol",
                        method = "AUCell")
saveRDS(tmp_obj, "./data/tmp/02.DEGs/Myeloid_obj_score.rds")

dir_for_tmp_obj <- "./data/tmp/02.DEGs/Myeloid_obj_score.rds"
tmp_obj <- readRDS(dir_for_tmp_obj)
DEGs <- irGSEA.integrate(object = tmp_obj,
                         group.by = "sub_celltype",
                         method = "AUCell")[["AUCell"]]

saveRDS(DEGs, "./data/DEGs/DEGs_pathway_Myeloid.rds")


# Barplot
dir_for_DEGs_pathway <- "./data/DEGs/DEGs_pathway_Myeloid.rds"
DEGs <-  readRDS(dir_for_DEGs_pathway)
DEGs_to_use <- c(
  "REACTOME-ALTERNATIVE-COMPLEMENT-ACTIVATION",
  "REACTOME-HDL-REMODELING",
  "WP-LIPID-PARTICLES-COMPOSITION",
  "WP-IL10-ANTIINFLAMMATORY-SIGNALING-PATHWAY",
  "REACTOME-CHYLOMICRON-CLEARANCE",
  "REACTOME-METALLOTHIONEINS-BIND-METALS",
  "REACTOME-LDL-REMODELING",
  "REACTOME-RESPONSE-TO-METAL-IONS",
  "REACTOME-PLASMA-LIPOPROTEIN-ASSEMBLY",
  "REACTOME-CHYLOMICRON-REMODELING",
  "LI-ADIPOGENESIS-BY-ACTIVATED-PPARG",
  "MYLLYKANGAS-AMPLIFICATION-HOT-SPOT-30",
  "WP-CYTOPLASMIC-RIBOSOMAL-PROTEINS",
  "TIAN-TNF-SIGNALING-NOT-VIA-NFKB",
  "NAGASHIMA-EGF-SIGNALING-UP",
  "YAMASHITA-LIVER-CANCER-WITH-EPCAM-UP",
  "TUOMISTO-TUMOR-SUPPRESSION-BY-COL13A1-UP")
top <- DEGs %>%
  filter(DEGs$Name %in% DEGs_to_use) %>% 
  filter(direction == "up") %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(avg_diff)) %>% 
  # slice_max(order_by = avg_diff, n = 10) %>%
  ungroup()

plot_df <- top %>%
  mutate(log10p = -log10(abs(p_val_adj)))
plot_df <- plot_df %>% 
  mutate(log10p = ifelse(cluster == "c06_Macrophage_", -log10p, log10p))

#ggplot
p <- ggplot(plot_df, aes(reorder(Name, log10p), log10p, fill = cluster))+
  geom_col() +
  theme_bw() +
  scale_fill_manual(name = "", values = c("c06_Macrophage_" = "#1084A4", "c07_Macrophage_MARCO" = "#8D4873")) + 
  ## label
  labs(x = '', y = bquote(-Log[10] * italic(P.adj)), 
       title = "TOP Pathways between c06 and c07 cells") +
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
ggsave("./fig/Fig02.Myeloid/Fig2G.Top10_DEGs_pathway_Myeloid.pdf", plot = p, width = 5, height = 4, dpi = 300)
ggsave("./fig/Fig02.Myeloid/Fig2G.Top10_DEGs_pathway_Myeloid.png", plot = p, width = 5, height = 4, dpi = 300)


