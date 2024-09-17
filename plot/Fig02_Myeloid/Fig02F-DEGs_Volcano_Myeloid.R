setwd("D:/Projects/CASH-Human/")
## 1. library
library(Seurat)
library(AUCell)
library(dplyr)
library(ggplot2)

# 2. params
dir_for_Myeloid_obj <- "./data/obj/Myeloid_obj.rds"
## load data
Myeloid_obj <- readRDS(dir_for_Myeloid_obj)
## subset
tmp_obj = subset(Myeloid_obj, sub_celltype %in% c("c06_Macrophage_", "c07_Macrophage_MARCO"))
levels(tmp_obj) <- unique(sort(levels(Myeloid_obj)))

# 3. Findmarkers
DEGs <- FindMarkers(tmp_obj, 
                    ident.1 = "c07_Macrophage_MARCO", 
                    ident.2 = "c06_Macrophage_", 
                    logfc.threshold = 0.1, 
                    min.pct = 0.01, 
                    only.pos = F)
DEGs <- DEGs %>% 
  arrange(desc(avg_log2FC))
write.csv(DEGs, "./data/tmp/02.DEGs/DEGs_Myeloid-c07_to_c06.csv")

# 4. Volcano
dir_for_DEGs <- "./data/tmp/02.DEGs/DEGs_Myeloid-c07_to_c06.csv"
DEGs <- read.csv(dir_for_DEGs, row.names = 1)
## setup threshold
log2FC = 1
padj = 0.05
DEGs$group <- ifelse(DEGs$p_val_adj < padj & abs(DEGs$avg_log2FC) >= log2FC,
                     ifelse(DEGs$avg_log2FC > log2FC,"c07_Macrophage_MARCO", "c06_Macrophage_"), "Ns")
DEGs$group <- factor(DEGs$group)
print(table(DEGs$group))
plot_df <- DEGs

p1 <- ggplot(plot_df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = group)) +
  geom_point(alpha = 0.8, size = 0.8) +
  geom_vline(xintercept = c(-log2FC, log2FC), linetype = 2, color = "grey") +
  labs(x = bquote(Log[2] * FoldChange), y = bquote(-Log[10] * italic(P.adj))) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c("c07_Macrophage_MARCO" = "red", 
                                "Ns" = "grey", 
                                "c06_Macrophage_" = "blue")) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.75))) +
  theme(
    ## axis
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12), 
    axis.text = element_text(size = 10),
    
    ## legend
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.position = "bottom"
  ) +
  scale_x_continuous(limits = c(-max(abs(plot_df$avg_log2FC)), max(abs(plot_df$avg_log2FC))))

## significant markers
sig_markers <- subset(plot_df, plot_df$p_val_adj < padj & abs(plot_df$avg_log2FC) >= log2FC)
top10_c06 <- sig_markers %>%
  filter(group == "c07_Macrophage_MARCO") %>%
  arrange(desc(avg_log2FC)) %>%
  head(10)
top10_c03 <- sig_markers %>%
  filter(group == "c06_Macrophage_") %>%
  arrange(avg_log2FC) %>%
  head(10)
top20 <- rbind(top10_c06, top10_c03)

## Add labels to the plot
p2 <- p1 + ggrepel::geom_text_repel(data = top20, aes(label = rownames(top20)), 
                                    size = 3,                                    
                                    color = "black",
                                    fontface = "bold.italic",
                                    box.padding = unit(0.5, "lines"),
                                    point.padding = unit(1, "lines"), 
                                    segment.color = "black", 
                                    show.legend = F,
                                    max.overlaps = 20
                                    )
print(p2)
ggsave("./fig/Fig02.Myeloid/Fig2F.DEGs_Volcano_Myeloid.pdf", plot = p2, width = 6, height = 6)
ggsave("./fig/Fig02.Myeloid/Fig2F.DEGs_Volcano_Myeloid.png", plot = p2, width = 6, height = 6)
