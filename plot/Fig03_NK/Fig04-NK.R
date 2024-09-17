setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(RColorBrewer)
library(ggplot2)


# 2. params
dir_for_obj <- "./data/obj/NK_obj.rds"
Subset_color_panel <- c(
  # NK cells
  "c28_NK_CD160" = "#E18727FF", 
  "c32_NK_CREM" = "#0072B5FF", 
  "c33_NK_IL32" = "#20854EFF",
  "c29_NK_NFKBIA" = "#7876B1FF",
  "c34_NK_CX3CR1" = "black", 
  "c30_NK_DNAJB1" = "red",
  "c31_NK_CCL4" = "#7D4E57"
  )
## load data
NK_obj <- readRDS(dir_for_obj)
NK_obj@meta.data[c("UMAP1", "UMAP2")] <- Embeddings(NK_obj, reduction = "umap.cca")
levels(NK_obj) <- unique(sort(levels(NK_obj)))
NK_obj$sub_celltype <- factor(NK_obj$sub_celltype, levels = unique(sort(levels(NK_obj))))
meta <- NK_obj@meta.data


# 3. markers
markers_list <- list(
  "CD56[bright]CD16[-]NK" = c("NCAM1", "CD160", "NFKBIA", "CCL3", "CCL4", "NR4A3", "DNAJB1", "HSPA1A"),
  "CD56[dim]CD16[+]NK" = c("FCGR3A", "CREM", "IL32", "KLRC2", "CX3CR1")
)


########## 4. UMAP_Dotplot ########## 
dir.create("./fig/Fig04.NK", showWarnings = F)

res <- "RNA_snn_res.0.4"
## UMAP
pdf(paste0("./fig/Fig04.NK/NK_UMAP_", res, ".pdf"), width = 8, height = 8)
p1 <- ggplot(meta, aes(x = UMAP1, y = UMAP2, color = sub_celltype)) +
  geom_point(size = 1, shape = 16, stroke = 0) +
  theme_void() +
  scale_color_manual(values = Subset_color_panel, name = '') +
  theme(aspect.ratio = 1, 
        legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 3)))
print(p1)
dev.off()

## Dotplot
pdf(paste0("./fig/Fig04.NK/NK_Dotplot_", res, ".pdf"), width = 10, height = 8)
p2 <- DotPlot(NK_obj, features = markers_list, dot.scale = 4) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_radius(range = c(0, 4)) +
  scale_color_gradientn(colors = brewer.pal(9, "RdBu")[9:1]) +
  guides(size = guide_legend(title = "Percent expressed"), color = guide_colorbar(title = "Average expression")) +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    text = element_text(size = 8),
    plot.margin = unit(c(1, 1, 1, 1), "char"),
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.grid = element_line(colour = "grey", linetype = "dashed", size = 0.2)
  ) +
  labs(x = "", y = "")
print(p2)
dev.off()


########## 5. Pct_bar_Plot ##########
## Pct bar plot
tmp <- meta[, c("sub_celltype", "orig.ident", "group")]
colnames(tmp) <- c("cluster", "orig.ident", "group")

num_cluster <- tmp %>%
  group_by(cluster, orig.ident) %>%
  summarise(n_clsuter = n(), .groups = 'drop')

num_total <- tmp %>%
  group_by(orig.ident) %>%
  summarise(n_total = n(), .groups = 'drop')

plot_df <- num_cluster %>%
  left_join(num_total, by = "orig.ident") %>%
  mutate(pct = (n_clsuter / n_total) * 100)
plot_df$group <- stringr::str_split(plot_df$orig.ident, "_", simplify = T)[, 1]
plot_df$cluster_short <- stringr::str_split(plot_df$cluster, "_", simplify = T)[,1]

plot <- ggplot(plot_df, aes(x = group, y = pct)) +
  geom_boxplot(aes(color = group), outlier.shape = NA, lwd = 1) +
  geom_jitter(aes(color = group), size = 1, width = 0.2) +
  facet_wrap(~ cluster_short, nrow = 1, labeller = label_wrap_gen(width = 10)) +
  ggsignif::geom_signif(
    comparisons = list(c("CASH", "Normal")),
    test = wilcox.test,
    map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
    step_increase = -0.1,
    textsize = 4,
    size = 1
  ) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c("CASH" = "#5bc0eb", "Normal" = "#9bc53d")) +
  labs(x = "group", y = "pct of cluster(%)", title = "proportion of each subpopulation cell in different samples", fill = "group") +
  theme(
    ## axis
    axis.text = element_text(size = 12, angle = 0, hjust = 1),
    axis.title = element_text(size = 14),
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 1),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_line(size = 1),
    axis.ticks.y = element_line(size = 1),
    ## legend
    legend.position = "bottom",
    legend.justification = "center", 
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    ## title
    plot.title = element_text(size = 14, hjust = 0.5),
    ## facet
    # strip.text.x = element_text(size = 10, angle = 0),
    strip.background = element_rect(color = "white", fill = NA, size = 1)
    ## panel
    # panel.spacing = unit(1, "lines"),
    # panel.border = element_rect(color = "black", fill = NA, size = 1)
  )
pdf("./fig/Fig04.NK/NK_Pct_bar_plot.pdf", width = 10, height = 10)
print(plot)
dev.off()