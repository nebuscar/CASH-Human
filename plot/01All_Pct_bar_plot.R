setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(dplyr)
library(ggplot2)

# 2. params
dir_for_obj <- "./data/obj/"

## load data
file_list <- list.files(dir_for_obj)
file_list <- file_list[!file_list %in% c("CASH-Normal-harmony.rds", "All_obj.rds")]

########## 5. Pct_bar_Plot ##########
for (i_file in file_list) {
  cell <- stringr::str_split(i_file, "_", simplify = T)[,1]
  sce_sub <- readRDS(file.path(dir_for_obj, i_file))
  meta <- sce_sub@meta.data
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
  
  p <- ggplot(plot_df, aes(x = group, y = pct)) +
    geom_boxplot(aes(color = group), outlier.shape = NA, lwd = 1) +
    geom_jitter(aes(color = group), size = 1, width = 0.2) +
    facet_wrap(~ cluster_short, nrow = 1, labeller = label_wrap_gen(width = 10)) +
    ggsignif::geom_signif(
      comparisons = list(c("CASH", "Normal")),
      test = wilcox.test,
      map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
      step_increase = -0.1,
      textsize = 6,
      size = 1
    ) +
    cowplot::theme_cowplot() +
    scale_color_manual(values = c("CASH" = "#5bc0eb", "Normal" = "#9bc53d")) +
    labs(x = "group", y = "pct of clusters(%)", title = "", fill = "group") +
    theme(
      ## axis
      axis.text = element_text(size = 12, angle = 0, hjust = 1),
      axis.title = element_text(size = 18, face = "bold"),
      axis.line = element_line(size = 1),
      axis.ticks = element_line(size = 1),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      
      ## legend
      legend.position = "right",
      legend.justification = "center", 
      legend.text = element_text(size = 18),
      legend.title = element_blank(),
      
      ## title
      plot.title = element_text(size = 14, hjust = 0.5),
      ## facet
      strip.text.x = element_text(size = 14, face = "bold"),
      strip.background = element_rect(color = NA, fill = NA, size = 1)
      ## panel
      # panel.spacing = unit(1, "lines"),
      # panel.border = element_rect(color = "black", fill = NA, size = 1)
    )
  print(p)
  ggsave(paste0("./fig/01Pct_bar_plot/", cell, "_Pct_bar_plot.png"), plot = p, width = 16, height = 8, dpi = 300)
}