yy_Dotplot <- function(seuratObj,
                       genes, # list (grouped) or vector (ungrouped)
                       group.by,
                       coord_flip = FALSE,
                       scale = TRUE,
                       dot.scale = 4,
                       gene_expr_cutoff = -Inf,
                       gene_pct_cutoff = 0,
                       cell_expr_cutoff = -Inf,
                       cell_pct_cutoff = 0,
                       panel.spacing_distance = 0.5,
                       return_data = FALSE) {
  # Required libraries
  require(dplyr)
  require(Seurat)
  require(RColorBrewer)
  require(ggplot2)
  require(cowplot)
  require(tidyverse)
  
  # Parameters
  col.min <- -2.5
  col.max <- 2.5
  dot.min <- 0
  
  # Reformat genes
  if (is.list(genes)) {
    if (is.null(names(genes))) {
      names(genes) <- 1:length(genes)
    }
    genes <- stack(genes)
    features <- genes$values
    features_group <- genes$ind
  } else {
    features <- genes
    features_group <- NULL
  }
  
  # Check for missing genes
  subset <- features %in% rownames(seuratObj)
  if (sum(!subset) > 0) {
    cat(paste0(paste(features[!subset], collapse = ", "), " is missing in gene list"))
  }
  if (length(features) == 0) {
    stop("No intersecting genes, please check gene name format.\n")
  }
  
  # Prepare plot input
  data.features <- FetchData(seuratObj, cells = colnames(seuratObj), vars = features, slot = "data")
  data.features$id <- seuratObj@meta.data[[group.by]]
  
  data.plot <- lapply(unique(data.features$id), function(ident) {
    data.use <- data.features[data.features$id == ident, 1:(ncol(data.features) - 1), drop = FALSE]
    avg.exp <- apply(data.use, 2, function(x) mean(expm1(x)))
    pct.exp <- apply(data.use, 2, function(x) sum(x > 0) / length(x))
    list(avg.exp = avg.exp, pct.exp = pct.exp)
  })
  names(data.plot) <- unique(data.features$id)
  data.plot <- lapply(names(data.plot), function(x) {
    data.use <- as.data.frame(data.plot[[x]])
    data.use$features.plot <- rownames(data.use)
    data.use$features.plot_show <- features
    data.use$id <- x
    data.use
  })
  data.plot <- do.call(rbind, data.plot)
  
  avg.exp.scaled <- sapply(unique(data.plot$features.plot), function(x) {
    data.use <- data.plot[data.plot$features.plot == x, "avg.exp"]
    if (scale) {
      data.use <- scale(data.use)
      MinMax(data.use, min = col.min, max = col.max)
    } else {
      log1p(data.use)
    }
  })
  avg.exp.scaled <- as.vector(t(avg.exp.scaled))
  
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(data.plot$features.plot, levels = unique(data.plot$features.plot))
  data.plot$features.plot_show <- factor(data.plot$features.plot_show, levels = unique(features))
  
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  
  if (is.factor(seuratObj@meta.data[[group.by]])) {
    data.plot$id <- factor(data.plot$id, levels = levels(seuratObj@meta.data[[group.by]]))
  }
  
  if (!is.null(features_group)) {
    data.plot <- data.plot %>%
      left_join(data.frame(features = features, features_group = features_group), by = c("features.plot_show" = "features")) %>%
      mutate(features_group = factor(features_group, levels = unique(features_group)))
  }
  
  # Filter genes/cells
  filter_cell <- data.plot %>%
    group_by(id) %>%
    summarise(max_exp = max(avg.exp.scaled), max_pct = max(pct.exp)) %>%
    filter(max_exp >= cell_expr_cutoff & max_pct >= cell_pct_cutoff)
  filter_gene <- data.plot %>%
    group_by(features.plot) %>%
    summarise(max_exp = max(avg.exp.scaled), max_pct = max(pct.exp)) %>%
    filter(max_exp >= gene_expr_cutoff & max_pct >= gene_pct_cutoff)
  
  data.plot <- data.plot %>%
    filter(features.plot %in% filter_gene$features.plot) %>%
    filter(id %in% filter_cell$id)
  
  # Plot
  plot <- ggplot(data.plot, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    scale_radius(range = c(0, dot.scale)) +
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
  
  if (coord_flip) {
    plot <- plot + coord_flip()
  }
  
  if (!is.null(features_group)) {
    plot <- plot +
      facet_grid(. ~ features_group, scales = "free_x", space = "free_x", switch = "y") +
      theme(panel.spacing = unit(panel.spacing_distance, "lines"), strip.background = element_blank())
  }
  
  if (return_data) {
    return(list(plot = plot, data = data.plot))
  } else {
    return(plot)
  }
}