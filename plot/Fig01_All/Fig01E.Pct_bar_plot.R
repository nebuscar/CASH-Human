setwd("D:/Projects/CASH-Human/")
# 1. library
library(Seurat)
library(dplyr)
library(ggplot2)

# 2. params
dir_for_all_obj <- "./data/obj/All_obj.rds"
celltype_color_panel <- c(
  "T cells" = "#98c9dd",       
  "NK" = "#a6d38e",           
  "Monocytes" = "#37a849",     
  "Macrophage" = "#207cb5", 
  "B cells" = "#f69595",      
  "Plasma cells" = "#eb2a2a", 
  "pDC" = "#9467bd",         
  "cDC1" = "#fcba71",         
  "cDC2" = "#f78200" 
)
## load data
all_obj <- readRDS(dir_for_all_obj)
meta <- all_obj@meta.data

########## Fig01E.Pct_bar_plot ##########
# 3. plot
dir.create("./fig/Fig01.All", showWarnings = F)

tmp <- meta[, c("orig.ident", "celltype")]

num_celltype <- tmp %>% 
  group_by(celltype, orig.ident) %>% 
  summarise(n_cluster = n(), .groups = "drop")
num_total <- tmp %>% 
  group_by(orig.ident)%>% 
  summarise(n_total = n(), .groups = "drop")
plot_df <- num_celltype %>%
  left_join(num_total, by = "orig.ident") %>%
  mutate(pct = (n_cluster / n_total) * 100)
plot_df$orig.ident <- factor(plot_df$orig.ident, levels = c("CASH_Patient1", "CASH_Patient2", "CASH_Patient3", "CASH_Patient4", "CASH_Patient5",
                                                            "Normal_Patient1", "Normal_Patient2", "Normal_Patient3", "Normal_Patient4", "Normal_Patient5",
                                                            "Normal_Patient6", "Normal_Patient7", "Normal_Patient8", "Normal_Patient9", "Normal_Patient10"))

p <- ggplot(plot_df, aes(x = orig.ident, y = pct, fill =  celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  cowplot::theme_cowplot() +
  scale_fill_manual(values = celltype_color_panel) + 
  labs(x = "", y = "Proportion", title = "") +
  theme(
    ## axis
    axis.text.x = element_text(size = 14, angle = 90),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 18, face = "bold"),
    
    ## legend
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.position = "right",
    
    ## panel
    axis.line.x = element_line(linewidth = 1),
    axis.line.y = element_line(linewidth = 1)) +
  coord_flip()
print(p)
ggsave("./fig/Fig01.All/Fig1E.Pct_bar_plot.pdf", plot = p, width = 8, height = 8, dpi = 300)


    
