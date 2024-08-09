# Figs for the paper
# 2024-08-08
library(patchwork)

tiff(paste("figs/avoidance_plot.tiff",sep = ""),
     height = 14, width = 14, units = "cm", res = 300, compression = "lzw",
     pointsize = 8)
p_avoidance_duiker + p_avoidance_impala +
  plot_annotation(tag_levels = c('a'), tag_prefix = '(',
                  tag_suffix = ')')
dev.off()
  
tiff(paste("figs/detection_plot.tiff",sep = ""),
     height = 14, width = 14, units = "cm", res = 300, compression = "lzw",
     pointsize = 8)
p_detection_duiker / p_detection_impala +
  plot_annotation(tag_levels = c('a'), tag_prefix = '(',
                  tag_suffix = ')')
dev.off()

tiff(paste("figs/cds_plot.tiff",sep = ""),
     height = 14, width = 14, units = "cm", res = 300, compression = "lzw",
     pointsize = 8)
p_cds_duiker / p_cds_impala +
  plot_annotation(tag_levels = c('a'), tag_prefix = '(',
                  tag_suffix = ')')
dev.off()
