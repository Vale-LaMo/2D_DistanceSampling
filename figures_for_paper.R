# Figs for the paper
# 2024-08-08
library(patchwork)


####---- Figure avoidance, detection e cds:  ----
# per creare le seguenti 3 figure,
# è necessario prima far girare il file 00_LT2D_[ultimaversione].Rmd
# sia per l'impala sia per il duiker

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


####---- Figura simulazioni  ----
# per creare la seguente figura, sono stati usati file creati con
# simulazioni.R (far rigirare il file)
tiff(paste("figs/simula_plot.tiff",sep = ""),
     height = 14, width = 14, units = "cm", res = 300, compression = "lzw",
     pointsize = 8)
p_simula_duiker / p_simula_impala +
  plot_annotation(tag_levels = c('a'), tag_prefix = '(',
                  tag_suffix = ')')
dev.off()


####---- Figura mappa DSM  ----
# per creare la seguente figura, sono stati usati file creati con
# 02_dsm_[ultimaversione].R (andare direttamente all'ultimo chunk)
tiff(paste("figs/dsm_plot.tiff",sep = ""),
     height = 14, width = 14, units = "cm", res = 300, compression = "lzw",
     pointsize = 8)
pNhat_duiker / pNhat_impala +
  plot_annotation(tag_levels = c('a'), tag_prefix = '(',
                  tag_suffix = ')')
dev.off()



