# visualization of proein distributions 
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(here))


figpath2 = here("misc/figures/"); dir.create(figpath2)
datapath2 = here("misc/cd8plotdata/"); dir.create(datapath2)

s = readRDS("data/processed_data/ThymomaSeurat.20200422.clustered.mdprocessed.rds")
d = cbind(s@meta.data, s@reductions$umap@cell.embeddings, as.data.frame(t(s@assays$CITE@data)))
d = droplevels(d)
d$celltype = as.character(d$celltype)

rm(s);gc()
#######
#saveRDS(d, file = paste0(datapath2, "cd8data.rds"))
#df = readRDS(here("misc/cd8plotdata/cd8data.rds"))
#cd8 = unique(d$celltype)[grepl(pattern = "CD8", unique(d$celltype))]
#cd8 = unique(d$celltype)[grepl(pattern = "CD8", unique(d$celltype))] %>% as.character()
d1 = d %>% 
  filter(celltype %in% c("CD8Tcell_TEMRA", "CD8Tcell_TEM", "CD8Tcell_naive"))
  
p = ggplot(d1 %>% filter(CD45RA < 10 & CD8 > 5), aes(x = CD27, y = CD45RA)) +
  theme_bw() + 
  ggsci::scale_color_jama() + 
  geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype  = "dashed") + 
  geom_density_2d(aes(color = celltype), size = 0.8) + 
  guides(color = guide_legend(override.aes = list(size=2))) + 
  theme(legend.key.size = unit(0.3, units = "cm")) + 
  theme(axis.title = element_text(size = 20)) + 
  ggtitle('dsb normalized protein expression') + 
  theme(legend.position = c(0.6, 0.92)) 
ggsave(p, filename = paste0(figpath2, "1_cd8_protein_gating.pdf"), width = 3.5, height = 3.5)

# raw 
p = ggplot(d1 %>% filter(CD45RA < 10 & CD8 > 5), aes(x = CD27, y = CD45RA)) +
  theme_bw() + 
  geom_bin2d(bins = 200) + 
  scale_fill_viridis_c(option = 'B') +
  geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype  = "dashed") + 
  theme(axis.title = element_text(size = 20)) + 
  ggtitle('dsb normalized protein expression') +
  theme(legend.position = 'none') 
ggsave(p, filename = paste0(figpath2, "1.raw_cd8_protein_gating.pdf"), width = 3.5, height = 3.5)


# CD8 subset 2 
d2 = d %>% 
  filter(celltype %in% c("CD8Tcell_naive.CD27neg", "CD8Tcell_Activated", "CD8Tcell_CD161pos"))

p = ggplot(d2 %>% filter(CD8 > 5), aes(x = CD161, y = CD38)) +
  theme_bw() + 
  ggsci::scale_color_jama() + 
  geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype  = "dashed") + 
  geom_density_2d(aes(color = celltype), size = 0.8) + 
  guides(color = guide_legend(override.aes = list(size=2))) + 
  theme(legend.key.size = unit(0.3, units = "cm")) + 
  theme(axis.title = element_text(size = 20)) + 
  ggtitle('dsb normalized protein expression') + 
  theme(legend.position = c(0.6, 0.8)) 
ggsave(p, filename = paste0(figpath2, "2_cd8_protein_gating.pdf"), width = 3.5, height = 3.5)
## 
# raw 
p = ggplot(d2 %>% filter(CD8 > 5), aes(x = CD161, y = CD38)) +
  theme_bw() + 
  geom_bin2d(bins = 150) + 
  scale_fill_viridis_c(option = 'B') +
  geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype  = "dashed") + 
  theme(axis.title = element_text(size = 20)) + 
  ggtitle('dsb normalized protein expression') +
  theme(legend.position = 'none') 
ggsave(p, filename = paste0(figpath2, "2.raw_cd8_protein_gating.pdf"), width = 3.5, height = 3.5)


### CD4 subsets 
d1 = d %>% 
  filter(celltype %in% c("CD4Tcell_TM", "CD4Tcell_naive", "CD4Tcell_CM", "CD4Tcell_EM"))
p = ggplot(d1, aes(x = CD45RO, y = CD62L)) +
  theme_bw() + 
  ggsci::scale_color_jama() + 
  geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype  = "dashed") + 
  geom_density_2d(aes(color = celltype), size = 0.8) + 
  guides(color = guide_legend(override.aes = list(size=2))) + 
  theme(legend.key.size = unit(0.3, units = "cm")) + 
  theme(axis.title = element_text(size = 20)) + 
  ggtitle('dsb normalized protein expression') + 
  theme(legend.position = c(0.8, 0.8)) 
ggsave(p, filename = paste0(figpath2, "1_cd4_protein_gating.pdf"), width = 3.5, height = 3.5)

# raw 
p = ggplot(d1, aes(x = CD45RO, y = CD62L)) +
  theme_bw() + 
  geom_bin2d(bins = 200) + 
  scale_fill_viridis_c(option = 'B') +
  geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype  = "dashed") + 
  theme(axis.title = element_text(size = 20)) + 
  ggtitle('dsb normalized protein expression') +
  theme(legend.position = 'none') 
ggsave(p, filename = paste0(figpath2, "1.raw_cd4_protein_gating.pdf"), width = 3.5, height = 3.5)




### DC subsets 
df2 = d %>% filter(celltype %in% c( "DC_CD1c", "DC_mDC"))

p = ggplot(data = df2, aes(x = CD11b, y = CD1c)) +
  theme_bw() + 
  ggsci::scale_color_jama() + 
  geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype  = "dashed") + 
  geom_density_2d(aes(color = celltype), size = 0.8) + 
  guides(color = guide_legend(override.aes = list(size=2))) + 
  theme(axis.title = element_text(size = 20)) + 
  ggtitle('dsb normalized protein expression') + 
  theme(legend.position = c(0.8, 0.8)) 
ggsave(p, filename = paste0(figpath2, "DC_protein_gating.pdf"),width = 3.5, height = 3.5)

# raw 
p = ggplot(data = df2, aes(x = CD11b, y = CD1c)) +
  theme_bw() + 
  geom_bin2d(bins = 150) + 
  scale_fill_viridis_c(option = 'B') +
  geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype  = "dashed") + 
  theme(axis.title = element_text(size = 20)) + 
  ggtitle('dsb normalized protein expression') +
  theme(legend.position = 'none') 
ggsave(p, filename = paste0(figpath2, "raw.DC_protein_gating.pdf"), width = 3.5, height = 3.5)


df3 = df %>% filter(celltype %in% c("Mono_CD16pos", "Mono_classical"))
p = ggplot(data = df3, aes(x = CD16, y = CD14)) +
  theme_bw() + 
  ggsci::scale_color_d3() + 
  geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype  = "dashed") + 
  geom_density_2d(aes(color = celltype), size = 0.8) + 
  guides(color = guide_legend(override.aes = list(size=2))) + 
  theme(axis.title = element_text(size = 20)) + 
  ggtitle('dsb normalized protein expression') + 
  theme(legend.position = c(0.7, 0.8)) 
p
ggsave(p, filename = paste0(figpath2, "monocyte_protein_gating.pdf"), width = 3.5, height = 3.5)


# R version 3.6.1 (2019-07-05)
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] here_0.1        Seurat_3.1.5    forcats_0.5.0   stringr_1.4.0   dplyr_1.0.3     purrr_0.3.3    
# [7] readr_1.3.1     tidyr_1.0.2     tibble_3.1.0    ggplot2_3.3.0   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15          colorspace_1.4-1    ellipsis_0.3.1      ggridges_0.5.2      rprojroot_1.3-2    
# [6] fs_1.3.2            rstudioapi_0.11     farver_2.0.3        leiden_0.3.3        listenv_0.8.0      
# [11] npsurv_0.4-0        ggrepel_0.8.2       fansi_0.4.2         lubridate_1.7.4     xml2_1.3.2         
# [16] codetools_0.2-16    splines_3.6.1       lsei_1.2-0          jsonlite_1.6.1      broom_0.7.4        
# [21] ica_1.0-2           cluster_2.1.0       dbplyr_1.4.2        png_0.1-7           uwot_0.1.8         
# [26] sctransform_0.2.1   compiler_3.6.1      httr_1.4.1          backports_1.2.1     assertthat_0.2.1   
# [31] Matrix_1.2-18       lazyeval_0.2.2      cli_2.3.1           htmltools_0.4.0     tools_3.6.1        
# [36] rsvd_1.0.3          igraph_1.2.5        gtable_0.3.0        glue_1.4.2          RANN_2.6.1         
# [41] reshape2_1.4.3      Rcpp_1.0.4          cellranger_1.1.0    vctrs_0.3.6         gdata_2.18.0       
# [46] ape_5.3             nlme_3.1-145        lmtest_0.9-37       globals_0.12.5      rvest_0.3.5        
# [51] lifecycle_1.0.0     irlba_2.3.3         gtools_3.8.1        future_1.16.0       MASS_7.3-51.5      
# [56] zoo_1.8-7           scales_1.1.0        hms_0.5.3           parallel_3.6.1      RColorBrewer_1.1-2 
# [61] reticulate_1.14     pbapply_1.4-2       gridExtra_2.3       stringi_1.4.6       caTools_1.18.0     
# [66] rlang_0.4.10        pkgconfig_2.0.3     bitops_1.0-6        lattice_0.20-40     ROCR_1.0-7         
# [71] labeling_0.3        patchwork_1.0.0     htmlwidgets_1.5.1   cowplot_1.0.0       tidyselect_1.1.0   
# [76] ggsci_2.9           RcppAnnoy_0.0.16    plyr_1.8.6          magrittr_2.0.2      R6_2.4.1           
# [81] gplots_3.0.3        generics_0.0.2      DBI_1.1.0           pillar_1.5.0        haven_2.2.0        
# [86] withr_2.4.0         fitdistrplus_1.0-14 survival_3.1-11     future.apply_1.4.0  tsne_0.1-3         
# [91] modelr_0.1.6        crayon_1.4.1        KernSmooth_2.23-16  utf8_1.1.4          plotly_4.9.2       
# [96] grid_3.6.1          readxl_1.3.1        isoband_0.2.0       data.table_1.12.8   reprex_0.3.0       
# [101] digest_0.6.27       munsell_0.5.0       viridisLite_0.3.0  

