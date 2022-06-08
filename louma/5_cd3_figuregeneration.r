# luoma figure generation 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(here))
suppressMessages(library(magrittr))
figpath = here("louma/figures_v4/"); dir.create(figpath)

# read processed object
s = readRDS(file = here('louma/generated_data/cd3_louma_Seurat_processed.rds'))

# umap  
boxbox = list(
  theme_bw(), 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank())
  )

# plot with vector graphics
p = AugmentPlot(
  DimPlot(s, reduction = 'umap', group.by = 'seurat_clusters',
          label = TRUE, repel = TRUE, label.size = 18.5,
          cols = ggsci::pal_d3(palette = 'category20',alpha = 0.8)(20))
  ) + 
  boxbox + NoLegend() + xlab('UMAP 1') + ylab('UMAP 2') 
ggsave(p, filename = paste0(figpath,'umap_luoma.pdf'), width = 2.5, height = 2.5)


# m.sig 
md = s@meta.data 

# distribution of mtor signature expression
p = ggplot(md, aes(x = reorder(seurat_clusters, mtor_thymoma), 
                   y = mtor_thymoma, fill = seurat_clusters, color = seurat_clusters)) +
  theme_bw() + 
  coord_flip() + 
  geom_boxplot(width=0.7, show.legend = F, outlier.size = 0.5) +
  scale_color_manual(values = ggsci::pal_d3(palette = 'category20')(20)) + 
  scale_fill_manual(values = ggsci::pal_d3(palette = 'category20',alpha = 0.5)(20)) + 
  xlab("cluster") + ylab("Thymoma mTOR Sig.\n CD8 TEMRA ")
ggsave(p,filename = paste0(figpath, 'temrasig.mtor_thymoma.pdf'), width = 2, height = 3)

# calculate number of cells per donor per cluster 
test = reorder(md$seurat_clusters, md$mtor_thymoma)
c.order= levels(test)
cell.n = 
  md %>% 
  group_by(seurat_clusters, sample) %>% 
  tally(name = 'cell.n') 
cell.n$seurat_clusters = factor(cell.n$seurat_clusters, levels = c.order)

p = ggplot(cell.n, aes(x = seurat_clusters, cell.n, y = log10(cell.n),
                       fill = seurat_clusters, color = seurat_clusters)) + 
  theme_bw() + 
  coord_flip() + 
  geom_boxplot(width=0.7, show.legend = F, outlier.size = 0.5) + 
  scale_color_manual(values = ggsci::pal_d3(palette = 'category20')(20)) + 
  scale_fill_manual(values = ggsci::pal_d3(palette = 'category20',alpha = 0.5)(20)) + 
  ylab('log10 cell number\n per donor') + xlab('cluster')
ggsave(p,filename = paste0(figpath, 'celln.luoma.pdf'), width = 2, height = 3)
 


# heatmap of genes from luoma by cluster 
source(here("louma/luomagenes.r"))
genesub = genecombined[genecombined %in% rownames(s@assays$RNA@data)]
dat = s@assays$RNA@data[genesub, ] %>% 
  as.matrix() %>% 
  t() %>%
  as.data.frame() %>%
  add_column(s@meta.data)

# plot 
dd_av = dat %>% 
  group_by(seurat_clusters) %>% 
  summarize_at(.vars = genesub, .funs = mean) %>% 
  column_to_rownames("seurat_clusters") %>%
  t()

library(scglmmr)
HeatmapDiag(dd_av, scale = "row",
            color = rev(pals::brewer.puor(25)), 
            width = 3.1, height = 3.1, 
            border_color = 'NA',
            fontsize_row = 7, fontsize_col = 7,
            filename = paste0(figpath,"f_luomageneheatmap.pdf")
            )



#### Cell freq in clusters 

md = md %>% mutate(sample_irae = paste(IRAE, subjectid))
freq_plot =   
  md %>% 
  group_by(seurat_clusters, sample_irae) %>% 
  tally() %>% 
  ungroup() %>% 
  group_by(sample_irae) %>% 
  mutate(percent = n/sum(n)) %>% 
  ungroup() %>% 
  separate(col = sample_irae, into = c('irae', 'sample'), sep = " ") %>% 
  mutate(irae = factor(irae, levels = c('control', 'no', 'colitis'))) 

# frequency bar plot of clusters of interest from pseudobulk analysis 
my_comp = list( c1 = c('colitis', 'no'), c2 = c('colitis', 'control'))
freq_plot$irae = factor(freq_plot$irae, c('control', 'no', 'colitis') )

source('util/MattPMutils.r')
cc = c("grey", "#2a579b", "red")
cc2 = sapply(cc, col.alpha, 0.2) %>% unname()
freq_plot$cluster = paste('C' ,freq_plot$seurat_clusters,sep =  '')
p = ggplot(freq_plot %>% filter(cluster %in% c('C10',  'C2')),
           aes(x = irae, y = percent, fill = irae, color = irae)) + 
  theme_bw() + 
  boxbox + 
  geom_boxplot(outlier.shape = NA,show.legend = FALSE) +
  ylab("fraction of subject's total cells") + xlab("") + 
  geom_jitter(shape = 21, size = 1, width = 0.1,stroke = 0.5,  show.legend = F) + 
  ggpubr::stat_compare_means(comparisons = my_comp, method = 'wilcox', size = 2) + 
  scale_color_manual(values = cc) + 
  scale_fill_manual(values = cc2) + 
  theme(strip.background = element_blank()) +
  theme(axis.text.y = element_text(size=6 ))  +
  facet_wrap(~cluster, scales = 'free_y') + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, size = 10))
ggsave(p, filename = paste0(figpath, 'frequency.per.donor.c10.c2.c4.pdf'), width = 2.2, height = 2.7)

# sanity check to make sure fraction of donor total calculated corectly 
test1 =
  md %>%
  group_by(seurat_clusters, sample_irae) %>%
  tally() %>%
  ungroup() %>%
  group_by(sample_irae) %>%
  mutate(percent = n/sum(n)) %>%
  filter(sample_irae == 'colitis C1') %>% 
  summarize(pt = sum(percent))
stopifnot(isTRUE(all.equal(1, test1$pt)))


######################################################################################

cc = c("grey", "#2a579b", "red")
cc2 = sapply(cc, col.alpha, 0.2) %>% unname()
freq_plot$cluster = paste('C' ,freq_plot$seurat_clusters,sep =  '')
p = ggplot(freq_plot,aes(x = irae, y = percent, fill = irae, color = irae)) + 
  theme_bw() + 
  boxbox + 
  geom_boxplot(outlier.shape = NA,show.legend = FALSE) +
  ylab("fraction of subject's total cells") + xlab("") + 
  geom_jitter(shape = 21, size = 1, width = 0.1,stroke = 0.5,  show.legend = F) + 
  ggpubr::stat_compare_means(comparisons = my_comp, method = 'wilcox', size = 2) + 
  scale_color_manual(values = cc) + 
  scale_fill_manual(values = cc2) + 
  theme(strip.background = element_blank()) +
  theme(axis.text.y = element_text(size=6 ))  +
  facet_wrap(~cluster, scales = 'free_y', nrow = 2) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, size = 10))
ggsave(p, filename = paste0(figpath, 'FULL.frequency.pdf'), width = 8, height = 4.7)
######################################################################################




sessionInfo()
# R version 3.6.1 (2019-07-05)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] scglmmr_0.1.0   magrittr_2.0.2  here_0.1        Seurat_3.1.5    forcats_0.5.0  
# [6] stringr_1.4.0   dplyr_1.0.3     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2    
# [11] tibble_3.1.0    ggplot2_3.3.0   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.1.4               reticulate_1.14          tidyselect_1.1.0        
# [4] lme4_1.1-21              RSQLite_2.2.0            AnnotationDbi_1.48.0    
# [7] htmlwidgets_1.5.1        grid_3.6.1               BiocParallel_1.20.1     
# [10] Rtsne_0.15               munsell_0.5.0            codetools_0.2-16        
# [13] ica_1.0-2                future_1.16.0            withr_2.4.0             
# [16] colorspace_1.4-1         GOSemSim_2.12.1          Biobase_2.46.0          
# [19] rstudioapi_0.11          stats4_3.6.1             ROCR_1.0-7              
# [22] ggsignif_0.6.0           DOSE_3.12.0              listenv_0.8.0           
# [25] emmeans_1.4.5            labeling_0.3             urltools_1.7.3          
# [28] polyclip_1.10-0          pheatmap_1.0.12          bit64_0.9-7             
# [31] farver_2.0.3             rprojroot_1.3-2          TH.data_1.0-10          
# [34] coda_0.19-4              vctrs_0.3.6              generics_0.0.2          
# [37] doParallel_1.0.15        R6_2.4.1                 graphlayouts_0.7.0      
# [40] rsvd_1.0.3               locfit_1.5-9.4           pals_1.6                
# [43] bitops_1.0-6             fgsea_1.12.0             gridGraphics_0.5-0      
# [46] assertthat_0.2.1         promises_1.1.0           scales_1.1.0            
# [49] multcomp_1.4-12          ggraph_2.0.2             enrichplot_1.6.1        
# [52] gtable_0.3.0             npsurv_0.4-0             egg_0.4.5               
# [55] globals_0.12.5           sandwich_2.5-1           tidygraph_1.1.2         
# [58] rlang_0.4.10             slanter_0.2-0            splines_3.6.1           
# [61] lazyeval_0.2.2           dichromat_2.0-0          broom_0.7.4             
# [64] europepmc_0.3            BiocManager_1.30.10      reshape2_1.4.3          
# [67] modelr_0.1.6             backports_1.2.1          httpuv_1.5.2            
# [70] qvalue_2.18.0            clusterProfiler_3.14.3   tools_3.6.1             
# [73] ggplotify_0.0.5          ellipsis_0.3.1           gplots_3.0.3            
# [76] RColorBrewer_1.1-2       BiocGenerics_0.32.0      ggridges_0.5.2          
# [79] Rcpp_1.0.4               plyr_1.8.6               progress_1.2.2          
# [82] RCurl_1.98-1.1           prettyunits_1.1.1        ggpubr_0.2.5            
# [85] pbapply_1.4-2            viridis_0.5.1            cowplot_1.0.0           
# [88] S4Vectors_0.24.3         zoo_1.8-7                haven_2.2.0             
# [91] ggrepel_0.8.2            cluster_2.1.0            colorRamps_2.3          
# [94] fs_1.3.2                 variancePartition_1.16.1 data.table_1.12.8       
# [97] DO.db_2.9                lmtest_0.9-37            triebeard_0.3.0         
# [100] reprex_0.3.0             RANN_2.6.1               mvtnorm_1.1-0           
# [103] fitdistrplus_1.0-14      hms_0.5.3                patchwork_1.0.0         
# [106] lsei_1.2-0               mime_0.9                 GSVA_1.34.0             
# [109] xtable_1.8-4             pbkrtest_0.4-8.6         XML_3.99-0.3            
# [112] readxl_1.3.1             IRanges_2.20.2           gridExtra_2.3           
# [115] compiler_3.6.1           maps_3.3.0               KernSmooth_2.23-16      
# [118] crayon_1.4.1             minqa_1.2.4              htmltools_0.4.0         
# [121] later_1.0.0              geneplotter_1.64.0       lubridate_1.7.4         
# [124] DBI_1.1.0                tweenr_1.0.1             dbplyr_1.4.2            
# [127] MASS_7.3-51.5            boot_1.3-24              Matrix_1.2-18           
# [130] cli_2.3.1                gdata_2.18.0             parallel_3.6.1          
# [133] igraph_1.2.5             pkgconfig_2.0.3          rvcheck_0.1.8           
# [136] plotly_4.9.2             foreach_1.4.8            xml2_1.3.2              
# [139] annotate_1.64.0          GeneOverlap_1.22.0       estimability_1.3        
# [142] rvest_0.3.5              digest_0.6.27            sctransform_0.2.1       
# [145] RcppAnnoy_0.0.16         tsne_0.1-3               graph_1.64.0            
# [148] cellranger_1.1.0         leiden_0.3.3             fastmatch_1.1-0         
# [151] uwot_0.1.8               edgeR_3.28.1             GSEABase_1.48.0         
# [154] shiny_1.4.0.2            gtools_3.8.1             nloptr_1.2.2.1          
# [157] lifecycle_1.0.0          nlme_3.1-145             jsonlite_1.6.1          
# [160] mapproj_1.2.7            limma_3.42.2             viridisLite_0.3.0       
# [163] fansi_0.4.2              pillar_1.5.0             ggsci_2.9               
# [166] lattice_0.20-40          fastmap_1.0.1            httr_1.4.1              
# [169] survival_3.1-11          GO.db_3.10.0             glue_1.4.2              
# [172] iterators_1.0.12         png_0.1-7                shinythemes_1.1.2       
# [175] bit_1.1-15.2             ggforce_0.3.1            stringi_1.4.6           
# [178] blob_1.2.1               org.Hs.eg.db_3.10.0      caTools_1.18.0          
# [181] memoise_1.1.0            irlba_2.3.3              future.apply_1.4.0      
# [184] ape_5.3   
