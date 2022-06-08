suppressMessages(library(scglmmr))
suppressMessages(library(here))
suppressMessages(library(tidyverse))

# set paths 
datapath = here('louma/generated_data_v4/'); dir.create(datapath)
figpath = here('louma/figures_v4/'); dir.create(figpath)

# plot global
boxbox = list(
  theme_bw(), 
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), 
        axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank())
)


# read results 
fit = readRDS(here('louma/generated_data_scglmmr/fit_object.rds'))
av = readRDS(file = here('louma/generated_data_scglmmr/av_object.rds'))

# calc average signature expression for each subject 
mtor.sigs2 = readRDS('cd8/mg_cd8/generated_dataV4/mtor.sigs2.rds')
d_mtor = lapply(av, function(x) Matrix::colMeans(x[intersect(mtor.sigs2$`TEMRA MTOR`, rownames(x)),  ])) 


d = do.call(rbind, d_mtor ) %>% 
  as.data.frame() %>% 
  rownames_to_column('cluster') %>% 
  gather(sample, av_exprs, C1:NC6) %>% 
  mutate(irae = if_else(str_sub(sample, 1,1) == 'N', 'no',
                        if_else(str_sub(sample, 2,2) == 'T', 'control', 'colitis')))
my_comp = list( c1 = c('colitis', 'no'), c2 = c('colitis', 'control'))

d$irae = factor(d$irae, c('control', 'no', 'colitis') )
d$cluster = str_sub(d$cluster, 2,-1)
source('util/MattPMutils.r')
cc = c("grey", "#2a579b", "red")
cc2 = sapply(cc, col.alpha, 0.2) %>% unname()
p = 
  ggplot(d %>% filter(cluster %in% c('C10','C4')), aes(x = irae, y = av_exprs, fill = irae, color = irae)) + 
  theme_bw() + 
  boxbox + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + 
  ylab("Thymoma IRAE mTOR-LE") + xlab("") + 
  geom_jitter(shape = 21, size = 1, width = 0.1,stroke = 0.5,  show.legend = F) + 
  ggpubr::stat_compare_means(comparisons = my_comp, method = 'wilcox', size = 2) + 
  scale_color_manual(values = cc) + 
  scale_fill_manual(values = cc2) + 
  theme(strip.background = element_blank()) +
  theme(axis.text.y = element_text(size=6 )) + 
  facet_wrap(~cluster, scales = 'free_y') + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, size = 10))
ggsave(p, filename = paste0(figpath, 'temra.sig.average.pdf'), width = 2.2, height = 2.7)


######################################################################################## 
p = 
  ggplot(d, aes(x = irae, y = av_exprs, fill = irae, color = irae)) + 
  theme_bw() + 
  boxbox + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + 
  ylab("Thymoma IRAE mTOR-LE") + xlab("") + 
  geom_jitter(shape = 21, size = 1, width = 0.1,stroke = 0.5,  show.legend = F) + 
  ggpubr::stat_compare_means(comparisons = my_comp, method = 'wilcox', size = 2) + 
  scale_color_manual(values = cc) + 
  scale_fill_manual(values = cc2) + 
  theme(strip.background = element_blank()) +
  theme(axis.text.y = element_text(size=6 )) + 
  facet_wrap(~cluster, scales = 'free_y', nrow = 2) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, size = 10))
ggsave(p, filename = paste0(figpath, 'FULLtemra.sig.average.pdf'), width = 8, height = 4.7)
########################################################################################






###  repeat for other mtor signatures 
p.aes = list(
  geom_jitter(shape = 21, size = 1, width = 0.1,stroke = 0.5,  show.legend = F), 
  ggpubr::stat_compare_means(comparisons = my_comp, method = 'wilcox', size = 2), 
    scale_color_manual(values = cc) ,
    scale_fill_manual(values = cc2) ,
    theme(strip.background = element_blank()),
    theme(axis.text.y = element_text(size=6 )) ,
    facet_wrap(~cluster, scales = 'free_y') ,
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5, size = 12))
)

# calc av score for each mtor module 
ms = lapply(av, function(x){ 
  r = list()
  for (i in 1:length(mtor.sigs2)) {
    r[[i]] = Matrix::colMeans(x[intersect(mtor.sigs2[[i]], rownames(x)),  ])
  }
  d = bind_rows(r) %>% as.data.frame()
  rownames(d) = names(mtor.sigs2)
  d = d %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column('sample') %>% 
    #gather(sample, av_exprs, C1:NC6) %>% 
    mutate(irae = if_else(str_sub(sample, 1,1) == 'N', 'no',
                          if_else(str_sub(sample, 2,2) == 'T', 'control', 'colitis')))
  return(d)
})
d2 = bind_rows(ms,.id = '.ids')
d2$irae = factor(d2$irae, c('control', 'no', 'colitis') )
d2$cluster = str_sub(d2$.ids, 2, -1)

# biocarta 
p = 
  ggplot(d2 %>%
           select(cluster, `Biocarta MTOR`, sample, irae) %>% 
           filter(cluster %in% c('C10','C4','C2')), 
         aes(x = irae, y = `Biocarta MTOR`, fill = irae, color = irae)) + 
  theme_bw() + 
  boxbox + 
  geom_boxplot(outlier.shape = NA,  show.legend = FALSE) + 
  ylab("Biocarta mTOR sig.") + xlab("") + 
  p.aes
ggsave(p, filename = paste0(figpath, 'biocarta.sig.average.pdf'), width = 2.7, height = 3)

######################################################################################## 
p = 
  ggplot(d2 %>%
           select(cluster, `Biocarta MTOR`, sample, irae),
         aes(x = irae, y = `Biocarta MTOR`, fill = irae, color = irae)) + 
  theme_bw() + 
  boxbox + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + 
  ylab("Thymoma IRAE mTOR-LE") + xlab("") + 
  geom_jitter(shape = 21, size = 1, width = 0.1,stroke = 0.5,  show.legend = F) + 
  ggpubr::stat_compare_means(comparisons = my_comp, method = 'wilcox', size = 2) + 
  scale_color_manual(values = cc) + 
  scale_fill_manual(values = cc2) + 
  theme(strip.background = element_blank()) +
  theme(axis.text.y = element_text(size=6 )) + 
  facet_wrap(~cluster, scales = 'free_y', nrow = 2) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, size = 10))
ggsave(p, filename = paste0(figpath, 'FULLBIOCARTA.sig.average.pdf'), width = 8, height = 4.7)
########################################################################################


# kegg 
p = 
  ggplot(d2 %>%
           select(cluster, `KEGG MTOR`, sample, irae) %>% 
           filter(cluster %in% c('C10','C4','C2')), 
         aes(x = irae, y = `KEGG MTOR`, fill = irae, color = irae)) + 
  theme_bw() + 
  boxbox + 
  geom_boxplot(outlier.shape = NA,  show.legend = FALSE) + 
  ylab("KEGG mTOR sig.") + xlab("") + 
  p.aes
ggsave(p, filename = paste0(figpath, 'kegg.sig.average.pdf'), width = 2.7, height = 3)

######################################################################################## 
p = 
  ggplot(d2 %>%
           select(cluster, `KEGG MTOR`, sample, irae),
         aes(x = irae, y = `KEGG MTOR`, fill = irae, color = irae)) + 
  theme_bw() + 
  boxbox + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + 
  ylab("Thymoma IRAE mTOR-LE") + xlab("") + 
  geom_jitter(shape = 21, size = 1, width = 0.1,stroke = 0.5,  show.legend = F) + 
  ggpubr::stat_compare_means(comparisons = my_comp, method = 'wilcox', size = 2) + 
  scale_color_manual(values = cc) + 
  scale_fill_manual(values = cc2) + 
  theme(strip.background = element_blank()) +
  theme(axis.text.y = element_text(size=6 )) + 
  facet_wrap(~cluster, scales = 'free_y', nrow = 2) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, size = 10))
ggsave(p, filename = paste0(figpath, 'FULLKEGG.sig.average.pdf'), width = 8, height = 4.7)
########################################################################################



# reactome
p = 
  ggplot(d2 %>%
           select(cluster,`Reactome MTOR`, sample, irae) %>% 
           filter(cluster %in% c('C10','C4','C2')), 
         aes(x = irae, y = `Reactome MTOR`, fill = irae, color = irae)) + 
  theme_bw() + 
  boxbox + 
  geom_boxplot(outlier.shape = NA,  show.legend = FALSE) + 
  ylab("Reactome mTOR sig.") + xlab("") + 
  p.aes
ggsave(p, filename = paste0(figpath, 'reactome.sig.average.pdf'), width = 2.7, height = 3)


######################################################################################## 
p = 
  ggplot(d2 %>%
           select(cluster, `Reactome MTOR`, sample, irae),
         aes(x = irae, y = `Reactome MTOR`, fill = irae, color = irae)) + 
  theme_bw() + 
  boxbox + 
  geom_boxplot(outlier.shape = NA, show.legend = FALSE) + 
  ylab("Thymoma IRAE mTOR-LE") + xlab("") + 
  geom_jitter(shape = 21, size = 1, width = 0.1,stroke = 0.5,  show.legend = F) + 
  ggpubr::stat_compare_means(comparisons = my_comp, method = 'wilcox', size = 2) + 
  scale_color_manual(values = cc) + 
  scale_fill_manual(values = cc2) + 
  theme(strip.background = element_blank()) +
  theme(axis.text.y = element_text(size=6 )) + 
  facet_wrap(~cluster, scales = 'free_y', nrow = 2) + 
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, size = 10))
ggsave(p, filename = paste0(figpath, 'FULLREACTOME.sig.average.pdf'), width = 8, height = 4.7)
########################################################################################


sessionInfo()
# R version 3.6.1 (2019-07-05)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] forcats_0.5.0   stringr_1.4.0   dplyr_1.0.3     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2    
# [7] tibble_3.1.0    ggplot2_3.3.0   tidyverse_1.3.0 here_0.1        scglmmr_0.1.0  
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1             backports_1.2.1          fastmatch_1.1-0          plyr_1.8.6              
# [5] igraph_1.2.5             GSEABase_1.48.0          splines_3.6.1            BiocParallel_1.20.1     
# [9] TH.data_1.0-10           urltools_1.7.3           digest_0.6.27            foreach_1.4.8           
# [13] htmltools_0.4.0          GOSemSim_2.12.1          viridis_0.5.1            GO.db_3.10.0            
# [17] gdata_2.18.0             fansi_0.4.2              magrittr_2.0.2           memoise_1.1.0           
# [21] doParallel_1.0.15        limma_3.42.2             annotate_1.64.0          graphlayouts_0.7.0      
# [25] modelr_0.1.6             GeneOverlap_1.22.0       sandwich_2.5-1           enrichplot_1.6.1        
# [29] prettyunits_1.1.1        colorspace_1.4-1         rvest_0.3.5              blob_1.2.1              
# [33] ggrepel_0.8.2            haven_2.2.0              crayon_1.4.1             RCurl_1.98-1.1          
# [37] jsonlite_1.6.1           graph_1.64.0             lme4_1.1-21              iterators_1.0.12        
# [41] survival_3.1-11          zoo_1.8-7                glue_1.4.2               slanter_0.2-0           
# [45] polyclip_1.10-0          gtable_0.3.0             emmeans_1.4.5            BiocGenerics_0.32.0     
# [49] scales_1.1.0             DOSE_3.12.0              pheatmap_1.0.12          mvtnorm_1.1-0           
# [53] DBI_1.1.0                edgeR_3.28.1             Rcpp_1.0.4               viridisLite_0.3.0       
# [57] xtable_1.8-4             progress_1.2.2           gridGraphics_0.5-0       bit_1.1-15.2            
# [61] europepmc_0.3            stats4_3.6.1             GSVA_1.34.0              httr_1.4.1              
# [65] fgsea_1.12.0             gplots_3.0.3             RColorBrewer_1.1-2       ellipsis_0.3.1          
# [69] pkgconfig_2.0.3          XML_3.99-0.3             farver_2.0.3             dbplyr_1.4.2            
# [73] locfit_1.5-9.4           utf8_1.1.4               labeling_0.3             ggplotify_0.0.5         
# [77] tidyselect_1.1.0         rlang_0.4.10             reshape2_1.4.3           later_1.0.0             
# [81] AnnotationDbi_1.48.0     cellranger_1.1.0         munsell_0.5.0            tools_3.6.1             
# [85] cli_2.3.1                generics_0.0.2           RSQLite_2.2.0            broom_0.7.4             
# [89] ggridges_0.5.2           fastmap_1.0.1            fs_1.3.2                 org.Hs.eg.db_3.10.0     
# [93] bit64_0.9-7              tidygraph_1.1.2          caTools_1.18.0           ggraph_2.0.2            
# [97] egg_0.4.5                nlme_3.1-145             mime_0.9                 DO.db_2.9               
# [101] xml2_1.3.2               rstudioapi_0.11          pbkrtest_0.4-8.6         compiler_3.6.1          
# [105] shinythemes_1.1.2        variancePartition_1.16.1 ggsignif_0.6.0           reprex_0.3.0            
# [109] tweenr_1.0.1             geneplotter_1.64.0       stringi_1.4.6            lattice_0.20-40         
# [113] Matrix_1.2-18            nloptr_1.2.2.1           vctrs_0.3.6              pillar_1.5.0            
# [117] lifecycle_1.0.0          BiocManager_1.30.10      triebeard_0.3.0          estimability_1.3        
# [121] data.table_1.12.8        cowplot_1.0.0            bitops_1.0-6             colorRamps_2.3          
# [125] httpuv_1.5.2             qvalue_2.18.0            R6_2.4.1                 promises_1.1.0          
# [129] KernSmooth_2.23-16       gridExtra_2.3            IRanges_2.20.2           codetools_0.2-16        
# [133] assertthat_0.2.1         boot_1.3-24              MASS_7.3-51.5            gtools_3.8.1            
# [137] rprojroot_1.3-2          withr_2.4.0              multcomp_1.4-12          S4Vectors_0.24.3        
# [141] parallel_3.6.1           hms_0.5.3                clusterProfiler_3.14.3   grid_3.6.1              
# [145] coda_0.19-4              minqa_1.2.4              rvcheck_0.1.8            ggpubr_0.2.5            
# [149] ggforce_0.3.1            lubridate_1.7.4          Biobase_2.46.0           shiny_1.4.0.2   