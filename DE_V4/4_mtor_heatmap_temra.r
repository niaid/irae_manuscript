# visulaziaiton of gsea contrast results 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))


figpath = here('DE_V4/figures/msig/'); dir.create(figpath, recursive = TRUE)

# leading edge visualization of tema mTOR leading edge signature genes. 
meta = readRDS(file = here('data/processed_data/full_remapped_metadata_thymoma.rds'))
hlmk.0 = readRDS(file = here('DE_V4/generated_data/hlmk.0.rds'))
av = readRDS(file = here('DE_V4/generated_data/av.rds'))
le_expr = LeadEdgeTidySampleExprs(av.exprs.list = av, 
                                  gsea.list = hlmk.0, 
                                  padj.filter = 0.2, 
                                  NES.filter = -Inf)

####################
# customize the base scglmmr plot above with addtiional variables 
mat = 
  LeadEdgeSampleHeatmap(tidy.exprs.list = le_expr,
                        modulename = "HALLMARK_MTORC1_SIGNALING",
                        celltype_plot = "CD8Tcell_TEMRA",
                        metadata = meta, 
                        metadata_annotate = c('irae', 'timepoint'),
                        sample_column = 'sample',
                        returnmat = TRUE)
                        
# make custom annotation 
heatmap_anno = meta[meta$celltype == "CD8Tcell_TEMRA", c('sample', 'irae', 'timepoint')] %>% 
  dplyr::group_by(sample) %>%
  dplyr::summarise_each(list(~unique(.))) %>%
  tibble::column_to_rownames('sample')
heatmap_anno$irae = ifelse(heatmap_anno$irae == 1, "Yes", "No")

anno_color = list(
  irae = c('No' = 'dodgerblue', 'Yes' = 'red'), 
  timepoint = c('0' = "orange", "1" = "white")
)

cu = c("#053061", "#1E61A5", "#3C8ABE", "#7CB7D6", "#BAD9E9", "#E5EEF3", 
       "#F9EAE1", "#F9C7AD", "#EB9273", "#CF5246", "#AB1529", "#67001F")

library(dendsort)
sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
mat_cluster_cols <- hclust(dist(t(mat)))
mat_cluster_cols <- sort_hclust(mat_cluster_cols)

pheatmap::pheatmap(mat, 
                   scale = "row",
                   border_color = NA,
                   treeheight_row = 0,
                   treeheight_col = 18,
                   #cutree_cols = 3,
                   annotation = heatmap_anno,
                   annotation_colors = anno_color,
                   color = cu,
                   width = 5, 
                   height = 6.45, 
                   filename = paste0(figpath, "CLUST_Ledge_average_mtor_temra_multitime_1.pdf"))



sessionInfo()
# R version 3.6.1 (2019-07-05)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] dendsort_0.3.4  scglmmr_0.1.0   magrittr_2.0.2  here_0.1        forcats_0.5.0   stringr_1.4.0  
# [7] dplyr_1.0.3     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_3.1.0    ggplot2_3.3.0  
# [13] tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
# [1] readxl_1.3.1             backports_1.2.1          fastmatch_1.1-0          plyr_1.8.6              
# [5] igraph_1.2.5             GSEABase_1.48.0          splines_3.6.1            BiocParallel_1.20.1     
# [9] TH.data_1.0-10           urltools_1.7.3           digest_0.6.27            foreach_1.4.8           
# [13] htmltools_0.4.0          GOSemSim_2.12.1          viridis_0.5.1            GO.db_3.10.0            
# [17] gdata_2.18.0             fansi_0.4.2              memoise_1.1.0            doParallel_1.0.15       
# [21] limma_3.42.2             annotate_1.64.0          graphlayouts_0.7.0       modelr_0.1.6            
# [25] GeneOverlap_1.22.0       sandwich_2.5-1           enrichplot_1.6.1         prettyunits_1.1.1       
# [29] colorspace_1.4-1         blob_1.2.1               rvest_0.3.5              ggrepel_0.8.2           
# [33] haven_2.2.0              crayon_1.4.1             RCurl_1.98-1.1           jsonlite_1.6.1          
# [37] graph_1.64.0             lme4_1.1-21              iterators_1.0.12         survival_3.1-11         
# [41] zoo_1.8-7                glue_1.4.2               slanter_0.2-0            polyclip_1.10-0         
# [45] gtable_0.3.0             emmeans_1.4.5            BiocGenerics_0.32.0      scales_1.1.0            
# [49] DOSE_3.12.0              pheatmap_1.0.12          mvtnorm_1.1-0            edgeR_3.28.1            
# [53] DBI_1.1.0                Rcpp_1.0.4               viridisLite_0.3.0        xtable_1.8-4            
# [57] progress_1.2.2           gridGraphics_0.5-0       bit_1.1-15.2             europepmc_0.3           
# [61] stats4_3.6.1             GSVA_1.34.0              httr_1.4.1               fgsea_1.12.0            
# [65] gplots_3.0.3             RColorBrewer_1.1-2       ellipsis_0.3.1           pkgconfig_2.0.3         
# [69] XML_3.99-0.3             farver_2.0.3             dbplyr_1.4.2             locfit_1.5-9.4          
# [73] utf8_1.1.4               ggplotify_0.0.5          tidyselect_1.1.0         rlang_0.4.10            
# [77] reshape2_1.4.3           later_1.0.0              AnnotationDbi_1.48.0     munsell_0.5.0           
# [81] cellranger_1.1.0         tools_3.6.1              cli_2.3.1                generics_0.0.2          
# [85] RSQLite_2.2.0            ggridges_0.5.2           broom_0.7.4              fastmap_1.0.1           
# [89] org.Hs.eg.db_3.10.0      bit64_0.9-7              fs_1.3.2                 tidygraph_1.1.2         
# [93] caTools_1.18.0           ggraph_2.0.2             egg_0.4.5                nlme_3.1-145            
# [97] mime_0.9                 DO.db_2.9                xml2_1.3.2               pbkrtest_0.4-8.6        
# [101] compiler_3.6.1           shinythemes_1.1.2        rstudioapi_0.11          variancePartition_1.16.1
# [105] ggsignif_0.6.0           reprex_0.3.0             tweenr_1.0.1             geneplotter_1.64.0      
# [109] stringi_1.4.6            lattice_0.20-40          Matrix_1.2-18            nloptr_1.2.2.1          
# [113] vctrs_0.3.6              pillar_1.5.0             lifecycle_1.0.0          BiocManager_1.30.10     
# [117] triebeard_0.3.0          estimability_1.3         data.table_1.12.8        cowplot_1.0.0           
# [121] bitops_1.0-6             colorRamps_2.3           httpuv_1.5.2             qvalue_2.18.0           
# [125] R6_2.4.1                 promises_1.1.0           KernSmooth_2.23-16       gridExtra_2.3           
# [129] IRanges_2.20.2           codetools_0.2-16         boot_1.3-24              MASS_7.3-51.5           
# [133] gtools_3.8.1             assertthat_0.2.1         rprojroot_1.3-2          withr_2.4.0             
# [137] multcomp_1.4-12          S4Vectors_0.24.3         parallel_3.6.1           hms_0.5.3               
# [141] clusterProfiler_3.14.3   grid_3.6.1               coda_0.19-4              minqa_1.2.4             
# [145] rvcheck_0.1.8            ggpubr_0.2.5             ggforce_0.3.1            Biobase_2.46.0          
# [149] shiny_1.4.0.2            lubridate_1.7.4   
# 
# 
