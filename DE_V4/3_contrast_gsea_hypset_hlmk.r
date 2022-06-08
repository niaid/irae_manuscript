# gene set enrichment on contrast results 
# R 3.6 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(variancePartition))
suppressMessages(library(scglmmr))

# save path (dir created in script 1)
datapath = here("DE_V4/generated_data/")
figpath = here('DE_V4/figures/')

# parallel options 
library(BiocParallel)
BiocParallel::register(BiocParallel::SnowParam(4))
pparam = BiocParallel::SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# hallmark and hypothesis sets 
hlmk = readRDS(file = here('signature_curation/hallmark.rds'))
hypothesis_set = readRDS("signature_curation/hypothesis_set_1_V2.rds")

# read model results 
cont0 = readRDS(file = here('DE_V4/generated_data/cont0.rds'))
fit1 = readRDS(file = here('DE_V4/generated_data/fit1.rds'))

# gsea: baseline differences between groups - hypothesis set and hallmark pathway 
r0 = scglmmr::ExtractResult(model.fit.list = cont0, what = 'gene.t.ranks',
                            coefficient.number = 1, coef.name = 'baseline_irae')
hlmk.0 = RunFgseaOnRankList(rank.list.celltype = r0, pathways = hlmk, BPPARAM = pparam)
hyp.0 = RunFgseaOnRankList(rank.list.celltype = r0, pathways = hypothesis_set, BPPARAM = pparam)

# gsea: treatment effect (across all donors) - hypothesis set and hallmark pathway 
rtreat = scglmmr::ExtractResult(model.fit.list = fit1, what = 'gene.t.ranks',
                                coefficient.number = 3,coef.name = 'treatment')
hlmk.treat = RunFgseaOnRankList(rank.list.celltype = rtreat, pathways = hlmk, BPPARAM = pparam)
hyp.treat = RunFgseaOnRankList(rank.list.celltype = rtreat, pathways = hypothesis_set, BPPARAM = pparam)

# gsea: difference in treatment effect between groups (irae / no irae) hypothesis set and hallmark pathway 
rdelta = scglmmr::ExtractResult(model.fit.list = fit1,what = 'gene.t.ranks', 
                                coefficient.number = 2, coef.name = 'treatment_delta')
hlmk.delta = RunFgseaOnRankList(rank.list.celltype = rdelta, pathways = hlmk, BPPARAM = pparam)
hyp.delta = RunFgseaOnRankList(rank.list.celltype = rdelta, pathways = hypothesis_set, BPPARAM = pparam)

# save objects 
saveRDS(hlmk.0, file = paste0(datapath, 'hlmk.0.rds'))
saveRDS(hyp.0, file = paste0(datapath, 'hyp.0.rds')) 
saveRDS(hlmk.treat, file = paste0(datapath, 'hlmk.treat.rds')) 
saveRDS(hyp.treat, file = paste0(datapath, 'hyp.treat.rds')) 
saveRDS(hlmk.delta, file = paste0(datapath, 'hlmk.delta.rds')) 
saveRDS(hyp.delta, file = paste0(datapath, 'hyp.delta.rds')) 


sessionInfo()
# R version 3.6.1 (2019-07-05)
# 
# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] BiocParallel_1.20.1      scglmmr_0.1.0            variancePartition_1.16.1 Biobase_2.46.0           BiocGenerics_0.32.0     
# [6] scales_1.1.0             foreach_1.4.8            limma_3.42.2             magrittr_2.0.2           here_0.1                
# [11] forcats_0.5.0            stringr_1.4.0            dplyr_1.0.3              purrr_0.3.3              readr_1.3.1             
# [16] tidyr_1.0.2              tibble_3.1.0             ggplot2_3.3.0            tidyverse_1.3.0         
# 
# loaded via a namespace (and not attached):
# [1] utf8_1.1.4             tidyselect_1.1.0       lme4_1.1-21            RSQLite_2.2.0          AnnotationDbi_1.48.0   grid_3.6.1            
# [7] devtools_2.2.2         munsell_0.5.0          codetools_0.2-16       withr_2.4.0            colorspace_1.4-1       GOSemSim_2.12.1       
# [13] rstudioapi_0.11        stats4_3.6.1           ggsignif_0.6.0         DOSE_3.12.0            emmeans_1.4.5          urltools_1.7.3        
# [19] polyclip_1.10-0        bit64_0.9-7            farver_2.0.3           pheatmap_1.0.12        rprojroot_1.3-2        coda_0.19-4           
# [25] vctrs_0.3.6            generics_0.0.2         TH.data_1.0-10         R6_2.4.1               doParallel_1.0.15      graphlayouts_0.7.0    
# [31] locfit_1.5-9.4         bitops_1.0-6           fgsea_1.12.0           gridGraphics_0.5-0     assertthat_0.2.1       promises_1.1.0        
# [37] multcomp_1.4-12        ggraph_2.0.2           enrichplot_1.6.1       gtable_0.3.0           egg_0.4.5              processx_3.4.2        
# [43] tidygraph_1.1.2        sandwich_2.5-1         rlang_0.4.10           slanter_0.2-0          splines_3.6.1          broom_0.7.4           
# [49] europepmc_0.3          modelr_0.1.6           BiocManager_1.30.10    reshape2_1.4.3         backports_1.2.1        httpuv_1.5.2          
# [55] qvalue_2.18.0          clusterProfiler_3.14.3 tools_3.6.1            usethis_1.5.1          ggplotify_0.0.5        ellipsis_0.3.1        
# [61] gplots_3.0.3           RColorBrewer_1.1-2     sessioninfo_1.1.1      ggridges_0.5.2         Rcpp_1.0.4             plyr_1.8.6            
# [67] progress_1.2.2         RCurl_1.98-1.1         ps_1.3.2               prettyunits_1.1.1      ggpubr_0.2.5           viridis_0.5.1         
# [73] cowplot_1.0.0          S4Vectors_0.24.3       zoo_1.8-7              haven_2.2.0            ggrepel_0.8.2          colorRamps_2.3        
# [79] fs_1.3.2               data.table_1.12.8      DO.db_2.9              reprex_0.3.0           triebeard_0.3.0        mvtnorm_1.1-0         
# [85] pkgload_1.0.2          hms_0.5.3              mime_0.9               GSVA_1.34.0            xtable_1.8-4           pbkrtest_0.4-8.6      
# [91] XML_3.99-0.3           readxl_1.3.1           IRanges_2.20.2         gridExtra_2.3          testthat_2.3.2         compiler_3.6.1        
# [97] KernSmooth_2.23-16     crayon_1.4.1           minqa_1.2.4            htmltools_0.4.0        later_1.0.0            snow_0.4-3            
# [103] geneplotter_1.64.0     lubridate_1.7.4        DBI_1.1.0              tweenr_1.0.1           dbplyr_1.4.2           MASS_7.3-51.5         
# [109] boot_1.3-24            Matrix_1.2-18          cli_2.3.1              gdata_2.18.0           igraph_1.2.5           pkgconfig_2.0.3       
# [115] rvcheck_0.1.8          xml2_1.3.2             annotate_1.64.0        GeneOverlap_1.22.0     estimability_1.3       rvest_0.3.5           
# [121] callr_3.4.2            digest_0.6.27          graph_1.64.0           cellranger_1.1.0       fastmatch_1.1-0        edgeR_3.28.1          
# [127] GSEABase_1.48.0        curl_4.3               shiny_1.4.0.2          gtools_3.8.1           nloptr_1.2.2.1         lifecycle_1.0.0       
# [133] nlme_3.1-145           jsonlite_1.6.1         desc_1.2.0             viridisLite_0.3.0      fansi_0.4.2            pillar_1.5.0          
# [139] lattice_0.20-40        fastmap_1.0.1          httr_1.4.1             pkgbuild_1.0.6         survival_3.1-11        GO.db_3.10.0          
# [145] glue_1.4.2             remotes_2.1.1          shinythemes_1.1.2      iterators_1.0.12       bit_1.1-15.2           ggforce_0.3.1         
# [151] stringi_1.4.6          blob_1.2.1             org.Hs.eg.db_3.10.0    caTools_1.18.0         memoise_1.1.0     
