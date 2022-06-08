# mtor  analysis 
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(here))
suppressMessages(library(scglmmr))

datapath = here("cd8/mtor_cd8/generated_datav4/"); dir.create(datapath)
figpath = here("cd8/mtor_cd8/figuresV4//"); dir.create(figpath)


s = readRDS(here("data/processed_data/ThymomaSeurat.20200422.clustered.mdprocessed.rds"))
ct = unique(s@meta.data$celltype)
tc = ct[grep(x = ct, pattern = "DPT|DNT|CD8|CD4")] %>% as.vector()
Idents(s) = "celltype"
s = SubsetData(s, ident.use = tc)
s = NormalizeData(s,normalization.method = 'LogNormalize',assay = 'RNA')

# mtor pathways 
mtor.sigs = readRDS(here('cd8/mg_cd8/generated_dataV4/mtor.sigs2.rds'))


# format metadata as factors group_id is order leveled for: 
# contrast_fit = contrast(emm1, method = list( (c21 - c20) - (c11 - c10) ))
md = s@meta.data %>% 
  filter(timepoint %in% c(0,1)) %>% 
  mutate(group_id = paste(irae, timepoint, sep = '_')) %>% 
  mutate(group_id = factor(group_id,  levels = c('0_0', '0_1', '1_0', '1_1'))) %>%
  mutate(subjectid = factor(sampleid)) %>% 
  select(celltype, subjectid, sex = gender, age, group_id) %>% 
  mutate(age = as.numeric(age)) %>% 
  droplevels()

# qc data to remove celltypes with no cells for some subhects at both timepoints 
# keeps MLE more stable for the estimates of random intercept 
ct.si = apply(table(md$celltype, md$subjectid) , 1, min) 
c.keep = names(ct.si[ct.si > 7])
md = md[md$celltype %in% c.keep, ]


# add single cell weighted module scores
# split to standardize within cell type 
ct.md = split(md, f = md$celltype)
mod_scores = lapply(ct.md, function(x){ 
  scglmmr::WeightedCellModuleScore(gene_matrix = s@assays$RNA@data[ ,rownames(x)], 
                                   module_list = mtor.sigs, 
                                   threshold = 0,
                                   # standardize within protein celltype
                                   cellwise_scaling = TRUE, 
                                   return_weighted = FALSE )
  })
ms = bind_rows(mod_scores)

# correctly order rows after the split. 
ms = ms[match(x = rownames(md), table = rownames(ms)), ]
stopifnot(all.equal(rownames(ms), rownames(md)))

# rm object 
rm(s); gc()

# Fit mixed model 
# init save path for figures 
plot_savepath = paste0(figpath, "/marginalmeans.m1/"); dir.create(plot_savepath)

# specify model 
f1 = 'modulescore ~ group_id + (1|subjectid)'

# fit sc mod mixed model on ewighted module scores. 
mm_res.m1 = scglmmr::FitLmerContrast(module_data_frame = ms, 
                                     celltype_column = 'celltype', 
                                     metadata = md, 
                                     lmer_formula = f1, 
                                     plotdatqc = TRUE, 
                                     fixed_effects = NULL,
                                     figpath = plot_savepath)

saveRDS(mm_res.m1,file = paste0(datapath, "mm_res.m1.rds"))


sessionInfo()
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scglmmr_0.1.0   here_0.1        Seurat_3.1.5    forcats_0.5.0   stringr_1.4.0   dplyr_1.0.3    
# [7] purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_3.1.0    ggplot2_3.3.0   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.1.4               reticulate_1.14          tidyselect_1.1.0         lme4_1.1-21             
# [5] RSQLite_2.2.0            AnnotationDbi_1.48.0     htmlwidgets_1.5.1        grid_3.6.1              
# [9] BiocParallel_1.20.1      Rtsne_0.15               munsell_0.5.0            codetools_0.2-16        
# [13] ica_1.0-2                future_1.16.0            withr_2.4.0              colorspace_1.4-1        
# [17] GOSemSim_2.12.1          Biobase_2.46.0           rstudioapi_0.11          stats4_3.6.1            
# [21] ROCR_1.0-7               ggsignif_0.6.0           DOSE_3.12.0              listenv_0.8.0           
# [25] labeling_0.3             emmeans_1.4.5            urltools_1.7.3           polyclip_1.10-0         
# [29] pheatmap_1.0.12          bit64_0.9-7              farver_2.0.3             rprojroot_1.3-2         
# [33] TH.data_1.0-10           coda_0.19-4              vctrs_0.3.6              generics_0.0.2          
# [37] doParallel_1.0.15        R6_2.4.1                 graphlayouts_0.7.0       rsvd_1.0.3              
# [41] locfit_1.5-9.4           bitops_1.0-6             fgsea_1.12.0             gridGraphics_0.5-0      
# [45] assertthat_0.2.1         promises_1.1.0           scales_1.1.0             multcomp_1.4-12         
# [49] ggraph_2.0.2             enrichplot_1.6.1         gtable_0.3.0             npsurv_0.4-0            
# [53] egg_0.4.5                globals_0.12.5           sandwich_2.5-1           tidygraph_1.1.2         
# [57] rlang_0.4.10             slanter_0.2-0            splines_3.6.1            lazyeval_0.2.2          
# [61] broom_0.7.4              europepmc_0.3            BiocManager_1.30.10      reshape2_1.4.3          
# [65] modelr_0.1.6             backports_1.2.1          httpuv_1.5.2             qvalue_2.18.0           
# [69] clusterProfiler_3.14.3   tools_3.6.1              ggplotify_0.0.5          ellipsis_0.3.1          
# [73] gplots_3.0.3             RColorBrewer_1.1-2       BiocGenerics_0.32.0      ggridges_0.5.2          
# [77] Rcpp_1.0.4               plyr_1.8.6               progress_1.2.2           RCurl_1.98-1.1          
# [81] prettyunits_1.1.1        ggpubr_0.2.5             pbapply_1.4-2            viridis_0.5.1           
# [85] cowplot_1.0.0            S4Vectors_0.24.3         zoo_1.8-7                haven_2.2.0             
# [89] ggrepel_0.8.2            cluster_2.1.0            colorRamps_2.3           fs_1.3.2                
# [93] variancePartition_1.16.1 magrittr_2.0.2           data.table_1.12.8        DO.db_2.9               
# [97] lmtest_0.9-37            triebeard_0.3.0          reprex_0.3.0             RANN_2.6.1              
# [101] mvtnorm_1.1-0            fitdistrplus_1.0-14      hms_0.5.3                patchwork_1.0.0         
# [105] lsei_1.2-0               mime_0.9                 GSVA_1.34.0              xtable_1.8-4            
# [109] pbkrtest_0.4-8.6         XML_3.99-0.3             readxl_1.3.1             IRanges_2.20.2          
# [113] gridExtra_2.3            compiler_3.6.1           KernSmooth_2.23-16       crayon_1.4.1            
# [117] minqa_1.2.4              htmltools_0.4.0          later_1.0.0              geneplotter_1.64.0      
# [121] lubridate_1.7.4          DBI_1.1.0                tweenr_1.0.1             dbplyr_1.4.2            
# [125] MASS_7.3-51.5            boot_1.3-24              Matrix_1.2-18            cli_2.3.1               
# [129] gdata_2.18.0             parallel_3.6.1           igraph_1.2.5             pkgconfig_2.0.3         
# [133] rvcheck_0.1.8            plotly_4.9.2             foreach_1.4.8            xml2_1.3.2              
# [137] annotate_1.64.0          GeneOverlap_1.22.0       estimability_1.3         rvest_0.3.5             
# [141] digest_0.6.27            sctransform_0.2.1        RcppAnnoy_0.0.16         tsne_0.1-3              
# [145] graph_1.64.0             cellranger_1.1.0         leiden_0.3.3             fastmatch_1.1-0         
# [149] uwot_0.1.8               edgeR_3.28.1             GSEABase_1.48.0          shiny_1.4.0.2           
# [153] gtools_3.8.1             nloptr_1.2.2.1           lifecycle_1.0.0          nlme_3.1-145            
# [157] jsonlite_1.6.1           viridisLite_0.3.0        limma_3.42.2             fansi_0.4.2             
# [161] pillar_1.5.0             lattice_0.20-40          fastmap_1.0.1            httr_1.4.1              
# [165] survival_3.1-11          GO.db_3.10.0             glue_1.4.2               iterators_1.0.12        
# [169] png_0.1-7                shinythemes_1.1.2        bit_1.1-15.2             ggforce_0.3.1           
# [173] stringi_1.4.6            blob_1.2.1               org.Hs.eg.db_3.10.0      caTools_1.18.0          
# [177] memoise_1.1.0            irlba_2.3.3              future.apply_1.4.0       ape_5.3  