suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(Seurat))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
figpath = here("cd8/mg_cd8/figuresV4/"); dir.create(figpath)
datapath = here("cd8/mg_cd8/generated_datav4/"); dir.create(datapath)

# h1 hds 
h1 = readRDS(file = here("data/starting_data/kotliarov_etal_CITE-seq_data/h1_d0_singlets_ann_Seurat3.1_dsbnorm.rds"))
d = cbind(h1@meta.data, as.data.frame(t(as.matrix(h1@assays$CITE@data))))

# set goal mg plot aes  
plot_aes = list(theme_bw(), scale_fill_viridis_c(option = "B"))

# healthy flu gates 
p1 = ggplot(d, aes(x = `CD8-PROT`, y = `CD3-PROT`)) + 
  geom_bin2d(bins = 200, show.legend = F) + plot_aes
p2 = ggplot(d %>% filter(`CD8-PROT` > 7 & `CD3-PROT` > 4 & `CD62L-PROT` < 3.5),
       aes(x = `CD27-PROT`, y = `CD45RA-PROT`)) + 
  geom_bin2d(bins = 200, show.legend = F) + plot_aes

# thymoma data 
th = readRDS(file = here("data/processed_data/ThymomaSeurat.20200422.clustered.mdprocessed.rds"))
d2 = cbind(th@meta.data, as.data.frame(t(as.matrix(th@assays$CITE@data))))

# thymoma gates (use same thresholds) 
p3 = ggplot(d2 %>% filter(timepoint == '0') %>% 
              filter(CD3 < 20), aes(x = CD8, y = CD3)) + 
  geom_bin2d(bins = 200, show.legend = F) + 
  plot_aes
p4 = ggplot(d2 %>% 
              filter(timepoint == '0') %>%
              filter(CD8 > 7 & CD3 > 4 & CD62L < 3.5),
            aes(x = CD27, y = CD45RA)) + 
  geom_bin2d(bins = 200, show.legend = F) + 
  plot_aes

# plot of temra population
p6 = p2 + ggtitle("Healthy Donor \nCD62L-CD3+CD8+") + p4 + ggtitle('Thymoma \nCD62L-CD3+CD8+')
ggsave(p6, filename = paste0(figpath, 'TEMRA_gated_cd8s_flu_thymoma.pdf'), width = 5.5, height = 3.5)

# gate the temra-like cells from the H1 and thymoma data 
temra_h1 = d %>% 
  filter(`CD8-PROT` > 7 & `CD3-PROT` > 4 & `CD45RA-PROT` > 5 & `CD62L-PROT` < 3.5 & `CD27-PROT` < 5) %>% 
  rownames()
temra_th = d2 %>% 
  filter(timepoint == '0') %>% 
  filter(CD8 > 7 & CD3 > 4 & CD45RA > 2.5 & CD62L < 3.5 & CD27 < 5) %>%
  rownames()

# quick check of overlap with gated and unsupervised classificaiton 
# based on this added an additional CD62L filter. 
d2_sub = d2[temra_th, ]
table(d2_sub$celltype)

# unify metadata columns
th_meta = th@meta.data %>% select(sampleid, batch = Batch, classificaiton = irae) 
th_meta$sampleid = paste('thymoma', th_meta$sampleid, sep = '_')
th_meta$classificaiton = ifelse(th_meta$classificaiton == 1, 'IRAE', 'NO_IRAE')
th_meta$sampleid = paste(th_meta$sampleid,th_meta$classificaiton, sep = '_')
h1_meta = h1@meta.data %>% select(sampleid, batch)
h1_meta$classificaiton = 'healthy'
h1_meta$sampleid = paste('healthy',h1_meta$sampleid,sep = '_')

# merge metadata  
meta = rbind(h1_meta, th_meta)

# combine counts 
th_rna = th@assays$RNA@counts[ ,temra_th]
h1_rna = h1@assays$RNA@counts[ ,temra_h1]
umi = cbind(h1_rna, th_rna)
meta = meta[rownames(meta) %in% c(temra_h1, temra_th), ]
umi = umi[ ,match(colnames(umi), rownames(meta))]
stopifnot(all.equal(colnames(umi), rownames(meta)))

# QC contingency of cells by subject for each celltype 
meta$celltype = 'gated_TEMRA'

# remove genes with no counts in all sx
umi = umi[Matrix::rowSums(umi) > 0, ]

# make pseudobulk 
scell = lapply(X = split(meta, f = meta$sampleid), FUN = rownames)
csample = lapply(scell, function(x) Matrix::rowSums(umi[ ,x]))
pb = as.matrix(do.call(cbind, csample))

# make desin mat 
design = as.data.frame.matrix(table(meta$sampleid, meta$classificaiton))
design[design>0] = 1

# QC design matrix rows match the pseudobulk data columns
stopifnot(isTRUE(all.equal(target = colnames(pb), current = rownames(design))))
stopifnot(Matrix::rankMatrix(design) == ncol(design))
stopifnot(any(colSums(design) == 0) == FALSE)

# make contrast matrix 
c_mat =  limma::makeContrasts(
  thymoma_v_healthy = ((IRAE + NO_IRAE) / 2) - healthy, 
  irae_vs_no_irae = (IRAE - NO_IRAE), 
  irae_vs_healthy = (IRAE  - healthy), 
  levels = colnames(design)
  )

# specify genes to test 
g.test = readRDS(file = here('DE_V4/generated_data/cont0.rds'))
temra.g.test = rownames(g.test$CD8Tcell_TEMRA$coefficients)
saveRDS(temra.g.test, file = paste0(datapath,'temra.g.test'))

# normalize and filter features 
dge = edgeR::DGEList(pb)
dge = edgeR::calcNormFactors(dge, method = "RLE")
dge = dge[temra.g.test, keep.lib.sizes=FALSE]
v = limma::voom(counts = dge, design = design, normalize.method = "none", save.plot = T, plot = T)
fit = limma::lmFit(v, design = design)
c_fit = limma::contrasts.fit(fit, contrasts = c_mat)
eb = limma::eBayes(c_fit)

# save model and data 
saveRDS(design, file = paste0(datapath, 'design.rds'))
saveRDS(eb,file = paste0(datapath, 'eb.rds'))
saveRDS(pb, file = paste0(datapath,'pb.rds'))
saveRDS(scell, file = paste0(datapath,'scell.rds'))
saveRDS(umi, file = paste0(datapath,'umi.rds'))


sessionInfo()
# R version 3.6.1 (2019-07-05)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] Seurat_3.1.5    magrittr_2.0.2  scglmmr_0.1.0   here_0.1        forcats_0.5.0   stringr_1.4.0  
# [7] dplyr_1.0.3     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_3.1.0    ggplot2_3.3.0  
# [13] tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.1.4               reticulate_1.14          tidyselect_1.1.0         lme4_1.1-21             
# [5] htmlwidgets_1.5.1        RSQLite_2.2.0            AnnotationDbi_1.48.0     grid_3.6.1              
# [9] BiocParallel_1.20.1      Rtsne_0.15               munsell_0.5.0            ica_1.0-2               
# [13] codetools_0.2-16         future_1.16.0            withr_2.4.0              colorspace_1.4-1        
# [17] GOSemSim_2.12.1          Biobase_2.46.0           rstudioapi_0.11          stats4_3.6.1            
# [21] ROCR_1.0-7               ggsignif_0.6.0           DOSE_3.12.0              listenv_0.8.0           
# [25] labeling_0.3             emmeans_1.4.5            urltools_1.7.3           polyclip_1.10-0         
# [29] bit64_0.9-7              farver_2.0.3             pheatmap_1.0.12          rprojroot_1.3-2         
# [33] coda_0.19-4              vctrs_0.3.6              generics_0.0.2           TH.data_1.0-10          
# [37] R6_2.4.1                 doParallel_1.0.15        graphlayouts_0.7.0       rsvd_1.0.3              
# [41] locfit_1.5-9.4           bitops_1.0-6             fgsea_1.12.0             gridGraphics_0.5-0      
# [45] assertthat_0.2.1         promises_1.1.0           scales_1.1.0             multcomp_1.4-12         
# [49] ggraph_2.0.2             enrichplot_1.6.1         gtable_0.3.0             npsurv_0.4-0            
# [53] globals_0.12.5           egg_0.4.5                tidygraph_1.1.2          sandwich_2.5-1          
# [57] rlang_0.4.10             slanter_0.2-0            splines_3.6.1            lazyeval_0.2.2          
# [61] broom_0.7.4              europepmc_0.3            BiocManager_1.30.10      reshape2_1.4.3          
# [65] modelr_0.1.6             backports_1.2.1          httpuv_1.5.2             qvalue_2.18.0           
# [69] clusterProfiler_3.14.3   tools_3.6.1              ggplotify_0.0.5          ellipsis_0.3.1          
# [73] gplots_3.0.3             RColorBrewer_1.1-2       BiocGenerics_0.32.0      ggridges_0.5.2          
# [77] Rcpp_1.0.4               plyr_1.8.6               progress_1.2.2           RCurl_1.98-1.1          
# [81] prettyunits_1.1.1        ggpubr_0.2.5             pbapply_1.4-2            viridis_0.5.1           
# [85] cowplot_1.0.0            S4Vectors_0.24.3         zoo_1.8-7                cluster_2.1.0           
# [89] haven_2.2.0              ggrepel_0.8.2            colorRamps_2.3           fs_1.3.2                
# [93] variancePartition_1.16.1 data.table_1.12.8        DO.db_2.9                lmtest_0.9-37           
# [97] triebeard_0.3.0          reprex_0.3.0             RANN_2.6.1               mvtnorm_1.1-0           
# [101] fitdistrplus_1.0-14      patchwork_1.0.0          lsei_1.2-0               hms_0.5.3               
# [105] mime_0.9                 GSVA_1.34.0              xtable_1.8-4             pbkrtest_0.4-8.6        
# [109] XML_3.99-0.3             readxl_1.3.1             IRanges_2.20.2           gridExtra_2.3           
# [113] compiler_3.6.1           KernSmooth_2.23-16       crayon_1.4.1             minqa_1.2.4             
# [117] htmltools_0.4.0          later_1.0.0              geneplotter_1.64.0       lubridate_1.7.4         
# [121] DBI_1.1.0                tweenr_1.0.1             dbplyr_1.4.2             MASS_7.3-51.5           
# [125] boot_1.3-24              Matrix_1.2-18            cli_2.3.1                gdata_2.18.0            
# [129] parallel_3.6.1           igraph_1.2.5             pkgconfig_2.0.3          rvcheck_0.1.8           
# [133] plotly_4.9.2             xml2_1.3.2               foreach_1.4.8            annotate_1.64.0         
# [137] GeneOverlap_1.22.0       estimability_1.3         rvest_0.3.5              digest_0.6.27           
# [141] tsne_0.1-3               sctransform_0.2.1        RcppAnnoy_0.0.16         graph_1.64.0            
# [145] leiden_0.3.3             cellranger_1.1.0         fastmatch_1.1-0          uwot_0.1.8              
# [149] edgeR_3.28.1             GSEABase_1.48.0          shiny_1.4.0.2            gtools_3.8.1            
# [153] nloptr_1.2.2.1           lifecycle_1.0.0          nlme_3.1-145             jsonlite_1.6.1          
# [157] viridisLite_0.3.0        limma_3.42.2             fansi_0.4.2              pillar_1.5.0            
# [161] lattice_0.20-40          fastmap_1.0.1            httr_1.4.1               survival_3.1-11         
# [165] GO.db_3.10.0             glue_1.4.2               png_0.1-7                shinythemes_1.1.2       
# [169] iterators_1.0.12         bit_1.1-15.2             ggforce_0.3.1            stringi_1.4.6           
# [173] blob_1.2.1               org.Hs.eg.db_3.10.0      caTools_1.18.0           memoise_1.1.0           
# [177] irlba_2.3.3              future.apply_1.4.0       ape_5.3    





