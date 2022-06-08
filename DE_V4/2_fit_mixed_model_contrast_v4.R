# DE analysis with contrasts 
# This script uses R version 4.0.5
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(variancePartition))
suppressMessages(library(scglmmr))

# save path (dir created in script 1)
datapath = here("DE_V4/generated_data/")
figpath = here('DE_V4/figures/')

# parallel options 
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# load data 
pb = readRDS(file = here('DE_V4/generated_data/pb.rds'))
sampleids = colnames(pb[[1]])

samplemd = 
  data.table::fread(file = here('data/starting_data/subject_data/meta_table.txt')) %>% 
  as.data.frame() %>% 
  mutate(sx = paste(sampleid, timepoint, sep = '_')) %>% 
  rename(subjectid = sampleid) %>% 
  mutate(irae_timepoint = paste(irae, timepoint, sep = '_')) %>% 
  filter(sx %in% sampleids) %>% 
  column_to_rownames('sx') 

# Set levels of cimbined factor 
# number_number is group_time; 1_1 = IRAE_POSTTREATMENT 0_0 is NO IRAE _ BASELINE 
samplemd$irae.time = factor(samplemd$irae_timepoint, levels = c("1_0", "1_1", "0_0", "0_1"))


# specify model 
f1 <- ~ 0 + irae.time + (1|subjectid) 

# make contrast matrix 
L2 = makeContrastsDream(
  formula = f1, 
  data = samplemd,
  contrasts = c(
    baseline = "irae.time1_0 - irae.time0_0", # note fixef model also fit for this single time point. 
    treatment_delta = "( irae.time1_1 - irae.time1_0 ) - ( irae.time0_1 - irae.time0_0 )",
    treatment = "( irae.time1_1 + irae.time0_1 ) / 2 - ( irae.time1_0 + irae.time0_0 ) / 2 "
    )
  )
plotContrasts(L2) + ggsave(filename = paste0(figpath,'contrastmodel.pdf'), width = 7, height = 4)


# fit model on each subset 
# init store 
fit1 = v1 = list()
for (i in 1:length(pb)) {

  # init data 
  meta = samplemd
  form = f1 
  contrast_matrix = L2
  counts = pb[[i]]
  
  # dge list 
  d = edgeR::DGEList(counts = counts, samples = meta)
  
  # filter cell type specific lowly expressed genes and calc norm factors 
  gtable = edgeR::filterByExpr(y = d$counts, 
                               min.count = 3,
                               design = as.factor(d$samples$irae.time)
                               )
  print(names(pb)[i]);print(table(gtable))
  d = d[gtable, keep.lib.sizes=FALSE]
  d = edgeR::calcNormFactors(object = d)
  
  # get voom observation level weights 
  v = voomWithDreamWeights(counts = d, 
                           formula = form,
                           data = meta, 
                           BPPARAM = pparam, 
                           plot = TRUE, save.plot = TRUE)
  # fit contrast mixed model 
  fitmm = dream(exprObj = v, 
                formula = form,
                data = meta,
                L = contrast_matrix,
                BPPARAM = pparam, 
                useWeights = TRUE, REML = TRUE)
  # save results 
  v1[[i]] = v
  fit1[[i]] = fitmm
}
names(v1) = names(fit1) = names(pb)



# baseline fixed effects model 
# contrast t0
irae.time = samplemd$irae.time
design.matrix = model.matrix( ~ 0 +  irae.time)
c0 = makeContrasts(
  baseline_irae = (irae.time1_0 - irae.time0_0),
  levels = colnames(design.matrix)
  )
# check contrast 
stopifnot(isTRUE(
  all.equal(
    as.numeric(L2[ ,'baseline']),
    as.numeric(c0)
  )))

fit0 = v0 = cont0 = list()
for (i in 1:length(pb)) {
  
  # init data 
  meta = samplemd
  contrast_matrix = c0
  counts = pb[[i]]
  
  # dge list 
  d = edgeR::DGEList(counts = counts, samples = meta)
  
  # filter cell type specific lowly expressed genes 
  gtable = edgeR::filterByExpr(y = d$counts, 
                               min.count = 3,
                               design = as.factor(d$samples$irae.time)
                               )
  table(gtable)
  d = d[gtable, keep.lib.sizes=FALSE]
  d = edgeR::calcNormFactors(object = d)
  
  # get voom observation level weights
  v = voom(counts = d,  design = design.matrix, save.plot = TRUE, plot = TRUE)
  
  # fit contrast mixed model 
  fit = limma::lmFit(object = v, design = design.matrix)
  cfit = contrasts.fit(fit = fit, contrasts = contrast_matrix)
  eb = limma::eBayes(fit = cfit)
  
  # store results 
  v0[[i]] = v
  fit0[[i]] = fit
  cont0[[i]] = eb
}
# name the elements of result lists 
names(v0) = names(cont0) = names(pb)


# save obj. used in model fits  
saveRDS(object = samplemd, file = paste0(datapath, 'samplemd.rds'))
saveRDS(object = f1, file = paste0(datapath, 'f1.rds'))
saveRDS(object = L2, file = paste0(datapath, 'L2.rds'))
saveRDS(object = c0, file = paste0(datapath, 'c0.rds'))

# save model fits
saveRDS(v1, file = paste0(datapath, 'v1.rds'))
saveRDS(fit1, file = paste0(datapath, 'fit1.rds'))
saveRDS(v0, file = paste0(datapath, 'v0.rds'))
saveRDS(cont0, file = paste0(datapath, 'cont0.rds'))

sessionInfo()
#######################################
# R version 4.0.5 Patched (2021-03-31 r80136)
# 
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] scglmmr_0.1.0            variancePartition_1.25.6 BiocParallel_1.24.1      limma_3.46.0            
# [5] magrittr_2.0.1           here_1.0.1               forcats_0.5.1            stringr_1.4.0           
# [9] dplyr_1.0.4              purrr_0.3.4              readr_1.4.0              tidyr_1.1.2             
# [13] tibble_3.0.6             ggplot2_3.3.3            tidyverse_1.3.0         
# 
# loaded via a namespace (and not attached):
# [1] tidyselect_1.1.0            lme4_1.1-26                 htmlwidgets_1.5.3           RSQLite_2.2.7              
# [5] AnnotationDbi_1.52.0        grid_4.0.5                  scatterpie_0.1.7            munsell_0.5.0              
# [9] codetools_0.2-18            statmod_1.4.35              withr_2.4.1                 colorspace_2.0-0           
# [13] GOSemSim_2.16.1             Biobase_2.50.0              knitr_1.31                  rstudioapi_0.13            
# [17] stats4_4.0.5                ggsignif_0.6.0              DOSE_3.16.0                 labeling_0.4.2             
# [21] MatrixGenerics_1.2.1        Rdpack_2.1.1                emmeans_1.5.4               GenomeInfoDbData_1.2.4     
# [25] polyclip_1.10-0             pheatmap_1.0.12             bit64_4.0.5                 farver_2.0.3               
# [29] rprojroot_2.0.2             downloader_0.4              coda_0.19-4                 vctrs_0.3.6                
# [33] generics_0.1.0              TH.data_1.0-10              xfun_0.24                   R6_2.5.0                   
# [37] doParallel_1.0.16           GenomeInfoDb_1.26.7         graphlayouts_0.7.2          locfit_1.5-9.4             
# [41] bitops_1.0-6                cachem_1.0.4                fgsea_1.16.0                DelayedArray_0.16.3        
# [45] assertthat_0.2.1            scales_1.1.1                nnet_7.3-15                 multcomp_1.4-16            
# [49] ggraph_2.0.5                enrichplot_1.10.2           gtable_0.3.0                egg_0.4.5                  
# [53] tidygraph_1.2.0             sandwich_3.0-0              rlang_0.4.10                slanter_0.2-0              
# [57] splines_4.0.5               rstatix_0.7.0               checkmate_2.0.0             broom_0.7.5                
# [61] abind_1.4-5                 BiocManager_1.30.10         reshape2_1.4.4              modelr_0.1.8               
# [65] backports_1.2.1             Hmisc_4.5-0                 qvalue_2.22.0               clusterProfiler_3.18.1     
# [69] tools_4.0.5                 ellipsis_0.3.1              gplots_3.1.1                RColorBrewer_1.1-2         
# [73] BiocGenerics_0.36.1         Rcpp_1.0.6                  plyr_1.8.6                  base64enc_0.1-3            
# [77] progress_1.2.2              zlibbioc_1.36.0             RCurl_1.98-1.3              prettyunits_1.1.1          
# [81] rpart_4.1-15                ggpubr_0.4.0                viridis_0.5.1               cowplot_1.1.1              
# [85] S4Vectors_0.28.1            zoo_1.8-8                   cluster_2.1.2               SummarizedExperiment_1.20.0
# [89] haven_2.3.1                 ggrepel_0.9.1               fs_1.5.0                    data.table_1.14.0          
# [93] lmerTest_3.1-3              DO.db_2.9                   openxlsx_4.2.3              reprex_1.0.0               
# [97] mvtnorm_1.1-1               matrixStats_0.58.0          hms_1.0.0                   GSVA_1.38.2                
# [101] xtable_1.8-4                pbkrtest_0.5-0.1            RhpcBLASctl_0.21-247.1      XML_3.99-0.6               
# [105] jpeg_0.1-8.1                rio_0.5.16                  readxl_1.3.1                IRanges_2.24.1             
# [109] gridExtra_2.3               compiler_4.0.5              KernSmooth_2.23-18          crayon_1.4.1               
# [113] shadowtext_0.0.9            htmltools_0.5.1.1           minqa_1.2.4                 ggfun_0.0.4                
# [117] Formula_1.2-4               snow_0.4-3                  lubridate_1.7.9.2           DBI_1.1.1                  
# [121] corrplot_0.84               tweenr_1.0.2                dbplyr_2.1.0                MASS_7.3-53.1              
# [125] boot_1.3-27                 Matrix_1.3-2                car_3.0-10                  cli_2.5.0                  
# [129] rbibutils_2.0               parallel_4.0.5              igraph_1.2.6                GenomicRanges_1.42.0       
# [133] pkgconfig_2.0.3             rvcheck_0.1.8               numDeriv_2016.8-1.1         foreign_0.8-81             
# [137] xml2_1.3.2                  foreach_1.5.1               annotate_1.68.0             XVector_0.30.0             
# [141] GeneOverlap_1.26.0          estimability_1.3            rvest_0.3.6                 digest_0.6.27              
# [145] graph_1.68.0                cellranger_1.1.0            fastmatch_1.1-0             htmlTable_2.1.0            
# [149] edgeR_3.32.1                GSEABase_1.52.1             curl_4.3                    gtools_3.8.2               
# [153] nloptr_1.2.2.2              lifecycle_1.0.0             nlme_3.1-152                jsonlite_1.7.2             
# [157] aod_1.3.1                   carData_3.0-4               viridisLite_0.3.0           pillar_1.4.7               
# [161] lattice_0.20-41             fastmap_1.1.0               httr_1.4.2                  survival_3.2-10            
# [165] GO.db_3.12.1                glue_1.4.2                  zip_2.1.1                   png_0.1-7                  
# [169] iterators_1.0.13            bit_4.0.4                   ggforce_0.3.3               stringi_1.5.3              
# [173] blob_1.2.1                  org.Hs.eg.db_3.12.0         latticeExtra_0.6-29         caTools_1.18.1             
# [177] memoise_2.0.0         

