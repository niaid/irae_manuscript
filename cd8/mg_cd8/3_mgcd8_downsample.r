### read bulk results and data and ddownsampling 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))

figpath = here("cd8/mg_cd8/figuresV4/")
datapath = here("cd8/mg_cd8/generated_dataV4/")



# full gene signatures enriched in temra baseline. 
g0.list.combined = readRDS(file = here('DE_V4/generated_data/temp/g0.list.combined.rds'))
sig = g0.list.combined$CD8Tcell_TEMRA %>% filter(NES > 0 & padj < 0.01)
csig = readRDS('signature_curation/combined_signatures.rds')
hlmk = readRDS('signature_curation/hallmark.rds')
csig = c(csig, hlmk)
s.test = csig[sig$pathway]

# mtor sigs 
mtor.sigs2 = readRDS(file = here('cd8/mg_cd8/generated_dataV4/mtor.sigs2.rds'))
s.test = c(s.test, mtor.sigs2)

# gated temra model fit - get ranks for contrasts being tested
# run gene set enrichment on mtor list 
eb = readRDS(file = here('cd8/mg_cd8/generated_dataV4/eb.rds'))
r.irae.noirae = ExtractResult(list('gated' = eb), what = 'gene.t.ranks',coefficient.number = 2, coef.name = 'irae_vs_no_irae')
r.irae.healthy = ExtractResult(list('gated' = eb), what = 'gene.t.ranks', coefficient.number = 3, coef.name = 'irae_vs_healthy')
g.irae.healthy = fgsea::fgsea(pathways = s.test, stats = r.irae.healthy$gated, nperm = 2.5e5)
g.irae.noirae = fgsea::fgsea(pathways = s.test, stats = r.irae.noirae$gated, nperm = 2.5e5)
# save 
saveRDS(g.irae.healthy, file = paste0(datapath,'g.irae.healthy.rds'))
saveRDS(g.irae.noirae, file = paste0(datapath,'g.irae.noirae.rds'))



# downsample cells to get confidence for mtor sig at baseline from the gated population 
Downsample = function(x, ncell){
  if (length(x) < ncell) {
    x = x
  } else{
    x = sample(x, size = ncell, replace = FALSE)
  }
}

# select reasonable number to downsample from thymoma samples based on 
# the number of cells from healthy donors 
# vector of cell ids from each donor from script 1 
# { scell = lapply(X = split(meta, f = meta$sampleid), FUN = rownames) }
scell = readRDS('cd8/mg_cd8/generated_dataV4/scell.rds')
ncells = unlist(lapply(scell, length))
med.healthy = median(ncells[grep('health', names(ncells))]) # 45 cells 


# downsample gated temra cells from each donor  
subcells = list()
for (i in 1:100) {
  subcells[[i]] = lapply(scell, Downsample , ncell = med.healthy)
}
saveRDS(subcells, file = paste0(datapath,'subcells.rds'))


# create a pseudobulk matrix for each downsampled cell set
# load umi matrix from gated cells 
umi = readRDS(file = here('cd8/mg_cd8/generated_dataV4/umi.rds'))
pb2 = list()
for (i in 1:length(subcells)) {
  pbulk_list = lapply(subcells[[i]], function(x) Matrix::rowSums(umi[ ,x]))
  pb2[[i]] = as.matrix(do.call(cbind, pbulk_list))
}
names(pb2) = paste('subsample', 1:100)


# specify contrasts
design = readRDS(file = here('cd8/mg_cd8/generated_data/design.rds'))
c_mat =  limma::makeContrasts(
  irae_vs_no_irae = (IRAE - NO_IRAE), 
  irae_vs_healthy = (IRAE  - healthy), 
  levels = colnames(design)
)


# pseudobulk pipeline for the list of matrices made from downsampling 
lmFitSamples = function(pbulk_list, genes.fit = 'filter', design){
  # format input data 
  # -- iterate over list of matrices 
  dat_return = lapply(pbulk_list, function(x){ 
    dat = x %>% 
      edgeR::DGEList() %>% 
      edgeR::calcNormFactors(method = "RLE")
    # fitler genes tested in fitted models 
    # either test a specific set of genes or filter 
    # using the edgeR criterion in filterByExpr based on design and min threshold  
    if (genes.fit == 'filter'){
      genes.use1 = edgeR::filterByExpr(dat, min.count = 1, design = design)
      dat = dat[genes.use1, keep.lib.sizes=FALSE]  
      }else { 
        dat = dat[genes.fit, keep.lib.sizes=FALSE]
      }
      # run model 
    fit = limma::voom(
      counts = dat, 
      design = design,
      normalize.method = "none",
      save.plot = FALSE, plot = FALSE) %>% 
      limma::lmFit(design = design) %>% 
      limma::contrasts.fit(contrasts = c_mat) %>% 
      limma::eBayes()
    return(fit)
    })
  return(dat_return)
}

# load genes tested in temra cluster 
temra.g.test = readRDS('cd8/mg_cd8/generated_data/temra.g.test')

# run 2 versions of downsampling analysis 
# one using fitted genes from full temra cluster 
# another using hte genes that meet the filtering criteria  
dge_ = lmFitSamples(pbulk_list = pb2, genes.fit = temra.g.test, design = design)
dge.filter = lmFitSamples(pbulk_list = pb2, genes.fit = 'filter', design = design)


# get ranks for filtered analysis and analysis using the itted genes from the TEMRA cluster. 
r.irae.filtered = ExtractResult(dge.filter, what = 'gene.t.ranks', coefficient.number = 1, coef.name = 'irae_vs_no_irae')
r.irae.subsample = ExtractResult(dge_, what = 'gene.t.ranks', coefficient.number = 1, coef.name = 'irae_vs_no_irae')
r.healthy.subsample = ExtractResult(dge_, what = 'gene.t.ranks', coefficient.number = 2, coef.name = 'irae_vs_healthy')
saveRDS(r.irae.filtered, file = paste0(datapath, 'r.irae.filtered.rds'))
saveRDS(r.irae.subsample, file = paste0(datapath, 'r.irae.subsample.rds'))
saveRDS(r.healthy.subsample, file = paste0(datapath, 'r.healthy.subsample.rds'))


# run gsea 
irae_gsea_filter = FgseaList(rank.list.celltype = r.irae.filtered, pathways = s.test, nperm = 2.5e5)
irae_gsea_samples = FgseaList(rank.list.celltype = r.irae.subsample, pathways = s.test, nperm = 2.5e5)
healthy_gsea_samples = FgseaList(rank.list.celltype = r.healthy.subsample, pathways = s.test, nperm = 2.5e5)

# save 
saveRDS(irae_gsea_filter, file = paste0(datapath,'irae_gsea_filter.rds'))
saveRDS(irae_gsea_samples, file = paste0(datapath,'irae_gsea_samples.rds'))
saveRDS(healthy_gsea_samples, file = paste0(datapath,'healthy_gsea_samples.rds'))


sessionInfo()
# R version 3.6.1 (2019-07-05)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] scglmmr_0.1.0   here_0.1        forcats_0.5.0   stringr_1.4.0   dplyr_1.0.3    
# [6] purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_3.1.0    ggplot2_3.3.0  
# [11] tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1             backports_1.2.1          fastmatch_1.1-0         
# [4] plyr_1.8.6               igraph_1.2.5             GSEABase_1.48.0         
# [7] splines_3.6.1            BiocParallel_1.20.1      TH.data_1.0-10          
# [10] urltools_1.7.3           digest_0.6.27            foreach_1.4.8           
# [13] htmltools_0.4.0          GOSemSim_2.12.1          viridis_0.5.1           
# [16] GO.db_3.10.0             gdata_2.18.0             fansi_0.4.2             
# [19] magrittr_2.0.2           memoise_1.1.0            doParallel_1.0.15       
# [22] limma_3.42.2             annotate_1.64.0          graphlayouts_0.7.0      
# [25] modelr_0.1.6             GeneOverlap_1.22.0       sandwich_2.5-1          
# [28] enrichplot_1.6.1         prettyunits_1.1.1        colorspace_1.4-1        
# [31] blob_1.2.1               rvest_0.3.5              ggrepel_0.8.2           
# [34] haven_2.2.0              crayon_1.4.1             RCurl_1.98-1.1          
# [37] jsonlite_1.6.1           graph_1.64.0             lme4_1.1-21             
# [40] iterators_1.0.12         survival_3.1-11          zoo_1.8-7               
# [43] glue_1.4.2               slanter_0.2-0            polyclip_1.10-0         
# [46] gtable_0.3.0             emmeans_1.4.5            BiocGenerics_0.32.0     
# [49] scales_1.1.0             DOSE_3.12.0              pheatmap_1.0.12         
# [52] mvtnorm_1.1-0            edgeR_3.28.1             DBI_1.1.0               
# [55] Rcpp_1.0.4               viridisLite_0.3.0        xtable_1.8-4            
# [58] progress_1.2.2           gridGraphics_0.5-0       bit_1.1-15.2            
# [61] europepmc_0.3            stats4_3.6.1             GSVA_1.34.0             
# [64] httr_1.4.1               fgsea_1.12.0             gplots_3.0.3            
# [67] RColorBrewer_1.1-2       ellipsis_0.3.1           pkgconfig_2.0.3         
# [70] XML_3.99-0.3             farver_2.0.3             dbplyr_1.4.2            
# [73] locfit_1.5-9.4           utf8_1.1.4               ggplotify_0.0.5         
# [76] tidyselect_1.1.0         rlang_0.4.10             reshape2_1.4.3          
# [79] later_1.0.0              AnnotationDbi_1.48.0     munsell_0.5.0           
# [82] cellranger_1.1.0         tools_3.6.1              cli_2.3.1               
# [85] generics_0.0.2           RSQLite_2.2.0            ggridges_0.5.2          
# [88] broom_0.7.4              fastmap_1.0.1            org.Hs.eg.db_3.10.0     
# [91] bit64_0.9-7              fs_1.3.2                 tidygraph_1.1.2         
# [94] caTools_1.18.0           ggraph_2.0.2             egg_0.4.5               
# [97] nlme_3.1-145             mime_0.9                 DO.db_2.9               
# [100] xml2_1.3.2               pbkrtest_0.4-8.6         compiler_3.6.1          
# [103] shinythemes_1.1.2        rstudioapi_0.11          testthat_2.3.2          
# [106] variancePartition_1.16.1 ggsignif_0.6.0           reprex_0.3.0            
# [109] tweenr_1.0.1             geneplotter_1.64.0       stringi_1.4.6           
# [112] desc_1.2.0               lattice_0.20-40          Matrix_1.2-18           
# [115] nloptr_1.2.2.1           ggsci_2.9                vctrs_0.3.6             
# [118] pillar_1.5.0             lifecycle_1.0.0          BiocManager_1.30.10     
# [121] triebeard_0.3.0          estimability_1.3         data.table_1.12.8       
# [124] cowplot_1.0.0            bitops_1.0-6             colorRamps_2.3          
# [127] httpuv_1.5.2             qvalue_2.18.0            R6_2.4.1                
# [130] promises_1.1.0           KernSmooth_2.23-16       gridExtra_2.3           
# [133] IRanges_2.20.2           codetools_0.2-16         pkgload_1.0.2           
# [136] boot_1.3-24              MASS_7.3-51.5            gtools_3.8.1            
# [139] assertthat_0.2.1         rprojroot_1.3-2          withr_2.4.0             
# [142] multcomp_1.4-12          S4Vectors_0.24.3         parallel_3.6.1          
# [145] hms_0.5.3                clusterProfiler_3.14.3   grid_3.6.1              
# [148] coda_0.19-4              minqa_1.2.4              rvcheck_0.1.8           
# [151] ggpubr_0.2.5             ggforce_0.3.1            Biobase_2.46.0          
# [154] shiny_1.4.0.2            lubridate_1.7.4    