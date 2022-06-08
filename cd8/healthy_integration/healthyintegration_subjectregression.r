# healthy integration script
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(here))

# define paths 
figpath = here("cd8/healthy_integration/figures/"); dir.create(figpath, recursive = T)
datapath = here("cd8/healthy_integration/generated_data/"); dir.create(datapath, recursive = T)

## load dataand subset to CD8 T cells 
s = readRDS(here("data/sg_data-PROJECTS_-RAE-200210_A00956_0050_BHNTNTDSXX_processed-R_analysis/ThymomaSeurat_20200422_clustered_Meta_added_SNG_add_hypmod_scores.Rds"))
ct = unique(s@meta.data$celltype)
tc = ct[grep(x = ct, pattern = "CD8Tcel")] %>% 
  droplevels() %>% 
  as.vector()
Idents(s) = "celltype"
s = SubsetData(s, ident.use = tc)

# Read flu data and convert to Seurat V3 object. 
flu = readRDS(file = here("data/baseline_flu/H1_day0_demultilexed_singlets.RDS"))
singlet = flu@data %>% colnames
flu_count = flu@raw.data[ ,singlet] 
flumd = readRDS(here("data/baseline_flu/md.rds"))
flumd = flumd %>% 
  filter(barcode_check %in% singlet & cohort =="H1N1" & timepoint == "d0") %>% 
  select(batch, age, gender, sampleid, celltype_joint, adjmfc.group, sampleid) 
h1 = CreateSeuratObject(counts = flu_count, meta.data = flumd, min.cells = 20, min.features = 200)


# subset to CD8 T cells 
Idents(h1) = "celltype_joint" 
ctflu = unique(h1@meta.data$celltype_joint)
tcflu = ctflu[grep(x = ctflu, pattern = "CD8")]  %>% as.vector()
h1 = SubsetData(h1, ident.use = tcflu)

# create, noramalize and score modules in list of objects 
sl = list("flu" = h1, "thymoma" = s)
sl = lapply(sl,  function(x){  
  SCTransform(object = x, vars.to.regress = c('sampleid')) 
  })

# define genes to integrate on 
integration_genes = SelectIntegrationFeatures(object.list = sl, nfeatures = 2000)

# prep object 
sl = PrepSCTIntegration(object.list = sl, 
                        anchor.features = integration_genes)

# set healthy donors as reference dataset 
reference_dataset = which(names(sl) == "flu")

# find integration anchors
anchors = FindIntegrationAnchors(object.list = sl, 
                                 normalization.method = "SCT",
                                 anchor.features = integration_genes,
                                 reference = reference_dataset)

# integrate 
int = IntegrateData(anchorset = anchors, 
                    normalization.method = "SCT")

# run pca 
int = RunPCA(int)
int = RunUMAP(int, dims = 1:30)

# cluster integrated data 
DefaultAssay(int) <- "integrated"
int = FindNeighbors(int,reduction = 'pca', dims = 1:30, nn.eps = 0.5)
int = FindClusters(int, resolution = 0.6, n.start = 10)

# save object 
saveRDS(object = int,file = paste0(datapath, 'integrated_thymoma_healthy_subjectregression.rds'))

sessionInfo()
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-conda_cos6-linux-gnu (64-bit)
# Running under: Red Hat Enterprise Linux Server 7.9 (Maipo)
# 
# Matrix products: default
# BLAS:   /sysapps/cluster/software/Anaconda2/5.3.0/lib/libblas.so.3.8.0
# LAPACK: /sysapps/cluster/software/Anaconda2/5.3.0/lib/liblapack.so.3.8.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Seurat_3.2.3    forcats_0.5.1   stringr_1.4.0   dplyr_1.0.4    
# [5] purrr_0.3.4     readr_1.4.0     tidyr_1.1.2     tibble_3.0.6   
# [9] ggplot2_3.3.3   tidyverse_1.3.0 here_1.0.1     
# 
# loaded via a namespace (and not attached):
#   [1] Rtsne_0.15           colorspace_2.0-0     deldir_0.2-9        
# [4] ellipsis_0.3.1       ggridges_0.5.3       rprojroot_2.0.2     
# [7] fs_1.5.0             spatstat.data_1.7-0  rstudioapi_0.13     
# [10] leiden_0.3.7         listenv_0.8.0        ggrepel_0.9.1       
# [13] lubridate_1.7.9.2    xml2_1.3.2           codetools_0.2-18    
# [16] splines_3.6.1        polyclip_1.10-0      jsonlite_1.7.2      
# [19] broom_0.7.4          ica_1.0-2            cluster_2.1.0       
# [22] dbplyr_2.0.0         png_0.1-7            uwot_0.1.10         
# [25] sctransform_0.3.2    shiny_1.6.0          compiler_3.6.1      
# [28] httr_1.4.2           backports_1.2.1      assertthat_0.2.1    
# [31] Matrix_1.3-2         fastmap_1.1.0        lazyeval_0.2.2      
# [34] cli_2.5.0            later_1.1.0.1        htmltools_0.5.1.1   
# [37] tools_3.6.1          rsvd_1.0.3           igraph_1.2.6        
# [40] gtable_0.3.0         glue_1.4.2           reshape2_1.4.4      
# [43] RANN_2.6.1           spatstat_1.64-1      Rcpp_1.0.6          
# [46] scattermore_0.7      cellranger_1.1.0     vctrs_0.3.6         
# [49] nlme_3.1-151         lmtest_0.9-38        globals_0.14.0      
# [52] ps_1.5.0             rvest_0.3.6          mime_0.9            
# [55] miniUI_0.1.1.1       lifecycle_1.0.0      irlba_2.3.3         
# [58] goftest_1.2-2        future_1.21.0        MASS_7.3-53         
# [61] zoo_1.8-8            scales_1.1.1         spatstat.utils_2.0-0
# [64] hms_1.0.0            promises_1.1.1       parallel_3.6.1      
# [67] RColorBrewer_1.1-2   gridExtra_2.3        reticulate_1.18     
# [70] pbapply_1.4-3        rpart_4.1-15         stringi_1.5.3       
# [73] rlang_0.4.10         pkgconfig_2.0.3      matrixStats_0.58.0  
# [76] lattice_0.20-41      tensor_1.5           ROCR_1.0-11         
# [79] patchwork_1.1.1      htmlwidgets_1.5.3    cowplot_1.1.1       
# [82] tidyselect_1.1.0     parallelly_1.23.0    RcppAnnoy_0.0.18    
# [85] plyr_1.8.6           magrittr_2.0.1       R6_2.5.0            
# [88] generics_0.1.0       DBI_1.1.1            mgcv_1.8-33         
# [91] pillar_1.4.7         haven_2.3.1          withr_2.4.1         
# [94] fitdistrplus_1.1-3   abind_1.4-5          survival_3.2-7      
# [97] future.apply_1.7.0   modelr_0.1.8         crayon_1.4.1        
# [100] KernSmooth_2.23-18   plotly_4.9.3         grid_3.6.1          
# [103] readxl_1.3.1         data.table_1.13.6    reprex_1.0.0        
# [106] digest_0.6.27        xtable_1.8-4         httpuv_1.5.5        
# [109] munsell_0.5.0        viridisLite_0.3.0   

