# add metadata to object 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(Seurat))

# set save path in data dir. 
datapath = file.path(here("data/processed_data/")); dir.create(datapath)

# read sample-level Metadata table
meta_table = read_delim(file = here("data/starting_data/subject_data/meta_table.txt"), delim = "\t")
meta_table = meta_table %>% 
  rename(htosample = 'sample') %>% 
  mutate(cohort_timepoint = paste(irae, timepoint, sep = "_")) %>% 
  mutate(sample = paste(sampleid, timepoint, sep = "_"))

# read clustered single cell data without sample and cluster annotations.   
s = readRDS(file = "data/starting_data/thymoma/ThymomaSeurat.20200422.clustered.Rds") 
s$htosample = s$sample
s$sample = NULL

# convert factors 
s@meta.data = 
  s@meta.data %>% 
  rownames_to_column("bc") %>% 
  mutate_if(is.factor, as.character) %>% 
  column_to_rownames("bc")
# save remapped metadata 
smet = s@meta.data

# add metadata based on hash id 
ids_map = names(meta_table)
md = s@meta.data
meta = list()
for (i in 1:length(ids_map)) {
  md[[ids_map[i]]] = plyr::mapvalues(md$htosample, from = meta_table$htosample, to = meta_table[[ids_map[i]]] )
  meta[[i]] = md %>% select(ids_map[i])
} 
meta = do.call(cbind, meta) 
meta = meta[ ,-1]
mdnew = cbind(meta, smet)

# smetadata to add to object 
saveRDS(meta, file = paste0(datapath, "metadata_add.rds"))
saveRDS(mdnew, file = paste0(datapath, "full_remapped_metadata_thymoma.rds"))

# add reformatted metadata ensure matched cols of assays with rows of md (automatic with AddMetaData)
stopifnot(isTRUE(all.equal(rownames(mdnew), colnames(s@assays$RNA@data) )))
s@meta.data[colnames(s@meta.data)] = NULL
s = AddMetaData(s,metadata = mdnew)

# save 
saveRDS(s,file = paste0(datapath, 'ThymomaSeurat.20200422.clustered.mdprocessed.rds'))


sessionInfo()
# R version 3.6.1 (2019-07-05)
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] Seurat_3.1.5    here_0.1        forcats_0.5.0   stringr_1.4.0   dplyr_1.0.3     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_3.1.0   
# [10] ggplot2_3.3.0   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] tsne_0.1-3          nlme_3.1-145        bitops_1.0-6        fs_1.3.2            lubridate_1.7.4     RcppAnnoy_0.0.16    RColorBrewer_1.1-2 
# [8] httr_1.4.1          rprojroot_1.3-2     sctransform_0.2.1   tools_3.6.1         backports_1.2.1     utf8_1.1.4          R6_2.4.1           
# [15] irlba_2.3.3         KernSmooth_2.23-16  uwot_0.1.8          lazyeval_0.2.2      DBI_1.1.0           colorspace_1.4-1    npsurv_0.4-0       
# [22] withr_2.4.0         gridExtra_2.3       tidyselect_1.1.0    compiler_3.6.1      cli_2.3.1           rvest_0.3.5         xml2_1.3.2         
# [29] plotly_4.9.2        caTools_1.18.0      scales_1.1.0        lmtest_0.9-37       pbapply_1.4-2       ggridges_0.5.2      digest_0.6.27      
# [36] htmltools_0.4.0     pkgconfig_2.0.3     dbplyr_1.4.2        htmlwidgets_1.5.1   rlang_0.4.10        readxl_1.3.1        rstudioapi_0.11    
# [43] generics_0.0.2      zoo_1.8-7           jsonlite_1.6.1      ica_1.0-2           gtools_3.8.1        magrittr_2.0.1      patchwork_1.0.0    
# [50] Matrix_1.2-18       Rcpp_1.0.4          munsell_0.5.0       fansi_0.4.2         ape_5.3             reticulate_1.14     lifecycle_1.0.0    
# [57] stringi_1.4.6       MASS_7.3-51.5       gplots_3.0.3        Rtsne_0.15          plyr_1.8.6          grid_3.6.1          parallel_3.6.1     
# [64] gdata_2.18.0        listenv_0.8.0       ggrepel_0.8.2       crayon_1.4.1        lattice_0.20-40     haven_2.2.0         cowplot_1.0.0      
# [71] splines_3.6.1       hms_0.5.3           pillar_1.5.0        igraph_1.2.5        reshape2_1.4.3      future.apply_1.4.0  codetools_0.2-16   
# [78] leiden_0.3.3        reprex_0.3.0        glue_1.4.2          lsei_1.2-0          data.table_1.12.8   modelr_0.1.6        png_0.1-7          
# [85] vctrs_0.3.6         cellranger_1.1.0    gtable_0.3.0        RANN_2.6.1          future_1.16.0       assertthat_0.2.1    rsvd_1.0.3         
# [92] broom_0.7.4         viridisLite_0.3.0   survival_3.1-11     cluster_2.1.0       globals_0.12.5      fitdistrplus_1.0-14 ellipsis_0.3.1     
# [99] ROCR_1.0-7   

