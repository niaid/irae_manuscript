Contrasting autoimmune and treatment effects reveals baseline set points
of immune toxicity following checkpoint inhibition
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

### Code to reproduce all manuscript results and figures

#### Matt Mulè

biorxiv preprint: 
https://www.biorxiv.org/content/10.1101/2022.06.05.494592v1

## Table of Contents

1.  [**Instructions for analysis workflow**](#instructions)
2.  [**Data processing**](#data%20processing)
3.  [**Assemble *a priori hypothesis* gene signatures for independent
    testing**](#modules)  
4.  [**Aggregated mixed effects models contrasting mRNA
    states**](#aggregated)
5.  [**Manual gating healthy donor and thymic cancer cells and
    robustness assessment of irAE signatures**](#mg)
6.  [**Single cell mixed effects modeling of T cell
    states**](#singlecell)
7.  [**Integration**](#integration)
8.  [**Additional analysis**](#misc)
9.  [**Analysis of tissue T cell data**](#luoma)

#### Instructions for analysis workflow <a name="instructions"></a>

This analysis was done using **R version 3.6.1**

The main packages used for this analysis and their versions are listed
below:  
**variancePartition** see below
[**scglmmr**](https://github.com/MattPM/scglmmr) version 0.1.0  
[**dsb**](https://github.com/niaid/dsb) version 0.1.0  
**lme4** version 1.1-21  
**limma** version 3.42.2  
**Seurat** version 3.1.5  
**tidyverse** version 1.3.0  
**here** version 0.1  
**dplyr** version 1.0.3  
**emmeans** version 1.4.5

To run the analysis, 1) download the repository 2) download starting
data and add the /data folder directly to the root directory containing
the .Rproj file; this directory should now contain irae.Rproj, the files
readme.md and readme.rmd, the directories containing analysis scripts
the directory you just added, data.

One can view the commented code and run each script in each subdirectory
as listed below, or source each R script in the order they appear below.
No file paths need to be specified or changed. Each R script is
self-contained, reading data from the /data folder and writing to
figures or results files within each analysis subdirectory relative to
the root directory using the R package `here`.

To install exact versions of packages used in the analysis run code
below.

``` r
# 
require(devtools)

# install variancePartition release 3 14 
devtools::install_github(repo = "https://github.com/GabrielHoffman/variancePartition@RELEASE_3_14")

# install scglmmr 
devtools::install_github(repo = "https://github.com/MattPM/scglmmr")

# install 
pkgs3.6 = c('dsb', 'lme4', 'limma', 'Seurat', 'tidyverse', 'here', 'dplyr', 'emmeans')
pkg.version = c('0.1.0', '1.1-21', '3.42.2', '3.1.5', '1.3.0', '0.1', '1.0.3', '1.4.5')

for (i in 1:length(pkgs3.6)) {
  devtools::install_version(package = pkgs3.6[[i]], version = pkg.version[i]) 
}
```

<br/> <br/>

#### Data processing <a name="data processing"></a>

Directory : **data**  
Input to this script downloaded from the data repository:  
1) “data/starting_data/subject_data/meta_table.txt”  
2) “data/starting_data/ThymomaSeurat.20200422.clustered.Rds”  
Output:  
“data/processed_data/ThymomaSeurat.20200422.clustered.mdprocessed.rds”  
**ThymomaSeurat.20200422.clustered.mdprocessed.rds** object is used
throughout analysis.

``` r
source(here('data/1_add_metadata_v3.r'))
```

<br/> <br/>

### Assemble *a priori hypothesis* gene signatures for independent testing <a name="modules"></a>

This compendium of gene signatures we pre selected to test is included.
The code below reproduces construction of the gene signature.  
Directory : **signature_curation**

``` r
source(here('signature_curation/2_create_gene_set_testing_V2.r'))
```

Subsets hypothesis gene sets from reactome, the Blood TRanscriptional
Modules (BTM) Li et al -
<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3946932/>  
Additional gene signatures included from:  
<https://pubmed.ncbi.nlm.nih.gov/31209400/>  
<https://pubmed.ncbi.nlm.nih.gov/23521917/>  
Directories with starting data: **Kotliarov_2020**  
**Shahabi_2013**  
**yao_ni_2019**

<br/> <br/>

#### Aggregated mixed effects models contrasting mRNA states <a name="aggregated"></a>

Directory: **DE_V4**

1_create_bulk_object_list.r - pool mRNA from each protein-defined subset
into pseudobulk library.  
2_fit_mixed_model_contrast_v4.r - fit models including the three main
contrasts (baseline, treatment effect, difference in treatment effects
between groups).  
3_contrast_gsea_hypset_hlmk - run gsea for each effect across cell types
on the hypothesis set of gene modules and the msigdb hallmark pathways.
Genes are ranked by t statistic from the modeled contrasts fit in script
2.  
4_mtor_heatmap_temra.r - generate figures  
4_partII_contrast_visualizations.r - generate figures  
5_baseline_gsea_visualization.r - generate figures  
6_sub_sig_venndiagram.r - generate figures

``` r
source(here('DE_V4/1_create_bulk_object_list.r'))
source(here('DE_V4/2_fit_mixed_model_contrast_v4.R'))
source(here('DE_V4/3_contrast_gsea_hypset_hlmk.r'))
source(here('DE_V4/4_mtor_heatmap_temra.r'))
source(here('DE_V4/4_partII_contrast_visualizations.r'))
source(here('DE_V4/5_baseline_gsea_visualization.r'))
source(here('DE_V4/6_sub_sig_venndiagram.r'))
source(here('DE_V4/Final_script_write.table.r'))
```

<br/> <br/>

### CD8 T cell analysis

#### Manual gating healthy donor and thymic cancer cells robustness assessment of signatures <a name="mg"></a>

1_mg_cd8.r - manually gate cells  
2_mg.cd8.avexprs.r - compare average signature expression in manually
gated cells.  
3_mgcd8_downsample.r - assess robustness of enrichment using k cell
downsampling and permutation  
4_mgcd8.vis.r - generate figures  
<br/>

``` r
source(here('cd8/mg_cd8/1_mg_cd8.r'))
source(here('cd8/mg_cd8/2_mg.cd8.avexprs.r'))
source(here('cd8/mg_cd8/3_mgcd8_downsample.r'))
source(here('cd8/mg_cd8/4_mgcd8.vis.r'))
```

<br/>

#### Single cell mixed effects modeling of T cell states <a name="singlecell"></a>

Directory: **mtor_cd8**

1_mtor_modules_scglmmr.r - score module activity and fit mixed effects
model within cell subets at single cell level.  
2_mixed_model_estimates_vis.r - generate figures (effect size
comparisons; baseline and treatment effects with CI)  
3_mtor_sig_venn.r - generate figures

``` r
source(here('cd8/mtor_cd8/1_mtor_modules_scglmmr.r'))
source(here('cd8/mtor_cd8/2_mixed_model_estimates_vis.r'))
source(here('cd8/mtor_cd8/3_mtor_sig_venn.r'))
```

<br/>

#### Integration with healthy T cells <a name="integration"></a>

Directory **healthy_integration**

healthyintegration_subjectregression.r - integrate thymic cancer cells
with healthy donor cells.
healthy_integration_subjectregression_figure_generation.r - generate
figures

``` r
source(here('cd8/healthy_integration/healthyintegration_subjectregression.r'))
source(here('cd8/healthy_integration/healthy_integration_subjectregression_figure_generation.r'))
```

<br/> <br/>

#### Additional analysis <a name="misc"></a>

Directory: **misc**  
Longitudinal clinical labs (CK) IRAE and CD8 T cell dsb normalized
protein visualizations, and umap / sample level stacked bar plots for
fig 1 visualization of dataset.  
clinical_labs_plot.r  
CD8subsets_proteinexprs.r  
umap_barplot_thymoma.r

``` r
source(here('misc/CD8subsets_proteinexprs.r'))
source(here('misc/clinical_labs_plot.r'))
source(here('misc/proteinmap_umap_samplebar.r'))
```

#### Analysis of tissue T cell data <a name="luoma"></a>

Analysos of data from [Luoma et.
al. 2020](https://pubmed.ncbi.nlm.nih.gov/32603654/)

``` r
source(here('louma/1_process_GSE144469.R'))
source(here('louma/2_luoma_pb_workflow.r'))
source(here('louma/4_cd3_cluster_mtor.r'))
source(here('louma/luoma_clustassociation.r'))
source(here('louma/luoma_cd8_genemap.r'))
source(here('louma/5_cd3_figuregeneration.r'))
```

<br/> <br/>

Some utility functions are in the directory `util` other custom
functions are included in the scglmmr package.

#### NIAID repository release notes

A review of this code has been conducted, no critical errors exist, and
to the best of the authors knowledge, there are no problematic file
paths, no local system configuration details, and no passwords or keys
included in this code.  
Primary author(s): Matt Mulè  
Organizational contact information: General: john.tsang AT nih.gov,
code: mulemp AT nih.gov  
Description: code to reproduce analysis of manuscript  
Usage instructions: Provided in this markdown

``` r
sessionInfo()


# R version 3.6.1 (2019-07-05)  
# Matrix products: default  
# 
# attached base packages:  
# [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#  [1] variancePartition_1.16.1 Biobase_2.46.0           BiocGenerics_0.32.0      scales_1.1.0            
#  [5] foreach_1.4.8            limma_3.42.2             magrittr_2.0.2           here_0.1                
#  [9] emmeans_1.4.5            lme4_1.1-21              Matrix_1.2-18            forcats_0.5.0           
# [13] stringr_1.4.0            dplyr_1.0.3              purrr_0.3.3              readr_1.3.1             
# [17] tidyr_1.0.2              tibble_3.1.0             ggplot2_3.3.0            tidyverse_1.3.0         
# [21] Seurat_3.1.5            
# 
# loaded via a namespace (and not attached):
#   [1] readxl_1.3.1        backports_1.2.1     plyr_1.8.6          igraph_1.2.5        lazyeval_0.2.2     
#   [6] splines_3.6.1       BiocParallel_1.20.1 listenv_0.8.0       usethis_1.5.1       TH.data_1.0-10     
#  [11] digest_0.6.27       htmltools_0.4.0     gdata_2.18.0        fansi_0.4.2         memoise_1.1.0      
#  [16] cluster_2.1.0       doParallel_1.0.15   ROCR_1.0-7          remotes_2.1.1       globals_0.12.5     
#  [21] modelr_0.1.6        sandwich_2.5-1      prettyunits_1.1.1   colorspace_1.4-1    rvest_0.3.5        
#  [26] ggrepel_0.8.2       haven_2.2.0         xfun_0.12           callr_3.4.2         crayon_1.4.1       
#  [31] jsonlite_1.6.1      survival_3.1-11     zoo_1.8-7           iterators_1.0.12    ape_5.3            
#  [36] glue_1.4.2          gtable_0.3.0        leiden_0.3.3        pkgbuild_1.0.6      future.apply_1.4.0 
#  [41] mvtnorm_1.1-0       DBI_1.1.0           Rcpp_1.0.4          viridisLite_0.3.0   xtable_1.8-4       
#  [46] progress_1.2.2      reticulate_1.14     rsvd_1.0.3          tsne_0.1-3          htmlwidgets_1.5.1  
#  [51] httr_1.4.1          gplots_3.0.3        RColorBrewer_1.1-2  ellipsis_0.3.1      ica_1.0-2          
#  [56] pkgconfig_2.0.3     uwot_0.1.8          dbplyr_1.4.2        utf8_1.1.4          tidyselect_1.1.0   
#  [61] rlang_0.4.10        reshape2_1.4.3      munsell_0.5.0       cellranger_1.1.0    tools_3.6.1        
#  [66] cli_2.3.1           generics_0.0.2      devtools_2.2.2      broom_0.7.4         ggridges_0.5.2     
#  [71] evaluate_0.14       yaml_2.2.1          npsurv_0.4-0        processx_3.4.2      knitr_1.28         
#  [76] fs_1.3.2            fitdistrplus_1.0-14 caTools_1.18.0      RANN_2.6.1          pbapply_1.4-2      
#  [81] future_1.16.0       nlme_3.1-145        xml2_1.3.2          compiler_3.6.1      pbkrtest_0.4-8.6   
#  [86] rstudioapi_0.11     curl_4.3            plotly_4.9.2        png_0.1-7           testthat_2.3.2     
#  [91] lsei_1.2-0          reprex_0.3.0        stringi_1.4.6       ps_1.3.2            desc_1.2.0         
#  [96] lattice_0.20-40     nloptr_1.2.2.1      vctrs_0.3.6         pillar_1.5.0        lifecycle_1.0.0    
# [101] lmtest_0.9-37       RcppAnnoy_0.0.16    estimability_1.3    data.table_1.12.8   cowplot_1.0.0      
# [106] bitops_1.0-6        irlba_2.3.3         patchwork_1.0.0     colorRamps_2.3      R6_2.4.1           
# [111] KernSmooth_2.23-16  gridExtra_2.3       sessioninfo_1.1.1   codetools_0.2-16    pkgload_1.0.2      
# [116] boot_1.3-24         MASS_7.3-51.5       gtools_3.8.1        assertthat_0.2.1    rprojroot_1.3-2    
# [121] withr_2.4.0         sctransform_0.2.1   multcomp_1.4-12     hms_0.5.3           grid_3.6.1         
# [126] coda_0.19-4         minqa_1.2.4         rmarkdown_2.1       Rtsne_0.15          lubridate_1.7.4  
```
