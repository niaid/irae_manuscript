# CD8 t cell states 
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(here))
library(scglmmr)
figpath = here("cd8/healthy_integration/figures_v2/"); dir.create(figpath)
datapath = here("cd8/healthy_integration/generated_data/")

# plot reduction 
int = readRDS("cd8/healthy_integration/generated_data/integrated_thymoma_healthy_subjectregression.rds")

d = cbind(int@meta.data, int@reductions$umap@cell.embeddings)
d$study = ifelse(d$cohort_timepoint %in% c('0_0', "0_1", "1_0", "1_1"), "Thymoma", "Healthy")
d$class = paste(d$study, d$irae, sep = "_")
d$class = str_replace(d$class, pattern = "_NA",replacement = "")
d$Batch = ifelse(is.na(d$Batch), yes = "B0", no = d$Batch)
d = d %>% filter(UMAP_2 < 6)

centers = d %>% dplyr::group_by(integrated_snn_res.0.6) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

# add that data layer to the plot 
boxbox = list(theme_bw(), 
              theme(panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(), 
                    axis.ticks.y = element_blank(),axis.ticks.x = element_blank(),
                    axis.text.x = element_blank(), axis.text.y = element_blank())
              ) 
p1 = ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = integrated_snn_res.0.6)) +
  boxbox + 
  geom_point(shape = 20, size = 0.1, show.legend = FALSE) + 
  ggsci::scale_color_d3(palette = 'category20') + 
  geom_point(data = centers, mapping = aes(x = UMAP_1, y = UMAP_2), size = 0, alpha = 0, show.legend = FALSE) + 
  ggrepel::geom_text_repel(data = centers,size = 6, color = 'black', mapping = aes(label = integrated_snn_res.0.6), show.legend = FALSE) 
ggsave(p1,filename = paste0(figpath,"integrated_clusters.png"),width = 4.5, height = 4)

# cohort split 
p2 = ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = study)) +
  theme_bw() + 
  theme(strip.background = element_blank()) + 
  facet_wrap(~study) + 
  geom_point(size = 0.2, show.legend = FALSE) + 
  ggsci::scale_color_nejm(alpha = 0.2)
ggsave(p2,filename = paste0(figpath,"integrated_clusters_cohortsplit.png"),width = 8, height = 4)


# TEMRA 
p = ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = integrated_snn_res.0.6)) +
  geom_point(size = 0.2, show.legend = F) + 
  geom_density_2d(data = d %>% filter(celltype == "CD8Tcell_TEMRA"), color = 'black') + 
  boxbox + 
  ggsci::scale_color_d3(palette = 'category20')
p
ggsave(p,filename = paste0(figpath, "TEMRA_density.png"), width = 4.5, height = 4)

# by subject 
d$sample[is.na(d$sample)] = "healthy"
d$adverse_event = ifelse(d$irae == '1', "IRAE", "NO IRAE")
d$sample_irae = paste(d$adverse_event, d$sample, sep = "_")
p = ggplot(d, aes(x = UMAP_1, y = UMAP_2, color = integrated_snn_res.0.6)) +
  geom_point(size = 0.2, show.legend = F) + 
  boxbox + 
  facet_wrap(~sample_irae) + 
  theme(strip.background = element_blank())  + 
  ggsci::scale_color_d3(palette = 'category20')
ggsave(p,filename = paste0(figpath, "subject_facet_integrated_clusters.png"), width = 11, height = 10)


# calculate mod score of enriched mTOR genes 
mtor_pathways = readRDS(file = here('DE_V3/data_highres/MTOR_SIG_LIST.rds'))
cycle_list = list(s_phase = cc.genes$s.genes, g2_m_phase = cc.genes$g2m.genes)
mtor_pathways = c(mtor_pathways, cycle_list)

# calc module scores 
ms = scglmmr::WeightedCellModuleScore(gene_matrix = int@assays$RNA@data, 
                                        module_list = mtor_pathways, 
                                        threshold = 0, return_weighted = F)

# add modulse scores (remove the outliers from umap as above)
d = d[ ,!colnames(d) %in% colnames(ms)]
ms = ms[rownames(ms) %in% rownames(d), ]
d = cbind(d, ms)

# hist by mtor score 
p = ggplot(d, aes(x = mtor_thymoma, color = cohort_timepoint)) + geom_density() 
ggsave(p, filename = paste0(figpath, 'temramtor_score_across_cohorts.pdf'), width = 5 , height = 4)

# calculate mtor expression threshold
med_mtor = median(d$mtor_thymoma)
mad_mtor = mad(d$mtor_thymoma)
mtor_threshold = med_mtor + 3.5*(mad_mtor)
#
p = ggplot(d, aes(x = UMAP_1, y = UMAP_2)) +
  boxbox + 
  geom_point(size = 0.2, show.legend = F) + 
  geom_density_2d(data = d %>% filter(mtor_thymoma > mtor_threshold), color = 'red', size = 1.2) 
ggsave(p,filename = paste0(figpath, "Mtor_gt_3.5mad_median_density.png"), width = 4.5, height = 4)

# calculate cell cycle expression threshold (s phase)
med_s = median(d$s_phase)
mad_s = mad(d$s_phase)
s_threshold = med_s + 3.5*(mad_s)
p = ggplot(d, aes(x = UMAP_1, y = UMAP_2)) +
  boxbox + 
  geom_point(size = 0.2, show.legend = F) + 
  geom_density_2d(data = d %>% filter(s_phase > s_threshold), color = 'yellow', size = 1.2) 
ggsave(p,filename = paste0(figpath, "sphase_gt_3.5mad_median_density.png"), width = 4.5, height = 4)

# calculate cell cycle expression threshold (g2m phase)
med_g2m = median(d$g2_m_phase)
mad_g2m = mad(d$g2_m_phase)
g2m_threshold = med_g2m + 3.5*(mad_g2m)
p = ggplot(d, aes(x = UMAP_1, y = UMAP_2)) +
  boxbox + 
  geom_point(size = 0.2, show.legend = F) + 
  geom_density_2d(data = d %>% filter(g2_m_phase > g2m_threshold), size = 1.2, color = 'purple') 
ggsave(p,filename = paste0(figpath, "g2mphase_gt_3.5mad_median_density.png"), width = 4.5, height = 4)

# Model the mtor score by integrated cluster 
d$cluster = d$integrated_snn_res.0.6
d$scaled_age = scale(as.numeric(d$age))
d$timepoint[(is.na(d$timepoint))] <- '0'
bl = d %>% filter(timepoint == 0)

f1 = as.formula(mtor_thymoma ~ cluster + scaled_age + Batch + (1|sampleid))
m1 = lme4::lmer(f1, data = bl)
emm1 = emmeans::emmeans(m1, specs = ~cluster:Batch)

# plot the marginal mean distribution over levels of cluster and batch. 
p = plot(emm1) + 
  theme_bw() + 
  ggtitle('mtor score model (baseline)') + 
  xlab('marginal means') + ylab('cluster:batch B0 = healthy donor B1 & B2 = Thymoma') + 
  geom_hline(yintercept = c(14,28)) + 
  geom_segment(aes(x = 0.1689, y = 0, xend = 0.1689, yend = 14), linetype = 'dashed') + 
  geom_segment(aes(x = 0.4108, y = 14, xend = 0.4108, yend = 28), linetype = 'dashed') + 
  geom_segment(aes(x = 0.3461, y = 28, xend = 0.3461, yend = 42), linetype = 'dashed')
ggsave(p,filename = paste0(figpath,'mtor_score_marginalmeans.pdf'), width = 4, height = 6)
d$sampleid %>% unique
d$cluster = d$integrated_snn_res.0.6

# subset freq barplot 
d$sample_irae
d$sample_irae = ifelse(d$sample_irae == 'NA_healthy', yes = paste(d$sampleid, 'healthy'), no = d$sample_irae)
dsum2 = d  %>%  
  group_by(sample_irae,cluster) %>% 
  tally()

p = 
  ggplot(dsum2, aes(x = sample_irae, y = n, fill = cluster)) + 
  geom_bar(position = 'fill', stat = 'identity',  show.legend = TRUE) +
  theme(axis.text.x = element_text(size = 8)) + 
  theme_bw() + 
  ggsci::scale_fill_d3(palette = 'category20') + 
  xlab("") + ylab("fraction of total cells") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0))
ggsave(p, filename = paste0(figpath, "integrated_cluster_barplot.pdf"), width = 8, height = 5)


# DEcluster markers 
int = FindVariableFeatures(int,assay = 'RNA',selection.method = 'vst')
vf = VariableFeatures(int)
gene_rm2 = vf[grep(pattern = "RP11|MT-|RPS|RPL", vf)]
vf = vf[!vf %in% gene_rm2]

Idents(int) = 'integrated_snn_res.0.6'
de_rna = FindAllMarkers(object = int,
                        assay = 'RNA',
                        features = vf,
                        test.use = 'roc',
                        only.pos = TRUE)
saveRDS(de_rna,file = paste0(datapath,'subjectregressionde_rna_inte.rds'))
de_rna = readRDS(file = here('cd8/healthy_integration/generated_data/subjectregressionde_rna_inte.rds'))

markers = de_rna %>% 
  group_by(cluster) %>% 
  arrange(desc(power)) %>% 
  top_n(1)
  
top = unique(markers$gene)

# simple heatmap 
genes = de_rna %>% 
  filter(!cluster  == '14') %>% 
  group_by(cluster) %>% 
  arrange(desc(power)) %>% 
  top_n(3) %$%  
  gene %>% 
  unique()

# heatmap
dat = int@assays$RNA@data[genes, ] %>% 
  as.matrix() %>% 
  t() %>%
  as.data.frame() %>%
  add_column(int@meta.data)

# plot 
dd_av = dat %>% 
  filter(!integrated_snn_res.0.6 == '14') %>% 
  group_by(integrated_snn_res.0.6) %>% 
  summarize_at(.vars = genes, .funs = mean) %>% 
  column_to_rownames("integrated_snn_res.0.6") %>%
  t()
pheatmap::pheatmap(dd_av,
                   clustering_method = 'single',
                   scale = "row",
                   treeheight_row = 10, treeheight_col = 10, 
                   color = rev(pals::brewer.puor(25)),
                   width = 4.2, height = 5.2,
                   filename = paste0(figpath,"integrated_geneheatmap.pdf"))

# proteins bet. clusters 
de_prot = FindAllMarkers(object = int,
                         assay = 'CITE',
                         test.use = 'roc',
                         only.pos = TRUE)
# arrange discriminitive proteins 
prots = de_prot %>% 
  filter(avg_diff > 1) %>% 
  filter(!cluster %in% c('14', '15')) %>%
  group_by(cluster) %>% 
  arrange(desc(power)) %>%
  arrange(cluster) %>% 
  top_n(3) %$%  
  gene %>% 
  unique()
prots = c(prots, 'CD45RO', 'CD45RA', 'CD279')
datp = int@assays$CITE@data[prots, ] %>% 
  as.matrix() %>% 
  t() %>%
  as.data.frame() %>%
  add_column(int@meta.data)
dd_av = datp %>% 
  filter(!integrated_snn_res.0.6 %in% c('14', '15')) %>% 
  group_by(integrated_snn_res.0.6) %>% 
  summarize_at(.vars = prots, .funs = mean) %>%  
  column_to_rownames("integrated_snn_res.0.6") %>%
  t()
dd_av
pheatmap::pheatmap(dd_av,
                   clustering_method = 'average',
                   #scale = "row",
                   treeheight_row = 10, treeheight_col = 10, 
                   color = viridis::inferno(n = 8),
                   #color = rev(pals::brewer.puor(25)),
                   width = 4.2, height = 2.2,
                   filename = paste0(figpath,"integrated_PROTheatmap.pdf"))


# cluster contingency 
cmat = table(d$celltype, d$integrated_snn_res.0.6) %>%
  log1p() 

cmat = cmat[ ,!colnames(cmat) == '14']
pheatmap::pheatmap(mat = cmat, color = viridis::inferno(n = 10),
                     treeheight_col = 10,treeheight_row = 10, 
                     main = 'log cell frequency', 
                     filename = paste0(figpath, 'integrated_cluster_contingency_thymoma.pdf'),
                     width = 5, height = 4)

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
#   [1] patchwork_1.0.0 scglmmr_0.1.0   here_0.1        Seurat_3.1.5    forcats_0.5.0   stringr_1.4.0   dplyr_1.0.3     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_3.1.0   
# [12] ggplot2_3.3.0   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.1.4               reticulate_1.14          tidyselect_1.1.0         lme4_1.1-21              RSQLite_2.2.0            AnnotationDbi_1.48.0     htmlwidgets_1.5.1       
# [8] grid_3.6.1               BiocParallel_1.20.1      Rtsne_0.15               munsell_0.5.0            codetools_0.2-16         ica_1.0-2                future_1.16.0           
# [15] withr_2.4.0              colorspace_1.4-1         GOSemSim_2.12.1          Biobase_2.46.0           rstudioapi_0.11          stats4_3.6.1             ROCR_1.0-7              
# [22] ggsignif_0.6.0           DOSE_3.12.0              listenv_0.8.0            labeling_0.3             emmeans_1.4.5            urltools_1.7.3           polyclip_1.10-0         
# [29] pheatmap_1.0.12          bit64_0.9-7              farver_2.0.3             rprojroot_1.3-2          TH.data_1.0-10           coda_0.19-4              vctrs_0.3.6             
# [36] generics_0.0.2           R6_2.4.1                 doParallel_1.0.15        graphlayouts_0.7.0       rsvd_1.0.3               isoband_0.2.0            locfit_1.5-9.4          
# [43] pals_1.6                 bitops_1.0-6             fgsea_1.12.0             gridGraphics_0.5-0       assertthat_0.2.1         promises_1.1.0           scales_1.1.0            
# [50] multcomp_1.4-12          ggraph_2.0.2             enrichplot_1.6.1         gtable_0.3.0             npsurv_0.4-0             egg_0.4.5                globals_0.12.5          
# [57] sandwich_2.5-1           tidygraph_1.1.2          rlang_0.4.10             splines_3.6.1            lazyeval_0.2.2           dichromat_2.0-0          broom_0.7.4             
# [64] europepmc_0.3            BiocManager_1.30.10      reshape2_1.4.3           modelr_0.1.6             backports_1.2.1          httpuv_1.5.2             qvalue_2.18.0           
# [71] clusterProfiler_3.14.3   tools_3.6.1              ggplotify_0.0.5          ellipsis_0.3.1           gplots_3.0.3             RColorBrewer_1.1-2       BiocGenerics_0.32.0     
# [78] ggridges_0.5.2           Rcpp_1.0.4               plyr_1.8.6               progress_1.2.2           RCurl_1.98-1.1           prettyunits_1.1.1        ggpubr_0.2.5            
# [85] pbapply_1.4-2            viridis_0.5.1            cowplot_1.0.0            S4Vectors_0.24.3         zoo_1.8-7                haven_2.2.0              ggrepel_0.8.2           
# [92] cluster_2.1.0            colorRamps_2.3           fs_1.3.2                 variancePartition_1.16.1 magrittr_2.0.1           data.table_1.12.8        DO.db_2.9               
# [99] lmtest_0.9-37            triebeard_0.3.0          reprex_0.3.0             RANN_2.6.1               mvtnorm_1.1-0            fitdistrplus_1.0-14      pkgload_1.0.2           
# [106] hms_0.5.3                lsei_1.2-0               mime_0.9                 GSVA_1.34.0              xtable_1.8-4             pbkrtest_0.4-8.6         XML_3.99-0.3            
# [113] readxl_1.3.1             IRanges_2.20.2           gridExtra_2.3            testthat_2.3.2           compiler_3.6.1           maps_3.3.0               KernSmooth_2.23-16      
# [120] crayon_1.4.1             minqa_1.2.4              htmltools_0.4.0          later_1.0.0              geneplotter_1.64.0       lubridate_1.7.4          DBI_1.1.0               
# [127] tweenr_1.0.1             dbplyr_1.4.2             MASS_7.3-51.5            boot_1.3-24              Matrix_1.2-18            cli_2.3.1                gdata_2.18.0            
# [134] parallel_3.6.1           igraph_1.2.5             pkgconfig_2.0.3          rvcheck_0.1.8            plotly_4.9.2             xml2_1.3.2               foreach_1.4.8           
# [141] annotate_1.64.0          estimability_1.3         rvest_0.3.5              digest_0.6.27            sctransform_0.2.1        RcppAnnoy_0.0.16         tsne_0.1-3              
# [148] graph_1.64.0             cellranger_1.1.0         leiden_0.3.3             fastmatch_1.1-0          edgeR_3.28.1             uwot_0.1.8               GSEABase_1.48.0         
# [155] shiny_1.4.0.2            gtools_3.8.1             nloptr_1.2.2.1           lifecycle_1.0.0          nlme_3.1-145             jsonlite_1.6.1           mapproj_1.2.7           
# [162] desc_1.2.0               limma_3.42.2             viridisLite_0.3.0        fansi_0.4.2              pillar_1.5.0             ggsci_2.9                lattice_0.20-40         
# [169] fastmap_1.0.1            httr_1.4.1               survival_3.1-11          GO.db_3.10.0             glue_1.4.2               png_0.1-7                shinythemes_1.1.2       
# [176] iterators_1.0.12         bit_1.1-15.2             ggforce_0.3.1            stringi_1.4.6            blob_1.2.1               org.Hs.eg.db_3.10.0      caTools_1.18.0          
# [183] memoise_1.1.0            irlba_2.3.3              future.apply_1.4.0       ape_5.3   