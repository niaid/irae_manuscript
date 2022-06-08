# visulaziaiton of gsea contrast results 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))

figpath = file.path(here('DE_V4/figures/contrasts/')); dir.create(figpath)
datapath = file.path(here('DE_V4/generated_data/contrasts/')); dir.create(datapath)
source(here('util/MattPMutils.r'))


# treatment contrast
# highlight results from treatment contrast with focus on results from T cell substs. 
hyp.treat = readRDS(here("DE_V4/generated_data/hyp.treat.rds"))
hlmk.treat = readRDS(here("DE_V4/generated_data/hlmk.treat.rds"))
g1.list.combined = Map(f = rbind, hyp.treat, hlmk.treat)


# base theme 
mtheme = list(
  theme_bw(), 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
)

d = RbindGseaResultList(g1.list.combined,NES_filter = -Inf,padj_filter = 0.01)
d2 = d %>% 
  filter(celltype %in% c('CD8Tcell_TEMRA', 'CD8Tcell_TEM', 'CD8Tcell_naive', 'DNT_CD45RAlo',
                         'CD8Tcell_Activated',  'DC_mDC', 'CD4Tcell_CM','CD4Tcell_TM', 'CD4Tcell_EM',
                         'CD4Tcell_naive', 'DNT_CD27neg' )) %>% 
  filter(! pathway %in% "HALLMARK_TNFA_SIGNALING_VIA_NFKB")

d2$celltype = str_replace_all(d2$celltype, pattern ="_", replacement = " ")
d2$pathway = str_replace_all(d2$pathway, pattern ="_", replacement = " ")
d2$pathway = str_replace_all(d2$pathway, pattern ="HALLMARK", replacement = "")
d2$pathway[d2$pathway == 'SLE SIG'] = 'SLE Interferon Sig'
p = 
ggplot(d2, aes(x = NES, y = reorder(pathway, NES),
               label = celltype, group = celltype, fill = celltype)) + 
  mtheme + 
  geom_linerange(aes(x = NES, color = celltype, xmin = 0, xmax = NES),
                 position = position_dodge(width = 1)) + 
  geom_point(aes(x = NES, color = celltype),position = position_dodge(width = 1)) + 
  geom_point(shape = 21, size = 2.5, position = position_dodge(width = 1)) + 
  ggsci::scale_fill_d3(alpha = 0.5) + 
  ggsci::scale_color_d3(alpha = 0.5) + 
  theme(axis.text = element_text(color = 'black')) + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  ylab('') + 
  xlab('Normalized Enrichment Score') + 
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 12))
ggsave(p,filename = paste0(figpath, 'treatment.pdf'), width = 7, height = 4)  

# visualize disribution of activated CD8 t cell phenotype genes 
# upset plot of mutual information between leading edge genes in enrichment. 
li.treat = LeadingEdgeIndexed(g1.list.combined, padj.threshold = 0.01)
names(li.treat$CD8Tcell_Activated) = str_replace_all(names(li.treat$CD8Tcell_Activated),
                                                      pattern ="_", replacement = " ")
names(li.treat$CD8Tcell_Activated) = str_replace_all(names(li.treat$CD8Tcell_Activated),
                                                      pattern ="HALLMARK", replacement = "")
li.treat$CD8Tcell_Activated = li.treat$CD8Tcell_Activated[-6]

# draw plot 
pdf(file = paste0(figpath, 'CD8Activated_tratment_.upsetplot.pdf'), width = 6, height = 4.5)
UpSetR::upset(UpSetR::fromList(li.treat$CD8Tcell_Activated), nsets = 8,
              sets.x.label = 'n Leading Edge Genes',  
              point.size=3.5,
              matrix.color="#e97572", 
              main.bar.color = '#e97572')
dev.off()


# average expression gene distribution for select genes in leading edge. 
# load averge data and sample metadata 
av = readRDS(file = here('DE_V4/generated_data/av.rds'))
samplemd = readRDS(file = here('DE_V4/generated_data/samplemd.rds')) %>% 
  rownames_to_column('sampleID')

dav2 = scglmmr::LeadEdgeTidySampleExprs(av.exprs.list = av,
                                        g1.list.combined, 
                                        padj.filter = 0.01, 
                                        NES.filter = 0)
fit1 = readRDS(file = here('DE_V4/generated_data/fit1.rds'))
fit1.res = ExtractResult(fit1, coefficient.number = 3, coef.name = 'treatment')
fit1.res$CD8Tcell_Activated %>% filter(gene %in% unique(unlist(li.treat$CD8Tcell_Activated)))
li.treat$CD8Tcell_Activated
dav2 = dav2$CD8Tcell_Activated
dav2 = dav2 %>% 
  mutate(timepoint = plyr::mapvalues( 
    x = dav2$sample,from = samplemd$sampleID, samplemd$timepoint))
dt = dav2 %>% filter(gene %in% c('PDCD1', 'LCK', 'MKI67', "CDC27", "CDC25B",
                                 "CKS2", 'OAS1', 'IFITM2', 'MCM3', 'GZMA'))
dt$timepoint[dt$timepoint == 0] = 'baseline'
dt$timepoint[dt$timepoint == 1] = 'post-treatment'
p = 
  ggplot(data = dt, aes(x = timepoint, fill = timepoint, color = timepoint, y = av_exp)) + 
  mtheme + 
  geom_boxplot(outlier.shape = NA) +  
  facet_wrap(~gene, scales = 'free_y', nrow = 1) + 
  ggsci::scale_fill_jama(alpha = 0.3) + 
  ggsci::scale_color_jama() + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 0)
        )  + 
  ylab('Average Expression')
ggsave(p,filename =  paste0(figpath, 'time.CD8activ.genes.pdf'), width = 8, height = 2)

################################################

# delta contrast 
## highlight select genes from Classical monocyte IFN enrichments 
hyp.delta = readRDS(here("DE_V4/generated_data/hyp.delta.rds"))
hlmk.delta = readRDS(here("DE_V4/generated_data/hlmk.delta.rds"))
gd.list.combined = Map(f = rbind, hyp.delta, hlmk.delta)
d = RbindGseaResultList(gd.list.combined,NES_filter = -Inf,padj_filter = 0.1)
d2 = d %>% 
  filter(celltype %in% c('Mono_classical', "CD8Tcell_TEM", "CD4Tcell_naive" ))

d2$celltype = str_replace_all(d2$celltype, pattern ="_", replacement = " ")
d2$pathway = str_replace_all(d2$pathway, pattern ="_", replacement = " ")
d2$pathway = str_replace_all(d2$pathway, pattern ="HALLMARK", replacement = "")
d2 = d2 %>% filter(! pathway %in% c(' DNA REPAIR', ' UV RESPONSE UP'))
d2$pathway[d2$pathway == "Reactome antigen processing and proteasome"] = "Reactome antigen processing"

d2$celltype = factor(d2$celltype, levels = c("Mono classical" , "CD4Tcell naive" , "CD8Tcell TEM" )) 
p = 
  ggplot(d2, aes(x = NES, y = reorder(pathway, NES), 
                 label = celltype, group = celltype, fill = celltype)) + 
  mtheme + 
  geom_linerange(aes(x = NES, color = celltype, xmin = 0, xmax = NES),
                 position = position_dodge(width = 1)) + 
  geom_point(aes(x = NES, color = celltype),position = position_dodge(width = 1)) + 
  geom_point(shape = 21, size = 2, position = position_dodge(width = 1)) + 
  ggsci::scale_fill_d3(alpha = 0.5) + 
  ggsci::scale_color_d3(alpha = 0.5) + 
  theme(axis.text = element_text(color = 'black')) + 
  geom_vline(xintercept = 0, linetype = 'dashed') + 
  ylab('') + 
  xlab('Normalized Enrichment Score') + 
  theme(axis.title.x = element_text(size = 14)) + 
  theme(axis.text.y = element_text(size = 12))
ggsave(p,filename = paste0(figpath, 'delta.pdf'), width = 7, height = 4)  

# upset plot of leading edge genes 
li.delta = LeadingEdgeIndexed(gd.list.combined,padj.threshold = 0.01)
li.delta$Mono_classical
names(li.delta$Mono_classical) = str_replace_all(names(li.delta$Mono_classical),
                                                 pattern ="_", replacement = " ")
names(li.delta$Mono_classical) = str_replace_all(names(li.delta$Mono_classical),
                                                 pattern ="HALLMARK", replacement = "")
pdf(file = paste0(figpath, 'mono.upsetplot.pdf'), width = 6, height = 3.5)
UpSetR::upset(UpSetR::fromList(li.delta$Mono_classical),sets.x.label = 'n Leading Edge Genes',  
                  point.size=5,
                  matrix.color="#72a8d1", 
                  main.bar.color = '#72a8d1')
dev.off()


# average monocyte exprs 
dav = scglmmr::LeadEdgeTidySampleExprs(av.exprs.list = av,
                                       gd.list.combined, 
                                       padj.filter = 0.01, 
                                       NES.filter = 0)
dav.mono = dav$Mono_classical
dav.mono = dav.mono %>% 
  mutate(group = plyr::mapvalues(x = dav.mono$sample,
                                 from = samplemd$sampleID, samplemd$irae_timepoint))
dt = dav.mono%>% filter(gene %in% c('B2M', 'CASP4', 'CD38', 'IL15', 'IL7' , 'NFKB1' ,'PSM8', 'STAT1', 'VAMP5'))
cu = c('#51a0fd', '#372282', '#ff290d', '#c1372f')
cu.alpha = sapply(cu, col.alpha, alpha = 0.2) %>% unname()
p =
  ggplot(data = dt, aes(x = group, fill = group,color = group, y = av_exp)) + 
  mtheme + 
  geom_boxplot(outlier.shape = NA) +  
  facet_wrap(~gene, scales = 'free_y', nrow = 1) + 
  scale_color_manual(values = cu) +
  scale_fill_manual(values = cu.alpha) + 
  theme(strip.background = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 0))  + 
  ylab('Average Expression')

ggsave(p,filename =  paste0(figpath, 'delta.mono.genes.pdf'), width = 6.5, height = 2)


sessionInfo()
# R version 3.6.1 (2019-07-05)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] scglmmr_0.1.0   magrittr_2.0.2  here_0.1        forcats_0.5.0   stringr_1.4.0   dplyr_1.0.3     purrr_0.3.3     readr_1.3.1     tidyr_1.0.2    
# [10] tibble_3.1.0    ggplot2_3.3.0   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.1.4               tidyselect_1.1.0         lme4_1.1-21              RSQLite_2.2.0            AnnotationDbi_1.48.0    
# [6] htmlwidgets_1.5.1        grid_3.6.1               BiocParallel_1.20.1      munsell_0.5.0            codetools_0.2-16        
# [11] withr_2.4.0              colorspace_1.4-1         GOSemSim_2.12.1          Biobase_2.46.0           knitr_1.28              
# [16] rstudioapi_0.11          stats4_3.6.1             ggsignif_0.6.0           DOSE_3.12.0              labeling_0.3            
# [21] emmeans_1.4.5            urltools_1.7.3           polyclip_1.10-0          bit64_0.9-7              farver_2.0.3            
# [26] pheatmap_1.0.12          rprojroot_1.3-2          coda_0.19-4              vctrs_0.3.6              generics_0.0.2          
# [31] TH.data_1.0-10           xfun_0.12                doParallel_1.0.15        R6_2.4.1                 graphlayouts_0.7.0      
# [36] locfit_1.5-9.4           bitops_1.0-6             fgsea_1.12.0             gridGraphics_0.5-0       assertthat_0.2.1        
# [41] promises_1.1.0           scales_1.1.0             multcomp_1.4-12          ggraph_2.0.2             nnet_7.3-13             
# [46] enrichplot_1.6.1         gtable_0.3.0             egg_0.4.5                tidygraph_1.1.2          sandwich_2.5-1          
# [51] rlang_0.4.10             slanter_0.2-0            splines_3.6.1            acepack_1.4.1            broom_0.7.4             
# [56] europepmc_0.3            checkmate_2.0.0          BiocManager_1.30.10      reshape2_1.4.3           modelr_0.1.6            
# [61] backports_1.2.1          httpuv_1.5.2             qvalue_2.18.0            Hmisc_4.4-0              clusterProfiler_3.14.3  
# [66] tools_3.6.1              ggplotify_0.0.5          ellipsis_0.3.1           gplots_3.0.3             RColorBrewer_1.1-2      
# [71] BiocGenerics_0.32.0      ggridges_0.5.2           Rcpp_1.0.4               plyr_1.8.6               base64enc_0.1-3         
# [76] progress_1.2.2           RCurl_1.98-1.1           prettyunits_1.1.1        ggpubr_0.2.5             rpart_4.1-15            
# [81] viridis_0.5.1            cowplot_1.0.0            S4Vectors_0.24.3         zoo_1.8-7                haven_2.2.0             
# [86] ggrepel_0.8.2            cluster_2.1.0            colorRamps_2.3           fs_1.3.2                 variancePartition_1.16.1
# [91] data.table_1.12.8        DO.db_2.9                triebeard_0.3.0          reprex_0.3.0             mvtnorm_1.1-0           
# [96] hms_0.5.3                mime_0.9                 GSVA_1.34.0              xtable_1.8-4             pbkrtest_0.4-8.6        
# [101] XML_3.99-0.3             jpeg_0.1-8.1             readxl_1.3.1             IRanges_2.20.2           gridExtra_2.3           
# [106] compiler_3.6.1           KernSmooth_2.23-16       crayon_1.4.1             minqa_1.2.4              htmltools_0.4.0         
# [111] later_1.0.0              Formula_1.2-3            geneplotter_1.64.0       lubridate_1.7.4          DBI_1.1.0               
# [116] tweenr_1.0.1             corrplot_0.90            dbplyr_1.4.2             MASS_7.3-51.5            boot_1.3-24             
# [121] Matrix_1.2-18            cli_2.3.1                gdata_2.18.0             parallel_3.6.1           igraph_1.2.5            
# [126] pkgconfig_2.0.3          rvcheck_0.1.8            foreign_0.8-76           foreach_1.4.8            xml2_1.3.2              
# [131] annotate_1.64.0          GeneOverlap_1.22.0       estimability_1.3         rvest_0.3.5              digest_0.6.27           
# [136] graph_1.64.0             cellranger_1.1.0         fastmatch_1.1-0          htmlTable_1.13.3         edgeR_3.28.1            
# [141] GSEABase_1.48.0          shiny_1.4.0.2            gtools_3.8.1             nloptr_1.2.2.1           lifecycle_1.0.0         
# [146] nlme_3.1-145             jsonlite_1.6.1           viridisLite_0.3.0        limma_3.42.2             fansi_0.4.2             
# [151] pillar_1.5.0             ggsci_2.9                lattice_0.20-40          fastmap_1.0.1            httr_1.4.1              
# [156] survival_3.1-11          GO.db_3.10.0             glue_1.4.2               UpSetR_1.4.0             iterators_1.0.12        
# [161] png_0.1-7                shinythemes_1.1.2        bit_1.1-15.2             ggforce_0.3.1            stringi_1.4.6           
# [166] blob_1.2.1               org.Hs.eg.db_3.10.0      latticeExtra_0.6-29      caTools_1.18.0           memoise_1.1.0 
