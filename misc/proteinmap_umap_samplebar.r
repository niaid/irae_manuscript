# sample visualizations 
# denoised / normalized proetin expression within protein based clusters 
# sample level frequency within clusters -- main lineages and CD8 subsets 
# umap plots for main lineage and CD8 cells 
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(BuenColors))



# seve path 
figpath = here("misc/umap/figures/"); dir.create(figpath)
datapath = here("misc/umap/generated_data/"); dir.create(datapath)

# process data 
# s = readRDS("data/processed_data/ThymomaSeurat.20200422.clustered.mdprocessed.rds")
s = readRDS("data/processed_data/ThymomaSeurat.20200422.clustered.mdprocessed.rds")
d = cbind(s@meta.data, s@reductions$umap@cell.embeddings, as.data.frame(t(s@assays$CITE@data)))
d = d %>% 
  filter(!celltype %in% c('dblt', 'CD8Tcell_IgDpos', 'mixed') ) %>% 
  filter(!timepoint == '2') 
rm(s); gc()

# update cluster annotation to include coarser  clusters and main lineage. 
cdat  = data.table::fread("git_ignore/cluster_annotation_05122020.txt")
d = d %>% mutate(lineage = plyr::mapvalues(x = celltype, from = cdat$Cluster_name, to = cdat$Category)) 
d$adverse_event = ifelse(d$irae == '1', "IRAE", "NO IRAE")
d$sample_irae = paste(d$adverse_event, d$sample, sep = "_")
d$lineage = factor(d$lineage)
saveRDS(d, file = paste0(datapath, 'd.rds'))


# sample level bar frequency - main lineage 
# sum counts per subject for proportion plot
d = readRDS(file = here("misc/umap/generated_data/d.rds"))
dsum2 = d  %>%
  group_by(sample_irae, lineage) %>% 
  tally()

# make sample labels clear 
dsum2$sample_irae = factor(dsum2$sample_irae)
sample_names = 
  c("P1 Baseline", "P1 Post Treatment", "P2 Baseline", "P2 Post Treatment",    
  "P3 baseline", "P3 Post Treatment", "P4 Baseline", "P4 Post Treatment", 
  "P5 Baseline", "P5 Post Treatment", "P6 Baseline", "P6 Post Treatment", 
  "P7 Baseline", "P7 Post Treatment", "P8 Baseline", "P8 Post Treatment", 
  "P9 Baseline", "P9 Post Treatment")
levels(dsum2$sample_irae) = sample_names


# Sample bar freq for main lineage - each sample 
p = ggplot(dsum2, aes(x = sample_irae, y = n, fill = lineage)) + 
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  geom_bar(position = 'fill', stat = 'identity',  show.legend = TRUE) +
  theme(axis.text.x = element_text(size = 0)) + 
  ggsci::scale_fill_igv(alpha = 0.8) + 
  theme(legend.key.size = unit(0.4, units = 'cm'), legend.position = 'top',
        legend.text = element_text(size = 11.5), legend.title = element_blank()) + 
  #theme(axis.text.x=element_text(angle = -90, hjust = 0, color = 'black', size = 12)) + 
  xlab("") + ylab("") +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave(p,filename = paste0(figpath, "sample_barplot_lineage.pdf"), width = 5.5, height = 6)



theme_bw() +
  theme(panel.grid = element_blank()) + 
  geom_bar(position = 'fill', stat = 'identity',  show.legend = TRUE) +
  theme(axis.text.x = element_text(size = 0)) + 
  scale_fill_manual(values = cu) + 
  xlab("") + ylab("") +
  theme(legend.key.size = unit(0.4, units = 'cm' ), legend.position = 'bottom',
        legend.text = element_text(size = 11.5), legend.title = element_blank()) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, color = 'black', size = 12)) + 
  guides(fill=guide_legend(ncol=2,byrow=TRUE))




# umap - all cells color main lineage. 
# stripped theme 
boxbox = list(theme_bw(), 
              theme(
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(), 
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.x = element_blank(), 
                axis.text.y = element_blank())
) 
# calc centers
centers = d %>% 
  dplyr::group_by(lineage) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
d$celltype = str_replace(d$celltype, pattern = "_",replacement = " ")

p = 
  ggplot(d, aes( x = UMAP_1 , y = UMAP_2, fill = lineage)) +
  boxbox + 
  ggsci::scale_fill_igv(alpha = 0.8) + 
  geom_point(size = 0.3, shape = 21, stroke = 0,  show.legend = FALSE) 
ggsave(p,filename = paste0(figpath, 'lineage_umap2.png'), width = 5, height = 4)



# protein heatmap of dsb normalized level exprs all clusters 
# While d has all cells (processed further below)
# vector of high information proteins to keep 
prot_sub = c(
              # T 
              'CD3', 'CD4', 'CD8', 'CD27',  
              # T diff 
              'CD62L', 'CD38', 'CD28', 'CD279', 
              'CD45RA', 'CD45RO', 
              'CD161',  'CD103', 
              # NK 
              'CD56', 
              # B 
              'CD19', 'IgD', 'IgM',
              # myel
              'CD34', 'CD71', 
              'CD14', 'CD16', 
              'CD1c','CD11b', 'CD11c', 'CD1d', 'CD64', 'HLA-DR', 
              'CD123', 
              # state
              'CD39', 'CD31', 'CD86', 'CD69'
              )
mtx = d %>% 
  group_by(celltype) %>% 
  summarize_at(.vars = prot_sub, .funs = mean) %>% 
  column_to_rownames('celltype')
xx = pheatmap::pheatmap(mtx, clustering_method = 'average')

# extract celltype order 
# prot_order = xx$tree_col$labels[xx$tree_col$order]
celltype_order = xx$tree_row$labels[xx$tree_row$order]

# summarize by percent expression 
index1 = prot_sub[1]; index2 = prot_sub[length(prot_sub)]

# calculate percent of cells above the 3.5 sd dsb threshold 
pgt = function(x){FSA::perc(x = x,dir = "gt",val = 3.5)}
dsummary = 
  d %>% 
  group_by(celltype) %>% 
  summarize_at(.vars = prot_sub, .funs = pgt) %>% 
  gather(prot, pct_express, index1:index2)

# calculate average dsb expression per cluster 
d_av = 
  d %>% 
  group_by(celltype) %>% 
  summarize_at(.vars = prot_sub, .funs = mean) %>% 
  gather(prot, gene_avg, index1:index2)

# merge summarized data 
d2 = full_join(dsummary, d_av)
d2$celltype = factor(d2$celltype, levels = rev(celltype_order))
d2$prot = factor(d2$prot, levels = rev(prot_sub))

# draw plot 
p = ggplot(d2 %>% filter(pct_express > 3.5), aes(x = celltype, y = prot, fill = gene_avg, size = pct_express)) +
  geom_point(shape = 21) + 
  scale_fill_gradientn(colours = BuenColors::jdb_palettes$solar_rojos,
                       limits = c(0,11), 
                       oob = scales::squish,
                       name = 'mean dsb\nnormalized \nvalue') + 
  scale_size_area(name = '% of cells \nabove noise \nthreshold') + 
  theme_bw() + 
  theme(axis.line  = element_blank()) +
  scale_x_discrete(position = 'top') + 
  theme(axis.text.x = element_text(angle = 58, vjust = 1, hjust=0, color = 'black', size = 11)) +
  theme(axis.text.y = element_text(color = 'black', size = 11)) +
  ylab('') + xlab('') + 
  theme(legend.position = 'bottom') + 
  theme(axis.ticks = element_blank())  +
  theme(plot.margin = unit(c(0.1, 1.5, 0.2,0.2), 'cm'))
ggsave(p,filename = paste0(figpath, 'protein_map_thymoma.pdf'), width = 9.1, height = 8.2)



# umap - T cell subsets 
# Process data and save sample bar plot and umap 
## Get vector of T cells to run umap and re load processed data 
tc = rownames(d[d$lineage %in% c("Tcell", "CD8Tcell"), ])
s = readRDS("data/processed_data/ThymomaSeurat.20200422.clustered.mdprocessed.rds")
s = subset(s, cells = tc); gc()

# define TC proteins to run umap on 
tcp = c("CD103", "CD127", "CD161", "CD25", "CD27", "CD273", "CD274", "CD278", "CD279",
        "CD28", "CD294", "CD3", "CD31", "CD38", "CD39", "CD4", "CD40", "CD45RA", "CD45RO",
        "CD5", "CD56", "CD57", "CD62L", "CD69", "CD7", "CD8", "CD86", "KLRG1")

# run umap
set.seed(1990)
Idents(s) <-  "celltype"
s = FindNeighbors(s, assay = "CITE",features = tcp, dims = NULL)
s = RunUMAP(s, assay = "CITE", graph = 'CITE_snn', n.neighbors = 40, min.dist = 0.2, n.epochs = 100)

# save embeddings and tcmd 
tcd = cbind(s@meta.data, s@reductions$umap@cell.embeddings, as.data.frame(t(s@assays$CITE@data)))
tct = tcd %>% mutate(lineage = plyr::mapvalues(x = celltype, from = cdat$Cluster_name, to = cdat$Category)) 
tcd$adverse_event = ifelse(tcd$irae == '1', "IRAE", "NO IRAE")
tcd$sample_irae = paste(tcd$adverse_event, tcd$sample, sep = "_")
tcd$celltype = factor(tcd$celltype)
saveRDS(tcd, file = paste0(datapath, 'tcd.rds'))

# read data saved above
tcd = readRDS(file = here('misc/umap/generated_data/tcd.rds'))
tcd$celltype = str_replace_all(tcd$celltype, pattern = "_", replacement = " ")


# Bar plot CD8  - sample level 
dsum3 = tcd  %>%
  group_by(sample_irae, celltype) %>% 
  tally()
cu = BuenColors::jdb_palettes$lawhoops
dsum3$sample_irae = factor(dsum3$sample_irae)
levels(dsum3$sample_irae) = sample_names
dsum3$celltype[dsum3$celltype == 'CD4CD8dpT CD161CD56pos'] = 'CD4CD8dpT CD161CD56+'
dsum3$celltype[dsum3$celltype == 'CD8Tcell naive.CD28neg'] = 'CD8Tcell naive.CD28-'
dsum3$celltype[dsum3$celltype == 'CD8Tcell naive.CD27neg'] = 'CD8Tcell naive.CD27-'
p = ggplot(dsum3, aes(x = sample_irae, y = n, fill = celltype)) + 
  theme_bw() +
  theme(panel.grid = element_blank()) + 
  geom_bar(position = 'fill', stat = 'identity',  show.legend = TRUE) +
  theme(axis.text.x = element_text(size = 8)) + 
  scale_fill_manual(values = cu) + 
  xlab("") + ylab("") +
  theme(legend.key.size = unit(0.4, units = 'cm' ), legend.position = 'bottom',
        legend.text = element_text(size = 11.5), legend.title = element_blank()) + 
  theme(axis.text.x=element_text(angle = -90, hjust = 0, color = 'black', size = 12)) + 
  guides(fill=guide_legend(ncol=2,byrow=TRUE))
ggsave(p,filename = paste0(figpath, "CD8sample_barplot_lineage.pdf"), width = 5, height = 9)


# umap CD8 
source('util/plot.utility.r')
cu2 = sapply(cu, col.alpha, 0.2) %>% unname()
centers = tcd %>% 
  dplyr::group_by(celltype) %>% 
  summarize(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
# draw plot 
p = 
  ggplot(tcd %>% filter(UMAP_1 > -10), aes( x = UMAP_1 , y = UMAP_2, fill = celltype)) +
  boxbox + 
  scale_fill_manual(values = cu2) + 
  geom_point(size = 1.2, shape = 21, stroke = 0,  show.legend = FALSE) +
  ggrepel::geom_label_repel(data = centers %>% filter(!celltype %in% 'Tcell CD71pos'),  
                            aes( x = UMAP_1 , y = UMAP_2, label = celltype),
                            force = 140, 
                            size = 4.3, segment.color = 'black', segment.alpha = 0.8,
                            label.padding = 0.15,
                            box.padding = 0,
                            max.iter = 80000,seed = 190, 
                            inherit.aes = FALSE, show.legend = FALSE) 
ggsave(p,filename = paste0(figpath, 'CD8lineage_umap3.png'), width = 6, height = 6)


sessionInfo()
# R version 3.6.1 (2019-07-05)
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] BuenColors_0.5.6 MASS_7.3-51.5    magrittr_2.0.2   Seurat_3.1.5     scglmmr_0.1.0    here_0.1         forcats_0.5.0   
# [8] stringr_1.4.0    dplyr_1.0.3      purrr_0.3.3      readr_1.3.1      tidyr_1.0.2      tibble_3.1.0     ggplot2_3.3.0   
# [15] tidyverse_1.3.0 
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
# [57] rlang_0.4.10             slanter_0.2-0            FSA_0.9.0                splines_3.6.1           
# [61] lazyeval_0.2.2           broom_0.7.4              europepmc_0.3            BiocManager_1.30.10     
# [65] reshape2_1.4.3           modelr_0.1.6             backports_1.2.1          httpuv_1.5.2            
# [69] qvalue_2.18.0            clusterProfiler_3.14.3   tools_3.6.1              ggplotify_0.0.5         
# [73] ellipsis_0.3.1           gplots_3.0.3             RColorBrewer_1.1-2       BiocGenerics_0.32.0     
# [77] ggridges_0.5.2           Rcpp_1.0.4               plyr_1.8.6               progress_1.2.2          
# [81] RCurl_1.98-1.1           prettyunits_1.1.1        ggpubr_0.2.5             pbapply_1.4-2           
# [85] viridis_0.5.1            cowplot_1.0.0            S4Vectors_0.24.3         zoo_1.8-7               
# [89] cluster_2.1.0            haven_2.2.0              ggrepel_0.8.2            colorRamps_2.3          
# [93] fs_1.3.2                 variancePartition_1.16.1 data.table_1.12.8        DO.db_2.9               
# [97] lmtest_0.9-37            triebeard_0.3.0          reprex_0.3.0             RANN_2.6.1              
# [101] mvtnorm_1.1-0            matrixStats_0.56.0       fitdistrplus_1.0-14      patchwork_1.0.0         
# [105] lsei_1.2-0               hms_0.5.3                mime_0.9                 GSVA_1.34.0             
# [109] xtable_1.8-4             pbkrtest_0.4-8.6         XML_3.99-0.3             readxl_1.3.1            
# [113] IRanges_2.20.2           gridExtra_2.3            compiler_3.6.1           KernSmooth_2.23-16      
# [117] crayon_1.4.1             minqa_1.2.4              htmltools_0.4.0          later_1.0.0             
# [121] snow_0.4-3               geneplotter_1.64.0       lubridate_1.7.4          DBI_1.1.0               
# [125] tweenr_1.0.1             corrplot_0.90            dbplyr_1.4.2             boot_1.3-24             
# [129] Matrix_1.2-18            cli_2.3.1                gdata_2.18.0             parallel_3.6.1          
# [133] igraph_1.2.5             pkgconfig_2.0.3          rvcheck_0.1.8            plotly_4.9.2            
# [137] xml2_1.3.2               foreach_1.4.8            annotate_1.64.0          GeneOverlap_1.22.0      
# [141] estimability_1.3         rvest_0.3.5              digest_0.6.27            tsne_0.1-3              
# [145] sctransform_0.2.1        RcppAnnoy_0.0.16         graph_1.64.0             leiden_0.3.3            
# [149] cellranger_1.1.0         fastmatch_1.1-0          uwot_0.1.8               edgeR_3.28.1            
# [153] GSEABase_1.48.0          shiny_1.4.0.2            gtools_3.8.1             nloptr_1.2.2.1          
# [157] lifecycle_1.0.0          nlme_3.1-145             jsonlite_1.6.1           viridisLite_0.3.0       
# [161] limma_3.42.2             fansi_0.4.2              pillar_1.5.0             ggsci_2.9               
# [165] lattice_0.20-40          fastmap_1.0.1            httr_1.4.1               survival_3.1-11         
# [169] GO.db_3.10.0             glue_1.4.2               png_0.1-7                shinythemes_1.1.2       
# [173] iterators_1.0.12         bit_1.1-15.2             ggforce_0.3.1            stringi_1.4.6           
# [177] blob_1.2.1               org.Hs.eg.db_3.10.0      caTools_1.18.0           memoise_1.1.0           
# [181] irlba_2.3.3              future.apply_1.4.0       ape_5.3   





