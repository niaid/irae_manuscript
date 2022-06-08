# baseline gsea contrast result.  
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))


figpath = file.path(here('DE_V4/figures/gsea/')); dir.create(figpath)
datapath = file.path(here('DE_V4/generated_data/temp/'))


# read curated results from g0.result.sort generated in script 4 
d2 = data.table::fread(file = here('DE_V4/generated_data/temp/g0.result.sort.curated.txt'),sep = '\t') %>% 
  filter(include == 1)

# load full gsea results in list format for PlotFgsea 
# filter to only include the temporally stable pathways 
g0.list.combined = readRDS(file = here('DE_V4/generated_data/temp/g0.list.combined.rds'))
test = list()
for (i in 1:length(g0.list.combined)) {
  test[[i]] =
    g0.list.combined[[i]] %>% 
    mutate(signal = paste(celltype, pathway, sep = '~')) %>% 
    filter(signal %in% d2$signal) %>% 
    mutate(class  = ifelse(str_sub(pathway, 1,5) == 'HALLM', 'hallmark', 'hypothesis'))
}
names(test) = names(g0.list.combined)
p = PlotFgsea(gsea_result_list = test, padj_filter = 0.01)

# manually edit plot to include color for rownames and shorten names. 
p$data = p$data %>% mutate(class  = ifelse(str_sub(pathway, 1,5) == 'HALLM', 'hallmark', 'hypothesis'))
p$data$pathway = str_replace_all(p$data$pathway, pattern ="_", replacement = " ")
p$data$pathway = str_replace_all(p$data$pathway, pattern ="HALLMARK", replacement = "")
p$data$celltype = str_replace_all(p$data$celltype, pattern ="_", replacement = " ")

# cluster row col in order to specify color vector of rownames. 
hcdat = p$data %>%
  select(celltype,pathway,NES) %>%
  spread(celltype, NES) %>%
  column_to_rownames("pathway") %>%
  as.matrix()
hcdat[is.na(hcdat)] = 0
xx = pheatmap::pheatmap(hcdat, silent = TRUE, clustering_method = "complete")
module_order = xx$tree_row$labels[xx$tree_row$order]
celltype_order = xx$tree_col$labels[xx$tree_col$order]
p$data$celltype = factor(p$data$celltype, levels = celltype_order)
p$data$pathway = factor(p$data$pathway, levels = rev(module_order))

# define color vector 
te = data.frame(pathway = rev(module_order)) 
te = te %>% mutate(class = plyr::mapvalues(x = pathway, from = p$data$pathway, to = p$data$class))
te$color = ifelse(te$class == 'hallmark', 'slateblue3', 'sienna3')
p = p + theme(axis.text.y = element_text(color = te$color))
ggsave(p,filename = paste0(figpath, 'gsea0.pdf'), width = 7, height = 5.5)

############################
# Correlation of average expression for baseline enrichedd signals 
############################

# correlation among baseline samples for specific leading edge signals 
# average expression of leading edge genes for each subset ; subset to baseline sx
av = readRDS(file = here('DE_V4/generated_data/av.rds'))
d0sx = colnames(av[[1]])[str_sub(colnames(av[[1]]),-1,-1) == '0']
av = lapply(av, function(x){ x[ ,d0sx] } )

li.g0 = LeadingEdgeIndexed(gsea.result.list = test, padj.threshold = 0.01)
li.g0 = base::Filter(length, li.g0)
av = av[names(li.g0)]

res = list()
for (i in 1:length(av)) {
  stopifnot(all.equal( names(av[i]), names(li.g0[i]) ))
  zscore = scglmmr::calc_avg_module_zscore(
    module.list = li.g0[[i]], average.data.frame = av[[i]]
    )
  rownames(zscore) = paste(rownames(zscore), names(av[i]), sep = '~')
  res[[i]] = zscore
}
rm(av); gc()
ds = do.call(rbind, res) %>% t()
colnames(ds) = str_replace_all(colnames(ds), pattern ="_", replacement = " ")
colnames(ds) = str_replace_all(colnames(ds), pattern ="HALLMARK", replacement = "")
saveRDS(ds, file = paste0(datapath, 'ds.rds'))

# correlation matrix of baseline expression across donors 
# of leading edge genes in temporally stable baseline enrichments 
ds = readRDS(here('DE_V4/generated_data/temp/ds.rds'))
pearsonmat = Hmisc::rcorr(ds, type = 'pearson')
spearmanmat = Hmisc::rcorr(ds, type = 'spearman')
rhomat.p = pearsonmat$r
rhomat.s = spearmanmat$r

# total correlation plot 

# make annotation 
ann = data.frame(signal = colnames(rhomat.p)) %>% 
  mutate(sig = signal) %>% 
  separate(signal, into = c('module','celltype'), sep = '~') %>% 
  select(sig, celltype) %>% 
  column_to_rownames('sig')

# use custom palette 
col.use = list(celltype =  unname(pals::polychrome(n = length(unique(ann$celltype)))))
names(col.use[[1]]) = unique(ann$celltype)


HeatmapDiag(rhomat.p, 
            annotation_col = ann, annotation_colors = col.use,
            border_color = NA,
            fontsize_col = 0.00001, fontsize_row = 0.00001, 
            treeheight_col = 0, treeheight_row = 0,
            filename = paste0(figpath, 'PEARSON_corrplot.g0.curated.pdf'),
            width = 6.5, height = 5)


# spearman vs pearson robustness check 
# show spearman rank stat vs the pearson rho 
rp = pheatmap::pheatmap(rhomat.p, 
                        fontsize_row = 3, fontsize_col = 3, 
                        treeheight_col = 0, treeheight_row = 0, 
                        filename = paste0(figpath, 'CHECK_PEARSON_corrplot.g0.curated.pdf'),
                        width = 4, height = 4)
rs = pheatmap::pheatmap(rhomat.s[rp$tree_row$labels[rp$tree_row$order], rp$tree_col$labels[rp$tree_col$order] ],
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        fontsize_row = 3, fontsize_col = 3, 
                        treeheight_col = 0, treeheight_row = 0, 
                        filename = paste0(figpath, 'CHECK_SPEARMAN_PEARSONORDER_corrplot.g0.curated.pdf'),
                        width = 4, height = 4)


# bonferonni adjust rcorr 
source(here('util/MattPMutils.r'))
p.adj = p.adjust.cormat(hmisc.cor = pearsonmat, method = 'fdr')

 
# fdr adjusted p blanked out  
cu = colorRampPalette(c("#4477AA","#77AADD","#FFFFFF", "#EE9988","#BB4444"))
pdf(file = paste0(figpath,"padj.FULLpearson_correlaiton_baseline_avgmodzscore.pdf"), width = 18,height = 14)
corrplot::corrplot(corr = pearsonmat$r,
                   p.mat = p.adj, 
                   insig = 'blank',
                   sig.level = c(0.02), 
                   pch.cex = 0.9, 
                   order = 'hclust',
                   method="color", 
                   col=cu(200),
                   tl.cex = 0.7, 
                   type="upper",
                   tl.col="black",
                   tl.srt=45, 
                   diag=FALSE, 
                   addgrid.col = "grey") 
dev.off()


### subset 
msub = c(
  " MTORC1 SIGNALING~CD8Tcell TEMRA", 
  " CHOLESTEROL HOMEOSTASIS~CD8Tcell TEMRA",
  " TNFA SIGNALING VIA NFKB~CD8Tcell TEMRA",
  " INFLAMMATORY RESPONSE~DC CD1c", 
  " HYPOXIA~DC CD1c", 
  " MTORC1 SIGNALING~DC mDC",
  " HYPOXIA~DC mDC",  
  " INTERFERON ALPHA RESPONSE~DC mDC",
  " INFLAMMATORY RESPONSE~DC mDC")

msub.short.name = c(
  "TEMRA - MTORC1", 
  "TEMRA - CHOLESTEROL",
  "TEMRA - TNFA VIA NFKB",
  "DC CD1c - INFLAMMATORY", 
  "DC CD1c - HYPOXIA", 
  "DC mDC - MTORC1", 
  "DC mDC - HYPOXIA",  
  "DC mDC - INTERFERON ALPHA",
  "DC mDC - INFLAMMATORY"
  )

msub.short.row = c(
    "MTORC1 - TEMRA", 
    "CHOLESTEROL - TEMRA",
    "TNFA VIA NFKB - TEMRA",
    "INFLAMMATORY - DC CD1c", 
    "HYPOXIA- DC CD1c", 
    "MTORC1 - DC mDC", 
    "HYPOXIA - DC mDC",  
    "INTERFERON ALPHA - DC mDC",
    "INFLAMMATORY - DC mDC"
  )

mtxn = pearsonmat$r[msub, msub]
colnames(mtxn) = msub.short.name
rownames(mtxn) = msub.short.row
padjn = p.adj[msub, msub]
colnames(padjn) =  msub.short.name
rownames(padjn) = msub.short.row


pdf(file = paste0(figpath,"padj.0.1.SUBSETpearson_correlaiton_baseline_avgmodzscore.pdf"), width = 5,height = 4)
p1=corrplot::corrplot(corr = mtxn,
                   p.mat = padjn, 
                   insig = 'label_sig',
                   sig.level = c(0.001, 0.01, 0.1), 
                   pch.cex = 0.9, 
                   order = 'hclust',
                   method="circle", 
                   col=cu(10),
                   cl.cex = 0.65,
                   tl.cex = 0.65, 
                   type="upper",
                   tl.col="black",
                   #tl.srt=30, 
                   diag=FALSE, 
                   addgrid.col = "grey") 
dev.off()


# unordered 
pdf(file = paste0(figpath,"padj.0.1.unordered.SUBSETpearson_correlaiton_baseline_avgmodzscore.pdf"), width = 5,height = 4)
p1=corrplot::corrplot(corr = mtxn,
                      p.mat = padjn, 
                      insig = 'label_sig',
                      sig.level = c(0.001, 0.01, 0.1), 
                      pch.cex = 0.9, 
                      order = 'original',
                      #order = 'hclust',
                      method="circle", 
                      col=cu(10),
                      cl.cex = 0.65,
                      tl.cex = 0.65, 
                      type="upper",
                      tl.col="black",
                      #tl.srt=30, 
                      diag=FALSE, 
                      addgrid.col = "grey") 
dev.off()




# make spearman correlation matrix 
p.adj.2 = p.adjust.cormat(hmisc.cor = spearmanmat, method = 'fdr')
pdf(file = paste0(figpath,"padj.0.1.SUBSETspearman_correlaiton_baseline_avgmodzscore.pdf"), width = 4.7,height = 4.7)
corrplot::corrplot(corr = spearmanmat$r[msub, msub],
                   p.mat = p.adj.2[msub, msub], 
                   insig = 'label_sig',
                   sig.level = c(0.001, 0.01, 0.1), 
                   pch.cex = 0.9, 
                   order = 'hclust',
                   method="circle", 
                   col=cu(10),
                   cl.cex = 0.4,
                   tl.cex = 0.4, 
                   type="upper",
                   tl.col="black",
                   tl.srt=30, 
                   diag=FALSE, 
                   addgrid.col = "grey") 
dev.off()


##### correlations 
dx = as.data.frame(ds) %>%
  rownames_to_column('sample') %>% 
  mutate(irae = ifelse(sample %in% c('104_0','201_0','102_0','103_0','101_0'),'irAE', 'noirAE')) %>%  
  select(sample, irae, msub)

pl = list(
  theme_bw(), 
  geom_point(shape = 21, size = 3, aes(fill = irae), show.legend = FALSE),
  theme(axis.title = element_text(size = 6)),
  scale_fill_manual(values = c(col.alpha('red', 0.4), col.alpha('blue', 0.4)))
)

mpl = setdiff(msub, " MTORC1 SIGNALING~CD8Tcell TEMRA")
for (i in 1:length(mpl)) {
  var  = as.name(mpl[i])
  p = ggplot(dx, aes(x = ` MTORC1 SIGNALING~CD8Tcell TEMRA`, y =  {{var}}  )) + pl
  ggsave(p, filename = paste0(figpath, var,'.pdf'), width = 2.5, height = 2.5)
}

mpl = setdiff(msub, " TNFA SIGNALING VIA NFKB~CD8Tcell TEMRA")
for (i in 1:length(mpl)) {
  var  = as.name(mpl[i])
  p = ggplot(dx, aes(x = ` TNFA SIGNALING VIA NFKB~CD8Tcell TEMRA`, y =  {{var}}  )) + pl
  ggsave(p, filename = paste0(figpath,'tnf.temra', var,'.pdf'), width = 2.5, height = 2.5)
}


p = ggplot(dx, aes(x = ` TNFA SIGNALING VIA NFKB~CD8Tcell TEMRA`, y =  ` INFLAMMATORY RESPONSE~DC mDC`  )) + pl
p = ggplot(dx, aes(x = ` TNFA SIGNALING VIA NFKB~CD8Tcell TEMRA`, y =  ` INFLAMMATORY RESPONSE~DC CD1c`  )) + pl
p
ggsave(p, filename = paste0(figpath, var,'.pdf'), width = 2.5, height = 2.5)

## legend 
p = ggplot(dx, aes(x = ` TNFA SIGNALING VIA NFKB~CD8Tcell TEMRA`, y =  ` INFLAMMATORY RESPONSE~DC mDC`  )) + 
  theme_bw() +
  geom_point(shape = 21, size = 3, aes(fill = irae), show.legend = TRUE)+ 
  theme(axis.title = element_text(size = 6))+ 
  scale_fill_manual(values = c(col.alpha('red', 0.4), col.alpha('blue', 0.4))) 
ggsave(p, filename = paste0(figpath,'LEGEND','.PDF'),width = 4 , height=4)



sessionInfo()
# R version 3.6.1 (2019-07-05)
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scglmmr_0.1.0   magrittr_2.0.2  here_0.1        forcats_0.5.0   stringr_1.4.0   dplyr_1.0.3    
# [7] purrr_0.3.3     readr_1.3.1     tidyr_1.0.2     tibble_3.1.0    ggplot2_3.3.0   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.1.4               tidyselect_1.1.0         lme4_1.1-21              htmlwidgets_1.5.1       
# [5] RSQLite_2.2.0            AnnotationDbi_1.48.0     grid_3.6.1               BiocParallel_1.20.1     
# [9] munsell_0.5.0            codetools_0.2-16         withr_2.4.0              colorspace_1.4-1        
# [13] GOSemSim_2.12.1          Biobase_2.46.0           knitr_1.28               rstudioapi_0.11         
# [17] stats4_3.6.1             ggsignif_0.6.0           DOSE_3.12.0              labeling_0.3            
# [21] emmeans_1.4.5            urltools_1.7.3           polyclip_1.10-0          bit64_0.9-7             
# [25] farver_2.0.3             pheatmap_1.0.12          rprojroot_1.3-2          coda_0.19-4             
# [29] vctrs_0.3.6              generics_0.0.2           TH.data_1.0-10           xfun_0.12               
# [33] R6_2.4.1                 doParallel_1.0.15        graphlayouts_0.7.0       locfit_1.5-9.4          
# [37] pals_1.6                 bitops_1.0-6             fgsea_1.12.0             gridGraphics_0.5-0      
# [41] assertthat_0.2.1         promises_1.1.0           scales_1.1.0             nnet_7.3-13             
# [45] multcomp_1.4-12          ggraph_2.0.2             enrichplot_1.6.1         gtable_0.3.0            
# [49] egg_0.4.5                tidygraph_1.1.2          sandwich_2.5-1           rlang_0.4.10            
# [53] slanter_0.2-0            splines_3.6.1            acepack_1.4.1            dichromat_2.0-0         
# [57] checkmate_2.0.0          broom_0.7.4              europepmc_0.3            BiocManager_1.30.10     
# [61] reshape2_1.4.3           modelr_0.1.6             backports_1.2.1          httpuv_1.5.2            
# [65] qvalue_2.18.0            Hmisc_4.4-0              clusterProfiler_3.14.3   tools_3.6.1             
# [69] ggplotify_0.0.5          ellipsis_0.3.1           gplots_3.0.3             RColorBrewer_1.1-2      
# [73] BiocGenerics_0.32.0      ggridges_0.5.2           Rcpp_1.0.4               plyr_1.8.6              
# [77] base64enc_0.1-3          progress_1.2.2           RCurl_1.98-1.1           prettyunits_1.1.1       
# [81] rpart_4.1-15             ggpubr_0.2.5             viridis_0.5.1            cowplot_1.0.0           
# [85] S4Vectors_0.24.3         zoo_1.8-7                cluster_2.1.0            haven_2.2.0             
# [89] ggrepel_0.8.2            colorRamps_2.3           fs_1.3.2                 variancePartition_1.16.1
# [93] data.table_1.12.8        DO.db_2.9                triebeard_0.3.0          reprex_0.3.0            
# [97] mvtnorm_1.1-0            hms_0.5.3                mime_0.9                 GSVA_1.34.0             
# [101] xtable_1.8-4             pbkrtest_0.4-8.6         XML_3.99-0.3             jpeg_0.1-8.1            
# [105] readxl_1.3.1             IRanges_2.20.2           gridExtra_2.3            compiler_3.6.1          
# [109] maps_3.3.0               KernSmooth_2.23-16       crayon_1.4.1             minqa_1.2.4             
# [113] htmltools_0.4.0          mgcv_1.8-31              later_1.0.0              Formula_1.2-3           
# [117] geneplotter_1.64.0       lubridate_1.7.4          DBI_1.1.0                corrplot_0.90           
# [121] tweenr_1.0.1             dbplyr_1.4.2             MASS_7.3-51.5            boot_1.3-24             
# [125] Matrix_1.2-18            cli_2.3.1                gdata_2.18.0             parallel_3.6.1          
# [129] igraph_1.2.5             pkgconfig_2.0.3          rvcheck_0.1.8            foreign_0.8-76          
# [133] xml2_1.3.2               foreach_1.4.8            annotate_1.64.0          GeneOverlap_1.22.0      
# [137] estimability_1.3         rvest_0.3.5              digest_0.6.27            graph_1.64.0            
# [141] cellranger_1.1.0         fastmatch_1.1-0          htmlTable_1.13.3         edgeR_3.28.1            
# [145] GSEABase_1.48.0          shiny_1.4.0.2            gtools_3.8.1             nloptr_1.2.2.1          
# [149] lifecycle_1.0.0          nlme_3.1-145             jsonlite_1.6.1           mapproj_1.2.7           
# [153] viridisLite_0.3.0        limma_3.42.2             fansi_0.4.2              pillar_1.5.0            
# [157] ggsci_2.9                lattice_0.20-40          fastmap_1.0.1            httr_1.4.1              
# [161] survival_3.1-11          GO.db_3.10.0             glue_1.4.2               UpSetR_1.4.0            
# [165] png_0.1-7                shinythemes_1.1.2        iterators_1.0.12         bit_1.1-15.2            
# [169] ggforce_0.3.1            stringi_1.4.6            blob_1.2.1               org.Hs.eg.db_3.10.0     
# [173] latticeExtra_0.6-29      caTools_1.18.0           memoise_1.1.0  
