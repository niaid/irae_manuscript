suppressMessages(library(tidyverse)) 
suppressMessages(library(here)) 
suppressMessages(library(Seurat))
figpath = here('louma/clust.exprsfigs/'); dir.create(figpath)
# data
s = readRDS(file = here('louma/generated_data/cd3_louma_Seurat_processed.rds'))
cd8.clusters = c('1','4','5','7','10','13','14','18')
Idents(s) <-  'seurat_clusters'
s$seurat_clusters = as.character(s$seurat_clusters)
s = subset(s,idents = cd8.clusters)
s = NormalizeData(s)
s = FindVariableFeatures(s, nfeatures = 5000)

# find markers 
de = FindAllMarkers(s, features = c(VariableFeatures(s), 'CD3E', 'CD8A', 'MKI67'), test.use = 'roc')

# group by top n 
de.sub = de %>% 
  group_by(cluster) %>%
  filter(avg_diff > 0) %>%  
  top_n(5, avg_diff) 
de.genes = c(unique(de.sub$gene), 'CD3E', 'CD8A', 'MKI67')
# average log exprs and frac non zero per cluster for de genes 
dat = cbind(as.data.frame(t(as.matrix(s@assays$RNA@data[de.genes, ]))),
                    s@meta.data)


ds= dat %>%
  select(barcode_full, seurat_clusters, de.genes) %>% 
  gather(gene, exprs, de.genes[1]:de.genes[length(de.genes)]) %>% 
  group_by(seurat_clusters, gene) %>% 
  summarize(
    gene_mean = mean(exprs),
    perc_exprs = FSA::perc(exprs, dir = 'gt',val = 0)
  )  


sub = c(

  # c4
  'CD38', 'GZMB', 'LAG3', 'HLA-DRA',
  # C10 
  'STMN1', 'TYMS', 'MKI67',
  # C7 
  'SH2D1A', 
  # c1
  'IL7R', 
  #c5 
  'KIR2DL4', 'FCER1G', 
  # c14 
  'TRAV1-2', 'NRC3', 'KLRB1', 
  # c13 
  'KLRC1', 'AREG', 
  # C18
  'TNFRSF18'
  
  
  
  )
gene.order = c(
  "TNFRSF18",  "GZMA"  ,   
  "TRBV19" ,   "TRAV29DV5", "TRGV2" ,    "AREG",      "TRGV8",    
  "KLRC1",     "TRBV28" ,   "TRAV39",    "NFKBIA",    "KLRB1",    
  "TRAV1-2",   "LTB" ,      "NCR3" ,     "TRDV1" ,    "TRGV4",    
  "FCER1G",    "KLRC2",     "KIR2DL4" ,  "ANXA1" ,    "FOS"  ,    
  "IL7R",      "SH2D1A" ,   "EOMES" ,    "ITGB2" ,    "GZMK",     
  "TUBB",      "HMGB2",     "TYMS" ,     "TUBA1B",    "STMN1" , "MKI67",          
  "HLA-DRA",   "LAG3",      "GZMB",      "LINC02446", "CD38")
  #"CD8A",      "CD3E" )
ds$gene = factor(ds$gene, levels = gene.order)

ds$seurat_clusters  = factor(ds$seurat_clusters, levels = unique(de.sub$cluster))
# p1 = 
#   ggplot(ds, aes(x = gene, y = seurat_clusters , size = perc_exprs, fill = gene_mean)) + 
#   theme_bw() + 
#   coord_flip() + 
#   geom_point(shape =21 ) + 
#   scale_fill_gradientn(colours = viridis::viridis(n = 10,option = 'B'),
#                         limits = c(0.7, 0.9), 
#                         oob = scales::squish,
#                         name = 'avg. expression') 
# p1
ds2 = ds %>% filter(gene %in% sub)
cu = pals::parula(n = 10)
cu = pals::tol.rainbow(n = 10)
cu = pals::brewer.reds(10)
p2 = 
  ggplot(ds2, aes(x = gene, y = seurat_clusters , size = perc_exprs, fill = gene_mean)) + 
  theme_bw() + 
  coord_flip() + 
  geom_point(shape =21) + 
  scale_size_area(max_size = 4) +
  theme(panel.grid.minor = element_line(size = 0.3), panel.grid.major = element_line(size = 0.3)) + 
  scale_fill_gradientn(colours = cu,limits = c(0, 1.6),oob = scales::squish,name = 'avg. expression') + 
  theme(axis.text = element_text(color = 'black', size = 8)) + 
  ylab('CD8 T cell clusters') + 
  xlab("") + 
  theme(axis.ticks = element_blank(),
        legend.key.size = unit(0.3, 'cm'), 
        legend.key.height = unit(0.3, 'cm'), 
        legend.key.width = unit(0.3, 'cm'),
        legend.title = element_text(size=8),  legend.text = element_text(size=8)) 
ggsave(p2,filename = paste0(figpath, 'cluster.exprs.cd8.pdf'), width = 4, height = 2.6)
# make heatmap diag heatmap 

# extract factor rows cols 

# 