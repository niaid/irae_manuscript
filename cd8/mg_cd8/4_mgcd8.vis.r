### read bulk results and data and do downsampling 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))

figpath = here("cd8/mg_cd8/figuresV4/")
datapath = here("cd8/mg_cd8/generated_dataV4/")


#sx.li = LeadingEdgeIndexed(gsea.result.list = irae_gsea_samples, padj.threshold = Inf)
irae_gsea_samples = readRDS(here('cd8/mg_cd8/generated_dataV4/irae_gsea_samples.rds'))
ds = do.call(rbind, irae_gsea_samples)
ds$pathway = str_remove_all(ds$pathway, pattern = 'HALLMARK')
ds$pathway = str_replace_all(ds$pathway, pattern = '_', replacement = ' ')
dens =  list(theme_bw(), 
             ggsci::scale_fill_npg(alpha = 0.3), ggsci::scale_color_npg(),
             geom_vline(xintercept = 0, linetype = 'dashed'), 
             xlab('Distribution of Normalized Enrichent Scores \n n=100 k=45 cell resampled model fits \n gsea on genes ranked by IRAE vs no IRAE coefficient')
             )
p = ggplot(ds, aes(x = NES, fill = pathway, color = pathway)) +
  geom_density(show.legend = FALSE) + dens + 
  theme(strip.text = element_text(size = 6.5), strip.background = element_blank()) +
  facet_wrap(~pathway, nrow = 1) 
ggsave(p,filename = paste0(figpath,'filterbyexp.temragated.subsamples.pdf'), width = 12, height = 2.5)






# irae vs healthy 
healthy_gsea_samples = readRDS(here('cd8/mg_cd8/generated_dataV4/healthy_gsea_samples.rds'))
ds = do.call(rbind, healthy_gsea_samples)
ds$pathway = str_remove_all(ds$pathway, pattern = 'HALLMARK')
ds$pathway = str_replace_all(ds$pathway, pattern = '_', replacement = ' ')

p = ggplot(ds, aes(x = NES, fill = pathway, color = pathway)) +
  geom_density(show.legend = FALSE) + dens + 
  facet_wrap(~pathway) + 
  theme(strip.text = element_text(size = 6.5), strip.background = element_blank())
ggsave(p,filename = paste0(figpath,'iraevshealthy.temragated.subsamples.pdf'), width = 5, height = 4.5)


# non filterbyexprs version 
irae_gsea_filter = readRDS(here('cd8/mg_cd8/generated_dataV4/irae_gsea_filter.rds'))

ds = do.call(rbind, irae_gsea_samples)
ds$pathway = str_remove_all(ds$pathway, pattern = 'HALLMARK')
ds$pathway = str_replace_all(ds$pathway, pattern = '_', replacement = ' ')

p = ggplot(ds, aes(x = NES, fill = pathway, color = pathway)) +
  geom_density(show.legend = FALSE) + dens + 
  facet_wrap(~pathway) + 
  theme(strip.text = element_text(size = 6.5), strip.background = element_blank())
#ggsave(p,filename = paste0(figpath,'iraevshealthy.temragated.subsamples.pdf'), width = 5, height = 4.5)



# run gsea on ranks 
# p1 =
#   ggplot(dsample, aes(y = pval, x = NES)) +
#   theme_bw() +
#   ylab('p value') + xlab('Normalized Enrichment Score') +
#   theme(strip.background = element_blank()) +
#   geom_point(size = 0.5, alpha = 0.2, shape = 21) +
#   facet_wrap(~pathway, nrow = 1) +
#   geom_hline(yintercept = 0.05, linetype = 'dashed')
# p1
# ggsave(p1,filename = paste0(figpath,'sample_multipath_irae_no_irae.pdf'), width = 7, height = 2.5)

# enrichment plot 
# original signal 
cont.0 = readRDS(file = here('DE_V4/generated_data/cont0.rds'))
r.0 = ExtractResult(cont.0,what = 'gene.t.ranks', coefficient.number = 1,coef.name = 'baseline_irae')
plist = list()
plist[[1]] = fgsea::plotEnrichment(pathway = s.test$HALLMARK_MTORC1_SIGNALING, stats = r.0$CD8Tcell_TEMRA)$data
plist[[2]] = fgsea::plotEnrichment(pathway = s.test$HALLMARK_MTORC1_SIGNALING, stats = r.irae.noirae$gated)$data
for (i in 1:length(r.irae.subsample)) {
  plist[[i + 2]] = fgsea::plotEnrichment(pathway = s.test$HALLMARK_MTORC1_SIGNALING, stats = r.irae.subsample[[i]])$data 
}
names(plist) = c('CD8 TEMRA Cluster', 'GATED CD8 TEMRA', names(irae_gsea_samples))
pdat = bind_rows(plist, .id = 'id')
pdat$class = ifelse(!pdat$id %in% c("CD8 TEMRA Cluster", "GATED CD8 TEMRA"), 
                    'k-cell samples', pdat$id)
p =
  ggplot(pdat,aes(x = x , y = y)) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 26), axis.text = element_text(size = 14)) + 
  xlab("gene rank") + ylab('Hallmark MTORC1\n Enrichment Score') +
  geom_line(data = pdat %>% filter(class == 'k-cell samples'), size = 0.2, alpha = 0.8, color = 'lightgrey') + 
  geom_line(data = pdat %>% filter(class == 'CD8 TEMRA Cluster'),size = 1.5, color = 'red') + 
  geom_line(data = pdat %>% filter(class == 'GATED CD8 TEMRA'),size = 1.5, color = 'black') 
ggsave(p,filename = paste0(figpath,'downsampled.gsea_irae_no_irae.pdf'), width = 8, height =5.5)

# draw box 
p =
  ggplot(pdat,aes(x = x , y = y)) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 26), axis.text = element_text(size = 14)) + 
  xlab("gene rank") + ylab('Hallmark MTORC1\n Enrichment Score') 
ggsave(p,filename = paste0(figpath,'downsampled.gsea_irae_no_irae_BOX.pdf'), width = 8, height =5.5)
# draw a legend   
gp = ggplot(pdat,aes(x = x , y = y, color = class)) +
  geom_line(size = 5) + scale_color_manual(values = c('red', 'black', 'lightgrey'))+ 
  theme(legend.title = element_blank())
glegend <- cowplot::get_legend(gp)  
pdf(file = paste0(figpath, 'legend.enrichment.pdf'), width = 3, height = 2)
grid::grid.draw(x = glegend)
dev.off()


#
dens =  list(theme_bw(), ggsci::scale_fill_jama(alpha = 0.5) , ggsci::scale_color_jama())

sx.li = LeadingEdgeIndexed(gsea.result.list = irae_gsea_samples, padj.threshold = Inf)
d = do.call(rbind, irae_gsea_samples)
ggplot(d, aes(x = NES, fill = pathway, color = pathway)) +  geom_density() + dens

d2 = do.call(rbind, healthy_gsea_samples)
ggplot(d2, aes(x = NES, fill = pathway, color = pathway)) + geom_density() + dens


# TNF
plist = list()
plist[[1]] = fgsea::plotEnrichment(pathway = s.test$HALLMARK_TNFA_SIGNALING_VIA_NFKB, stats = r.0$CD8Tcell_TEMRA)$data
plist[[2]] = fgsea::plotEnrichment(pathway = s.test$HALLMARK_TNFA_SIGNALING_VIA_NFKB, stats = r.irae.noirae$gated)$data
for (i in 1:length(r.irae.subsample)) {
  plist[[i + 2]] = fgsea::plotEnrichment(pathway = s.test$HALLMARK_TNFA_SIGNALING_VIA_NFKB,
                                         stats = r.irae.subsample[[i]])$data 
}
names(plist) = c('CD8 TEMRA Cluster', 'GATED CD8 TEMRA', names(irae_rank_samples))
pdat = bind_rows(plist, .id = 'id')
pdat$class = ifelse(!pdat$id %in% c("CD8 TEMRA Cluster", "GATED CD8 TEMRA"), 
                    'k-cell samples', pdat$id)
p =
  ggplot(pdat,aes(x = x , y = y)) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 26), axis.text = element_text(size = 14)) + 
  xlab("gene rank") + ylab('Hallmark TNFA \n Enrichment Score') +
  geom_line(data = pdat %>% filter(class == 'CD8 TEMRA Cluster'),size = 3, color = 'black') + 
  geom_line(data = pdat %>% filter(class == 'GATED CD8 TEMRA'),size = 3, color = 'deepskyblue4') +
  geom_line(data = pdat %>% filter(class == 'k-cell samples'), size = 0.2, alpha = 0.5, color = 'deepskyblue3')
p
ggsave(p,filename = paste0(figpath,'TNFAdownsampled.gsea_irae_no_irae.pdf'), width = 8, height =5.5)


# cholesterol
plist = list()
plist[[1]] = fgsea::plotEnrichment(pathway = s.test$HALLMARK_CHOLESTEROL_HOMEOSTASIS, stats = r.0$CD8Tcell_TEMRA)$data
plist[[2]] = fgsea::plotEnrichment(pathway = s.test$HALLMARK_CHOLESTEROL_HOMEOSTASIS, stats = r.irae.noirae$gated)$data
for (i in 1:length(r.irae.subsample)) {
  plist[[i + 2]] = fgsea::plotEnrichment(pathway = s.test$HALLMARK_CHOLESTEROL_HOMEOSTASIS,
                                         stats = r.irae.subsample[[i]])$data 
}
names(plist) = c('CD8 TEMRA Cluster', 'GATED CD8 TEMRA', names(irae_rank_samples))
pdat = bind_rows(plist, .id = 'id')
pdat$class = ifelse(!pdat$id %in% c("CD8 TEMRA Cluster", "GATED CD8 TEMRA"), 
                    'k-cell samples', pdat$id)

pdat %>% head 
p =
  ggplot(pdat,aes(x = x , y = y)) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 26), axis.text = element_text(size = 14)) + 
  xlab("gene rank") + ylab('Hallmark Cholesterol \n Enrichment Score') +
  geom_line(data = pdat %>% filter(class == 'CD8 TEMRA Cluster'),size = 3, color = 'black') + 
  geom_line(data = pdat %>% filter(class == 'GATED CD8 TEMRA'),size = 3, color = 'deepskyblue4') +
  #geom_area(data = pdat %>% filter(class == 'k-cell samples'), size = 0.2, alpha = 0.5, color = 'deepskyblue3')
  geom_line(data = pdat %>% filter(class == 'k-cell samples'), size = 0.2, alpha = 0.5, color = 'deepskyblue3')
ggsave(p,filename = paste0(figpath,'CHOdownsampled.gsea_irae_no_irae.pdf'), width = 8, height =5.5)






# 
# 
# 
# # irae vs healthy donors 
# dsample2 = do.call(rbind, iraehealthy_gsea_samples)
# p1 = 
#   ggplot(dsample2, aes(y = pval, x = NES)) + 
#   theme_bw() +   
#   ylab('p value') + xlab('Normalized Enrichment Score') + 
#   theme(strip.background = element_blank()) + 
#   geom_point(size = 0.5, alpha = 0.2, shape = 21) + 
#   facet_wrap(~pathway, nrow = 1) + 
#   geom_hline(yintercept = 0.05, linetype = 'dashed')
# p1
# 
# ggsave(p1,filename = paste0(figpath,'sample_multipath_irae_no_irae.pdf'), width = 7, height = 2.5)
# plist2 = list()
# plist2[[1]] = fgsea::plotEnrichment(pathway = sig[[8]], stats = rank_irae_healthy$gated)$data
# for (i in 1:length(irae_healthy_samples)) {
#   # make first list element the full data 
#   plist2[[i + 1]] = fgsea::plotEnrichment(pathway = sig[[8]], stats = irae_healthy_samples[[i]]) 
#   plist2[[i + 1]] = plist2[[i + 1]]$data
# }
# names(plist2) = c('full data', names(irae_healthy_samples))
# pdat2 = bind_rows(plist2, .id = 'id')
# p = ggplot(pdat2, aes(x = x , y = y)) + 
#   theme_bw() + 
#   ggtitle('thymoma IRAE group vs healthy donors\ngated TEMRA-like cells') +
#   xlab("gene rank") + ylab('Enrichment Score') + 
#   geom_line(data = pdat2 %>% filter(!id == 'full data'), size = 0.15, alpha = 0.4, color = 'black') + 
#   geom_line(data = pdat2 %>% filter(id == 'full data'), size = 1, color = 'deepskyblue3') 
# ggsave(p,filename = paste0(figpath,'sample_gsea_irae_healthy.pdf'), width = 4, height = 3.5)
# 
# # combined p value plot for both contrasts facetted by module
# dsample$model = 'IRAE vs no IRAE'
# dsample2$model = 'IRAE vs healthy'
# dsample3 = rbind(dsample, dsample2)
# p1 = 
#   ggplot(dsample3, aes(y = pval, x = NES, color = model)) + 
#   theme_bw() +   
#   ylab('k-cell resampled GSEA p value') + xlab('k-cell resampled Normalized Enrichment Score') + 
#   theme(strip.background = element_blank()) + 
#   geom_point(size = 1.5, shape = 16, alpha = 0.2) + 
#   facet_wrap(~pathway, nrow = 1) + 
#   scale_color_manual(values = c("deepskyblue3", 'red')) + 
#   geom_hline(yintercept = 0.05, linetype = 'dashed') + 
#   guides(color = guide_legend(override.aes = list(size=2, alpha = 1))) +
#   theme(legend.position = c(0.88,0.7)) 
# ggsave(p1,filename = paste0(figpath,'sample_multipath_multicontrast.pdf'), width = 7, height = 2.5)

