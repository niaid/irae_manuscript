# visulaziaiton of gsea contrast results 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))

####################################
# derive temporally stable gsea for hallmark and hypothesis set 
####################################

figpath.temp = file.path(here('DE_V4/figures/temp/')); dir.create(figpath.temp)
datapath = file.path(here('DE_V4/generated_data/temp/')); dir.create(datapath)

# cobinee baseline hallmark and hypothesis set gene set enrichment. 
hyp.0 = readRDS(here("DE_V4/generated_data/hyp.0.rds"))
hlmk.0 = readRDS(here("DE_V4/generated_data/hlmk.0.rds"))
g0.list.combined = Map(f = rbind, hyp.0, hlmk.0)
saveRDS(g0.list.combined, file = paste0(datapath,'g0.list.combined.rds'))


# repeat for treatment contrast 
hyp.treat = readRDS(here("DE_V4/generated_data/hyp.treat.rds"))
hlmk.treat = readRDS(here("DE_V4/generated_data/hlmk.treat.rds"))
g1.list.combined = Map(f = rbind, hyp.treat, hlmk.treat)
saveRDS(g1.list.combined, file = paste0(datapath,'g1.list.combined.rds'))

# repeat for treatment delta contrast 
hyp.delta = readRDS(here("DE_V4/generated_data/hyp.delta.rds"))
hlmk.delta = readRDS(here("DE_V4/generated_data/hlmk.delta.rds"))
gd.list.combined = Map(f = rbind, hyp.delta, hlmk.delta)
saveRDS(gd.list.combined, file = paste0(datapath,'gd.list.combined.rds'))

# retain temporally stable enrichments 
# subtract treatment and delta enrichments from baseline enrichments 
# define G0(j, i) in each cell type i for a given baseline G0 enrichmed pathway j at adjusted p value 0.01: 
#    if the same signal enriched at adjusted p 0.1 (0.1 not 0.01) in the treatment Gt(j,i) or delta contrast Gd(j,i) 
#    remove this signal  from the baseline enr 
#    The remaining modules are more consistent with temporal stability. 
#    G0stable = G0 - ( intersect(G0, Gt) + intersect(G0, Gd) ) 
d.0 = RbindGseaResultList(g0.list.combined, padj_filter = 0.01) %>% mutate(contrast = 'baseline')
d.treat = RbindGseaResultList(g1.list.combined, padj_filter = 0.1) %>% mutate(contrast = 'treatment')
d.delta = RbindGseaResultList(gd.list.combined, padj_filter = 0.1) %>% mutate(contrast = 'treatment_delta')
dl = list(d.0, d.treat, d.delta)
dl = lapply(dl, function(x){ x %>% mutate(signal = paste(pathway, celltype, sep = '~')) })
d = do.call(rbind, dl)
G.global = table(d$signal) 
non.unique = G.global[G.global>1] %>% names()

# curated baseline enrichment based on tamporal stsbility 
g0.li = LeadingEdgeIndexed(gsea.result.list = g0.list.combined, padj.threshold = 0.01)
g0.list.sub = lapply(g0.list.combined, function(x) x %>% filter(padj < 0.01))
g0.JI = EnrichmentJaccard(gsealist = g0.list.sub, indexedgenes = g0.li, 
                          saveplot = TRUE, returnJaccardMtx = TRUE,
                          figpath = figpath.temp)
g0.result.sort = g0.JI$sortedgsea %>%
  mutate(signal = paste(celltype, pathway, sep = '~')) %>% 
  filter(!signal %in% non.unique) # filter out non temporal stable signals 
data.table::fwrite(g0.result.sort, file = paste0(datapath, 'g0.result.sort.txt'),sep = "\t")



