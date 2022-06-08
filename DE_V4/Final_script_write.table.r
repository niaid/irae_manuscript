# make tables 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))

datapath = file.path(here('DE_V4/generated_data/tables/')); dir.create(datapath)

# specify order of variables in the output table for readability 
var.order = c('contrast', 'celltype', 'pathway', 'NES',
              'padj', 'av_jaccard', 'leadingEdge')

format.result = function(x) { 
  x %>% 
    select(all_of(var.order), everything()) %>%
    arrange(celltype, desc(NES)) %>% 
    tibble::remove_rownames()
  }

# Baseline temporally stable 
# filter to only include the temporally stable pathways 
g0.list.combined = readRDS(file = here('DE_V4/generated_data/temp/g0.list.combined.rds'))
g0.list.combined = lapply(g0.list.combined, function(x) x %>% filter(padj < 0.01))
li.g0 = LeadingEdgeIndexed(g0.list.combined, padj.threshold = 0.01)
d0.res = EnrichmentJaccard(gsealist = g0.list.combined, 
                           indexedgenes = li.g0, 
                           saveplot = FALSE, FALSE)
d0.res$contrast = 'baseline_irae'
d0.res = format.result(d0.res)

  


#treatment
g1.list.combined = readRDS(file = here('DE_V4/generated_data/temp/g1.list.combined.rds'))
g1.list.combined = lapply(g1.list.combined, function(x) x %>% filter(padj < 0.01))
li.g1 = LeadingEdgeIndexed(g1.list.combined, padj.threshold = 0.01)
d1.res = EnrichmentJaccard(gsealist = g1.list.combined,
                           indexedgenes = li.g1,
                           saveplot = FALSE, FALSE)
d1.res$contrast = 'treatment'
d1.res = format.result(d1.res)



# treatment_delta  
gd.list.combined = readRDS(file = here('DE_V4/generated_data/temp/gd.list.combined.rds'))
gd.list.combined = lapply(gd.list.combined, function(x) x %>% filter(padj < 0.05))
li.gd = LeadingEdgeIndexed(gd.list.combined, padj.threshold = 0.05)
gd.res = EnrichmentJaccard(gsealist = gd.list.combined,
                           indexedgenes = li.gd,
                           saveplot = FALSE, FALSE)
gd.res$contrast = 'treatment_delta'
gd.res = format.result(gd.res)

# combine 
d = rbind(d1.res, gd.res, d0.res)

data.table::fwrite(d,
                   file = paste0(datapath,'irae.combined.pseudobulk.results.txt'), 
                   sep = '\t')


