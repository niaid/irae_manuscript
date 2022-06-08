# make gene set table 
library(here)
library(magrittr)
hyp = readRDS(file = here('signature_curation/hypothesis_set_1_V2.rds'))
mtx = sapply(hyp, "length<-", max(lengths(hyp)))
data.table::fwrite(mtx, file = paste0(here('git_ignore/hypset.txt')),sep = '\t')
