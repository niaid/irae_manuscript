### read bulk results and data and do downsampling 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(scglmmr))
figpath = here("cd8/mg_cd8/figuresV4/")
datapath = here("cd8/mg_cd8/generated_dataV4/")


# read sig 
csig = readRDS('signature_curation/combined_signatures.rds')
hlmk = readRDS('signature_curation/hallmark.rds')
csig = c(csig, hlmk)
mt.sig = csig[names(csig)[grepl('MTOR', names(csig))]]

biocarta.sig = list('Biocarta MTOR' = dplyr::union(mt.sig$BIOCARTA_IGF1MTOR_PATHWAY, mt.sig$BIOCARTA_MTOR_PATHWAY))
kegg.sig = list('KEGG MTOR' = mt.sig$KEGG_MTOR_SIGNALING_PATHWAY)
reactome.sig =  list('Reactome MTOR' = dplyr::union(mt.sig$REACTOME_ENERGY_DEPENDENT_REGULATION_OF_MTOR_BY_LKB1_AMPK,
                                                   mt.sig$REACTOME_MTORC1_MEDIATED_SIGNALLING))

mtor.sigs = c(biocarta.sig, kegg.sig, reactome.sig)
saveRDS(mtor.sigs, file = paste0(datapath, 'mtor.sigs.rds'))

g0.list.combined = readRDS(file = here('DE_V4/generated_data/temp/g0.list.combined.rds'))
temra.mtor = g0.list.combined$CD8Tcell_TEMRA %>% filter(pathway == 'HALLMARK_MTORC1_SIGNALING')
temra.mtor =  temra.mtor$leadingEdge
names(temra.mtor) = 'TEMRA MTOR'
temra.mtor

# add temramtor 
mtor.sigs2 = c(biocarta.sig, kegg.sig, reactome.sig, temra.mtor)
saveRDS(mtor.sigs2, file = paste0(datapath, 'mtor.sigs2.rds'))


# make upset plot 
ml = UpSetR::fromList( mtor.sigs2 )
pdf(file = paste0(figpath, 'mtorupset.pdf'), width = 5, height = 3.5)
UpSetR::upset(data = ml) 
dev.off()

# save genes 
d = lapply(mtor.sigs2, as.data.frame) 
d = do.call(rbind, d) %>% tibble::rownames_to_column('sig')
data.table::fwrite(d,file = paste0(datapath,'d.txt'),sep = '\t')


