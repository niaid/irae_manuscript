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

le1 = GetLeadingEdgeGenes(gsea.result.list = test,celltype.index = 16,module.name = 'HALLMARK_MTORC1_SIGNALING')
le2 = GetLeadingEdgeGenes(gsea.result.list = test,celltype.index = 14,module.name = 'HALLMARK_MTORC1_SIGNALING')
le3 = GetLeadingEdgeGenes(gsea.result.list = test,celltype.index = 16,module.name = 'HALLMARK_TNFA_SIGNALING_VIA_NFKB')
le4 = GetLeadingEdgeGenes(gsea.result.list = test,celltype.index = 15,module.name = 'HALLMARK_INFLAMMATORY_RESPONSE')


genes = list(
  'mDC mTOR' = le1, 
  'TEMRA mTOR' = le2, 
  'mDC TNFA' = le3, 
  'CD1c DC Inflammatory' = le4
)


temp <- 
  VennDiagram::venn.diagram(genes, 
                            # color 
                            # fill = pals::brewer.seqseq2(4),
                            fill = c('#f7931e', '#22b573', '#f7931e', '#0071bc'),
                            # text
                            cat.cex = 0.5, cex = 2.5,
                            cex.fontfamily ="Helvetica",
                            cat.default.pos = "outer",
                            cat.fontfamily = "Helvetica",
                            # file 
                            imagetype="png", 
                            height = 900, 
                            width = 900, 
                            resolution = 1200, 
                            compression = "lzw",
                            filename = NULL
                            # paste0(figpath,"mtor_signaturesvenn.png")
  )
pdf(file = paste0(figpath,"correlated.signatures.pdf"), width = 3.9, height = 4)
grid::grid.draw(temp)
dev.off()



genes.2 = list(
  #'mDC mTOR' = le1, 
  'TEMRA mTOR' = le2, 
  #'mDC TNFA' = le3, 
  'CD1c DC Inflammatory' = le4
)
## 
temp <- 
  VennDiagram::venn.diagram(genes.2, 
                            # color 
                            # fill = pals::brewer.seqseq2(4),
                            fill = c('#f7931e', '#22b573'),
                            # text
                            cat.cex = 0.5, cex = 2.5,
                            cex.fontfamily ="Helvetica",
                            cat.default.pos = "outer",
                            cat.fontfamily = "Helvetica",
                            # file 
                            imagetype="png", 
                            height = 900, 
                            width = 900, 
                            resolution = 1200, 
                            compression = "lzw",
                            filename = NULL
                            # paste0(figpath,"mtor_signaturesvenn.png")
  )
pdf(file = paste0(figpath,"correlated.signatures.sub.pdf"), width = 3.9, height = 4)
grid::grid.draw(temp)
dev.off()


