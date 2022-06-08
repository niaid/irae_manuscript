suppressMessages(library(here))
suppressMessages(library(grDevices))
figpath = here("cd8/mtor_cd8/figuresV4//")

# read mtor signatures with the thymoma irae temra mtor leading edge 
sig = readRDS(file = here('DE_V3/data_highres/MTOR_SIG_LIST.rds'))

# subset of curated mtor signatures to plot 
sigsub = sig[c(8,7, 3, 5)]
names(sigsub) = c("TM Sig","Reactome", "KEGG", "Biocarta")
# draw set intersection of gene membership in signatures. 
temp <- 
  VennDiagram::venn.diagram(sigsub, 
                          # color 
                          fill = pals::cols25(4),
                          #fill = c('#7f8d93', '#ffe18e', '#eb7160', 'red'),
                          # text
                          cat.cex = 1, cex = 2.5, cex.fontfamily ="Helvetica",
                          cat.default.pos = "inner",
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
pdf(file = paste0(figpath,"mtor_signaturesvenn.pdf"), width = 3.9, height = 4)
grid::grid.draw(temp)
dev.off()





