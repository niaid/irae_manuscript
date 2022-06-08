### cluster markers from original paper 
cd8gene = c("KIR2DL4", "TRDC", "ANXA1", "IL7R","GZMK","EOMES", "SLC4A10",
            "KLRB1", "MIAT", "RNF213","CCR7", "KLI-2","UCP2", "GZMB",
            "TYMS", "MKI67", "TFC1" ,"TCF7", "ITK", "CD8A", "CD3E")
cd4gene = c("IL4I1", "IL23R", "NR4A1", "EGR1", "CCL5", "ANXA1", "S1PR1", "KLF2", 
            "CCR7", "TOX2", "DRAIC1","DRAIC", "IKZF2", "FOXP3", "GZMH", "GBP5", "GPR25", "CD4")
genecombined = c(cd4gene,cd8gene) %>% unique()