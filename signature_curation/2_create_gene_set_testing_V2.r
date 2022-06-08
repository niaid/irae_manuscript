# test specific hypothesis related to IRAE biomarkers and plausable immune related processes 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(here))


# load BTM and use select BTM pathways 
btm = readRDS(file = "signature_curation/BTM_li.rds")
test_btm = btm[c("LI.M21 cell adhesion (lymphocyte homing)", 
                 "LI.M30 cell movement, Adhesion & Platelet activation",
                 "LI.M51 cell adhesion", 
                 "LI.M84 integrins and cell adhesion", "LI.M91 adhesion and migration, chemotaxis", 
                 "LI.M117 cell adhesion (GO)" , 
                 "LI.M133.1 cell cell adhesion", 
                 "LI.M127 type I interferon response",     
                 "LI.M158.0 interferon alpha response (I)" ,    
                 "LI.M158.1 interferon alpha response (II)"  ,  
                 "LI.M111.0 viral sensing & immunity; IRF2 targets network (I)" ,            
                 "LI.M111.1 viral sensing & immunity; IRF2 targets network (II)", 
                 # state generic 
                 "LI.M167 enriched in cell cycle" ,
                 "LI.M160 leukocyte differentiation" ,
                 "LI.M7.1 T cell activation (I)", 
                 "LI.M52 T cell activation (IV)",
                 "LI.M204.0 chaperonin mediated protein folding (I)",
                 # state specific 
                 "LI.M4.12 C-MYC transcriptional network",
                 "LI.M20 AP-1 transcription factor network", 
                 "LI.M206 Wnt signaling pathway",
                 "LI.M18 T cell differentiation via ITK and PKC", 
                 "LI.M19 T cell differentiation (Th2)",
                 "LI.M36 T cell surface, activation",
                 "LI.M165 enriched in activated dendritic cells (II)",
                 "LI.S6 CD4 T cell surface signature Th1-stimulated",
                 "LI.M4.6 cell division in stimulated CD4 T cells",
                 # other 
                 "LI.M209 lysosome",
                 "LI.M164 xenobiotic metabolism", 
                 "LI.M225 metabolism of steroids")] 

# load reactome and use select reactome pathways 
reactome = readRDS(file = "signature_curation/reactome.rds")
test_reactome = reactome[c("REACTOME_TRANSLOCATION_OF_ZAP_70_TO_IMMUNOLOGICAL_SYNAPSE", 
                           "REACTOME_IL_7_SIGNALING",
                           "REACTOME_PD1_SIGNALING", 
                           "REACTOME_PECAM1_INTERACTIONS", 
                           "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION", 
                           "REACTOME_ANTIGEN_PROCESSING_UBIQUITINATION_PROTEASOME_DEGRADATION")]

names(test_reactome) = c("Reactome ZAP70 immunological synapse", 
                         "Reactome IL7 signaling", 
                         "Reactome PD1 signaling", 
                         "Reactome PECAM1 interactions",
                         "Reactome Cross presentation",
                         "Reactome antigen processing and proteasome")


# baseline signatuers we described in Kotliarov et al 2020. 
nmed = readRDS("signature_curation/Kotliarov_2020/baseline_natmed_signatures.rds")


# irae associated genes from shahabi et. al. 
# https://pubmed.ncbi.nlm.nih.gov/23521917/
hypset1 = read_delim("signature_curation/Shahabi_2013/shahabi_0.txt", delim = "\t")$Shahabi_0 %>% as.vector()
hypset2 = read_delim("signature_curation/Shahabi_2013/shahabi_1.txt", delim = "\t")$Shahabi_1 %>% as.vector()
hyp_test = list(hypset1, hypset2) 
names(hyp_test) = c("IRAE_Ipilumimab_Baseline_shahabi", "IRAE_Ipilumimab_3wk_shahbai")
                    

### Data from Yao et. al. https://pubmed.ncbi.nlm.nih.gov/31209400/
prog_sig = data.table::fread("signature_curation/yao_ni_2019/progenitor_like_sig_yao_NI_2019_tox_Stable2.txt")
exh = data.table::fread("signature_curation/yao_ni_2019/terminally_exhausted_sig_yao_NI_2019_tox_Stable2.txt")
mem_sig = data.table::fread("signature_curation/yao_ni_2019/memory_precursor_sig_yao_NI_2019_tox_Stable2.txt")
yao_wherry = list("Progenitor-like Yao Nat. Imm. 2019" = prog_sig,
           "Exhausted Yao Nat. Imm. 2019" = exh,
           "Memory-precursor Yao Nat. Imm. 2019" = mem_sig)

convertMouseGeneList <- function(x){
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = biomaRt::getLDS(attributes = c("mgi_symbol"), 
                            filters = "mgi_symbol", 
                            values = x ,
                            mart = mouse, attributesL = c("hgnc_symbol"),
                            martL = human, uniqueRows=TRUE)
  return(genesV2)
}
yao_human = lapply(yao, function(x) convertMouseGeneList(x = x$`Gene Symbol`))
yao_human = lapply(yao_human, function(x) x$HGNC.symbol)

# combine signatures and save 
hyp_test1 = c(hyp_test, test_btm, test_reactome, nmed, yao_human) 
saveRDS(hyp_test1, file = "signature_curation/hypothesis_set_1_V2.rds")

sessionInfo()
# R version 3.6.1 (2019-07-05)
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] magrittr_2.0.2  here_0.1        forcats_0.5.0   stringr_1.4.0   dplyr_1.0.3     purrr_0.3.3     readr_1.3.1    
# [8] tidyr_1.0.2     tibble_3.1.0    ggplot2_3.3.0   tidyverse_1.3.0
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.4           lubridate_1.7.4      prettyunits_1.1.1    rprojroot_1.3-2      assertthat_0.2.1    
# [6] digest_0.6.27        utf8_1.1.4           BiocFileCache_1.10.2 R6_2.4.1             cellranger_1.1.0    
# [11] backports_1.2.1      reprex_0.3.0         stats4_3.6.1         RSQLite_2.2.0        httr_1.4.1          
# [16] pillar_1.5.0         rlang_0.4.10         progress_1.2.2       curl_4.3             readxl_1.3.1        
# [21] rstudioapi_0.11      data.table_1.12.8    blob_1.2.1           S4Vectors_0.24.3     bit_1.1-15.2        
# [26] biomaRt_2.42.0       munsell_0.5.0        broom_0.7.4          compiler_3.6.1       modelr_0.1.6        
# [31] pkgconfig_2.0.3      askpass_1.1          BiocGenerics_0.32.0  openssl_1.4.1        tidyselect_1.1.0    
# [36] IRanges_2.20.2       XML_3.99-0.3         fansi_0.4.2          withr_2.4.0          crayon_1.4.1        
# [41] dbplyr_1.4.2         rappdirs_0.3.1       grid_3.6.1           jsonlite_1.6.1       gtable_0.3.0        
# [46] lifecycle_1.0.0      DBI_1.1.0            scales_1.1.0         cli_2.3.1            stringi_1.4.6       
# [51] fs_1.3.2             xml2_1.3.2           ellipsis_0.3.1       generics_0.0.2       vctrs_0.3.6         
# [56] tools_3.6.1          bit64_0.9-7          Biobase_2.46.0       glue_1.4.2           hms_0.5.3           
# [61] parallel_3.6.1       AnnotationDbi_1.48.0 colorspace_1.4-1     rvest_0.3.5          memoise_1.1.0       
