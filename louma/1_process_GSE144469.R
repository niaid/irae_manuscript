# mtor  analysis 
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(here))

datapath = here("data/louma/seurat_objects/"); dir.create(datapath)

# read in all datasets from Louma 
cd3_files = list.files(path = here("data/louma/GSE144469_RAW/cd3/"), pattern = "CD3", full.names = TRUE) %>% as.list()
cd3_names = list.files(path = here("data/louma/GSE144469_RAW/cd3/"), pattern = "CD3", full.names = FALSE)
cd3 = lapply(cd3_files, function(x){ Read10X(data.dir = x, strip.suffix = TRUE) })
for (i in 1:length(cd3)) {
  colnames(cd3[[i]]) = paste(colnames(cd3[[i]]),cd3_names[i], sep = "_") 
}
cd3_m = do.call(cbind, cd3)

# create cell metadata 
cell_met = colnames(cd3_m) %>% 
  as.data.frame() %>% 
  rename(barcode_full = '.') %>% 
  mutate(sample = str_sub(barcode_full, -7, -1)) %>% 
  mutate(sample = if_else(str_sub(sample,1,1) == "_", true = str_sub(sample,2,7 ),false = sample ))

# add metadata based on map values list from GSE file  
meta = read_delim(file = here("louma/louma_metadata.txt"), delim = "\t")
meta = meta %>% mutate(subject_sample = paste(subjectid,sampleid, sep = "-"))
idsmap = colnames(meta)
met = list()
for (i in 1:length(idsmap)) {
  md = cell_met
  md[[idsmap[i]]] = plyr::mapvalues(md$sample, from = meta$subject_sample, to = meta[[idsmap[i]]] )
  met[[i]] = md %>% select(idsmap[i])
} 
met = do.call(cbind, met) %>% as.data.frame() 
met = cbind(cell_met, met) 
met = met %>% mutate(barcode = barcode_full) %>% 
  column_to_rownames("barcode")

# create a single merged CD3 object and save 
s = CreateSeuratObject(counts = cd3_m, meta.data = met, assay = "RNA", min.cells = 10, min.features = 200)
saveRDS(s, file = paste0(datapath, "cd3_louma_Seurat.rds"))

