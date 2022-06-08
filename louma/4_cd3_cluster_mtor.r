# luoma data processing script cluster all cd3 + sorted cells 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(Seurat))
suppressMessages(library(here))
suppressMessages(library(magrittr))
'%ni%' = Negate('%in%')

# save paths 
figpath = here("louma/figures/"); dir.create(figpath)
datapath = here("louma/generated_data/"); dir.create(datapath)

# load louma CD3 object. 
s = readRDS(file = here("data/starting_data/louma/seurat_objects/cd3_louma_Seurat.rds"))

# run simple cluster pipeline 
s = FindVariableFeatures(s, selection.method = "vst", nfeatures = 2000)
s = ScaleData(s, features =  VariableFeatures(s))
s = RunPCA(s, features = VariableFeatures(s), npcs = 50)
# ElbowPlot(s, ndims = 50)
s = FindNeighbors(s, dims = 1:30)
s = FindClusters(s, resolution = 0.5)
s = RunUMAP(s, dims = 1:30)

# add mod score for mtor sigs 
suppressMessages(library(scglmmr))
sig = readRDS(file = here('DE_V3/data_highres/MTOR_SIG_LIST.rds'))
mtor_scores = WeightedCellModuleScore(seurat_object = NULL, gene_matrix = s@assays$RNA@data, 
                                      module_list = sig, threshold = 0, cellwise_scaling = F)

# save a stripped down object with cluster info and metadata 
# this step saves 2GB memory bc internal graphs etc are not saved in the object
s2 = readRDS(file = here("data/louma/seurat_objects/cd3_louma_Seurat.rds"))

# add umap as dr embedding  
umap = as.matrix(s@reductions$umap@cell.embeddings)
s2[["umap"]] <- Seurat::CreateDimReducObject(embeddings = umap, key = "UMAP_", assay = "RNA")

# add cluster info
cmd =
  cbind(seurat_clusters = s@meta.data$seurat_clusters, 
        bc = s@meta.data$barcode_full,
        mtor_scores) %>% 
  remove_rownames() %>% 
  column_to_rownames('bc')
s2 = AddMetaData(object = s2, metadata = cmd)

# save stripped down object wit umap embeddings adn newcluster and mtor score metadata 
saveRDS(object = s2, file = paste0(datapath,'cd3_louma_Seurat_processed.rds'))

