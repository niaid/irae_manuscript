suppressMessages(library(scglmmr))
suppressMessages(library(here))

# set paths 
datapath = here('louma/generated_data_scglmmr/'); dir.create(datapath)
figpath = here('louma/figures_scglmmr/'); dir.create(figpath)

# load data 
s = readRDS(file = here('data/louma/seurat_objects/cd3_louma_Seurat.rds'))
meta = readRDS(file = here('louma/cd3_clustering/generated_data/louma_clustermetadata.rds'))
meta$sampleid = meta$subjectid

# define counts and metadata and subset to cells above rm seurat object from workspace 
meta$seurat_clusters = paste("C", meta$seurat_clusters, sep = "")
umi = s@assays$RNA@counts

# QC contingency of cells by subject for each celltype 
tab = scglmmr::SubjectCelltypeTable(metadata = meta, celltype_column = "seurat_clusters", sample_column = "sampleid")
tab$celltypes_remove; tab$`low representation celltypes`; tab$table

# remove cells prior to pseudobulk analysis 
meta = meta[!meta$seurat_clusters %in% tab$celltypes_remove, ]

# subset data 
umi = umi[ ,rownames(meta)]

# pseudobulk workflow 
pb = scglmmr::PseudobulkList(rawcounts = umi, metadata = meta, sample_col = "sampleid",
                             celltype_col = "seurat_clusters", avg_or_sum = "sum")
designmat = scglmmr::BulkDesignMatrix(metadata = meta, sample_column = "sampleid",
                                      variable_column = "IRAE", pseudobulklist = pb)
dge = scglmmr::NormalizePseudobulk(pseudobulklist = pb, design_matrix = designmat,
                                   minimum.gene.count = 5)


# qc model rank and detect unused factor levels. 
stopifnot(Matrix::rankMatrix(designmat) == ncol(designmat))
stopifnot(any(colSums(designmat) == 0) == FALSE)


# custom a priori contrasts
c_mat = makeContrasts(
  colitis_vs_nocolitis = (IRAEcolitis - IRAEno), 
  colitis_vs_healthy = IRAEcolitis - IRAEcontrol, 
  nocolitis_vs_healthy = IRAEno - IRAEcontrol, 
  #control_adjusted_colitis_vs_nocolitis = (IRAEcolitis - IRAEcontrol) - (IRAEno - IRAEcontrol),
  levels = colnames(designmat)
)

# fit mixed model for the multi timepoint contrasts 
fit = RunVoomLimma(dgelists = dge, design_matrix = designmat, do_contrast_fit = TRUE, my_contrast_matrix = c_mat)


# make average matrix for visualization
av = scglmmr::PseudobulkList(rawcounts = umi, metadata = meta, sample_col = "sampleid",
                             celltype_col = "seurat_clusters", avg_or_sum = "average")

# save 
saveRDS(fit, file = paste0(datapath, "fit_object.rds"))
saveRDS(pb, file = paste0(datapath, "pb_object.rds"))
saveRDS(av, file = paste0(datapath, "av_object.rds"))
saveRDS(dge, file = paste0(datapath, "dge_object.rds"))
saveRDS(meta, file = paste0(datapath, "meta_object.rds"))

