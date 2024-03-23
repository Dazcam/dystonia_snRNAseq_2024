#--------------------------------------------------------------------------------------
#
#    Dystonia - prepare adult brain scRNAseq data
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Stiletti data: https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443
#  Note that the Stiletti data has counts stored in seurat_obj[["RNA"]]$data
#  Using Seurat 5 primarily, but also bioconductor packages for QC 

##  Load Packages, functions and variables  -------------------------------------------
message('Setting environment variables ...')
if (Sys.info()[["nodename"]] == "Darrens-iMac-2.local") {
  
  library(yaml)
  root_dir <- '~/Desktop/dystonia_snRNAseq_2024/'
  yaml_file <- yaml.load_file(paste0(root_dir, 'config/config.yaml'))
  region <- yaml.load(yaml_file$region)

  source(paste0(root_dir, 'workflow/scripts/dystonia_functions.R'))
  source(paste0(root_dir, 'workflow/scripts/dystonia_Renvs.R'))
  source(paste0(root_dir, 'workflow/scripts/dystonia_gene_lists.R'))

} else {

  source('scripts/dystonia_functions.R')
  source('scripts/dystonia_Renvs.R')
  source('scripts/dystonia_gene_lists.R')
  
}

## Load Data --------------------------------------------------------------------------
seurat_object <- readRDS(paste0(R_dir, '01seurat_', region, '.rds'))

# Initial counts and qc plots  --------------------------------------------------------
qc_plot_noFilt <- create_basic_qc_plots(seurat_object, meta_id = sample_split) 
counts <- tibble('Cells' = c(ncol(seurat_object)),
                 'Genes' = c(nrow(seurat_object)))
dist_counts_noFilt <- rbind(summary(seurat_object$nCount_RNA),
                  summary(seurat_object$nFeature_RNA)) %>%
  t() %>%
  as_tibble(rownames = 'measure') %>%
  dplyr::rename(cells = V1,
                genes = V2)

### Filtering  -----------------------------------------------------------
# Join layers to apply filters
seurat_object <- JoinLayers(seurat_object)

# Convert to single cell experiment
sce_obj <- SingleCellExperiment(list(counts = as(seurat_object[["RNA"]]$counts, "dgCMatrix")),
                                colData = seurat_object@meta.data)

# Identify cell outliers and apply filters to sce object
sce_obj <- get_cell_outliers(sce_obj, 3, 'higher', 5, 5)
cell_outlier_plot <- create_outlier_plots(sce_obj, sce_obj$sum_outlier, sce_obj$detected_outlier, 
                                    sce_obj$mito_outlier, sce_obj$ribo_outlier)
seurat_object <- subset_seurat_object(seurat_object, cell_outliers = sce_obj$cell_outlier)
cell_outlier_cnts_tbl <- tibble(
  measure = c('umi', 'genes', 'mito', 'ribo', 'total'), 
  count = c(sum(sce_obj$sum_outlier), sum(sce_obj$detected_outlier), sum(sce_obj$mito_outlier),
            sum(sce_obj$ribo_outlier), sum(sce_obj$cell_outlier))  
) 

# Post filter counts and QCs
qc_plot_post_Filt <- create_basic_qc_plots(seurat_object, meta_id = sample_split) 
counts_post_Filt <- rbind(counts,
                         tibble('Cells' = c(ncol(seurat_object)),
                          'Genes' = c(nrow(seurat_object)))) %>%
  mutate(Measure = c('Before', 'After')) %>%
  relocate(Measure)
dist_counts_post_Filt <- rbind(summary(seurat_object$nCount_RNA),
                            summary(seurat_object$nFeature_RNA)) %>%
  t() %>%
  as_tibble(rownames = 'measure') %>%
  dplyr::rename(cells = V1,
                genes = V2) %>%
  inner_join(dist_counts_noFilt, by = 'measure', suffix = c('_before', '_after')) %>%
  relocate(measure, cells_before, cells_after, genes_before, genes_after)

# Remove sce object to save space 
rm(sce_obj)

## Create sketch object  --------------------------------------------------------------
# Sketching is not currently compatible with SCT 
# normalisation: https://github.com/satijalab/seurat/issues/7336
meta_col_orig_cnts <- get_meta_col_counts(seurat_object, 'orig.ident')
meta_col_sample_cnts <- get_meta_col_counts(seurat_object, 'sample_id')
seurat_object <- split(seurat_object, seurat_object@meta.data[[sample_split]]) # Split first by lane
seurat_object <- create_sketch_object(seurat_object, dims = pc_thresh, norm_method = 'log_sketch',
                                      resolution = resolution_set, cell_num = sketch_num)

# Check PCs
pca_plot <- DimHeatmap(seurat_object, dims = 30:pc_thresh, cells = 500, balanced = TRUE)

# Create resolution plots for comparison
res_plotlist <- create_resolution_plotlist(seurat_object, resolution = resolution_set,
                                           meta_id = sample_split)

# Plot QCs
cluster_qc_plot <- create_cluster_qc_plot(seurat_object, pc_thresh, sample_split)
meta_col_cluster_cnts <- get_meta_col_counts(seurat_object, 'seurat_clusters')

## Integration - using Harmony
seurat_object <- run_integration(seurat_object, 'harmony', pc_thresh, resolution_set)
integration_plotlist <- create_integration_plotlist(seurat_object, 
                                                    meta_id = sample_split, 
                                                    dims = pc_thresh,
                                                    reduction = resolution_set) 

## Project sketch data to entire object - still having issues with this:
seurat_object <- project_sketch_data(seurat_object,
                                     pc_thresh,
                                     'harmony',
                                     'umap.harmony',
                                     'harmony_clusters_0.1')

DefaultAssay(seurat_object) <- "sketch"

## Create markdown doc  ---------------------------------------------------------------
rmarkdown::render(markdown_prep_doc, output_file = markdown_prep_html, output_dir = R_dir)

# Join layers - Required for vln plotting or diff expression - FCX crashes locally
seurat_object <- JoinLayers(seurat_object)

# Save Seurat object
message('\nWriting Seurat object ...\n')
saveRDS(seurat_object, paste0(R_dir, '02seurat_', region, '.rds'))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

