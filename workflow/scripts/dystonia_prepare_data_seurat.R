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
source('~/Desktop/dystonia_snRNAseq_2024/workflow/scripts/dystonia_functions.R')
source('~/Desktop/dystonia_snRNAseq_2024/workflow/scripts/Renvs.R')
source('~/Desktop/dystonia_snRNAseq_2024/workflow/scripts/dystonia_gene_lists.R')
  
## Local Parallelisation
# future::plan()
# future::plan("multiprocess", workers = 3) # depreciated
# future::plan()

## Load Data --------------------------------------------------------------------------
seurat_object <- readRDS(paste0(results_dir, 'prelim/seurat_', region, '.rds'))

# Initial counts and qc plots  --------------------------------------------------------
qc_plot_noFilt <- create_basic_qc_plots(seurat_object) 
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
sce_fcx <- SingleCellExperiment(list(counts = as(seurat_object[["RNA"]]$counts, "dgCMatrix")),
                                colData = seurat_object@meta.data)

# Identify cell outliers and apply filters to sce object
sce_fcx <- get_cell_outliers(sce_fcx, 3, 'higher', 5, 5)
cell_outlier_plot <- create_outlier_plots(sce_fcx, sce_fcx$sum_outlier, sce_fcx$detected_outlier, 
                                    sce_fcx$mito_outlier, sce_fcx$ribo_outlier)
seurat_object <- subset_seurat_object(seurat_object, cell_outliers = sce_fcx$cell_outlier)
cell_outlier_cnts_tbl <- tibble(
  measure = c('umi', 'genes', 'mito', 'ribo', 'total'), 
  count = c(sum(sce_fcx$sum_outlier), sum(sce_fcx$detected_outlier), sum(sce_fcx$mito_outlier),
            sum(sce_fcx$ribo_outlier), sum(sce_fcx$cell_outlier))  
) 

# Post filter counts and QCs
qc_plot_post_Filt <- create_basic_qc_plots(seurat_object) 
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
rm(sce_fcx)

### ISSUE -------------------------------------------------------------------------------
# Creating sketch object using sample_id is too much for local resources for cer and fcx
# Create sketch object and run default Seurat processing 
# Slow takes ages to calculate the Leverage score
# Works for orig ident
# Split by sample ID for downstream analyses
meta_col_orig <- get_meta_col_counts(seurat_object, 'orig.ident')
meta_col_sample <- get_meta_col_counts(seurat_object, 'sample_id')

seurat_object <- split(seurat_object, seurat_object@meta.data[[sample_split]])
seurat_object <- create_sketch_object(seurat_object, pc_thresh)

# Plot QCs
cluster_qc_plot <- create_cluster_qc_plot(seurat_object, pc_thresh, sample_split)

# Cluster counts - Need to change to full dataset first data first
#seurat_sk_fcx$seurat_clusters

# Choose PC threshold
pca_plot <- DimHeatmap(seurat_sk_fcx, dims = 30:pc_thresh, cells = 500, balanced = TRUE)

# Consider using PCA tools - but how to pull out 
#https://github.com/kevinblighe/PCAtools/issues/36
# VarFeatNorm = seurat_sk_str@assays$sketch@layers$scale.data
# colnames(VarFeatNorm) <- colnames(seurat_sk_str)
# p = PCAtools::pca(VarFeatNorm, rank = 30, BSPARAM = IrlbaParam(), metadata = seurat_sk_str@meta.data)
# scree_plot <- PCAtools::screeplot(seurat_sk_str@reductions$pca, axisLabSize = 18, titleLabSize = 22)

# Switch between dataset types
# Analyze the full dataset (on-disk)
# DefaultAssay(obj) <- "RNA"
# # Analyse the sketched dataset (in-memory)
# DefaultAssay(obj) <- "sketch"

## Integration - using harmony
seurat_object <- run_integration(seurat_object)
integration_plot <- create_integration_plot(seurat_object, meta_id = sample_split) 

## This needs to be run after integration --------------------------------
## Not ran this yet - may take a while!!
# project_sketch_data(seurat_object, pc_thresh) 

# Join layers - Need to do this before stacked vln plotting of diff expression
seurat_object <- JoinLayers(seurat_object)

# ### NOTES ON PRELIMINARY CLUST LABELS  ------------------------------------------------
# if (region == 'str') {
#   
#   # Striatum
#   # 13 and 25 experesses GAD 1 and the OLIGS
#   # 22 does not express GAD but exp InN transporters
#   clusters_id_recode <- str_clusters_recode
#   vln_clr_recode <- str_vln_cols_recode
#   umap_clr_recode <- str_umap_cols_recode
# 
# } else if (region == 'cer') {
#   
#   pc_thresh <- 30
#   sample_split <- 'orig.ident'
#   
# } else {
#   
#   # FCX
#   seurat_sk_fcx$harmony_clusters_recode <- fcx_clusters_recode
#   umap_recode <- DimPlot(seurat_sk_fcx, reduction = 'umap.harmony', 
#                          group.by = 'harmony_clusters_recode',
#                          cols = fcx_umap_cols_recode)
#   vln_plot_recode <- create_stacked_vln_plot(seurat_sk_fcx, 'harmony_clusters_recode', 
#                                              fcx_genes, toupper(region), 
#                                              fcx_vln_cols_recode)
#   
#   # Dystonia gene check
#   vln_plot_recode_dyst <- create_stacked_vln_plot(seurat_sk_fcx, 'harmony_clusters_recode', dystonia_genes, 
#                                                   toupper(region), fcx_vln_cols_recode)
#   
# }
# 
# # Plot umap and vlns with cell types and colours specified
# # Recode cell IDs
# seurat_object$harmony_clusters_recode <- clusters_id_recode
# 
# # UMAP with recode cell IDs
# umap_recode <- DimPlot(seurat_object, reduction = 'umap.harmony', 
#                        group.by = 'harmony_clusters_recode',
#                        cols = umap_clr_recode)
# 
# # Region specific gene check vln plot
# vln_plot_recode <- create_stacked_vln_plot(seurat_object, 'harmony_clusters_recode', 
#                                            general_genes, toupper(region), 
#                                            vln_clr_recode)
# # Dystonia gene check vln plot
# vln_plot_recode_dyst <- create_stacked_vln_plot(seurat_object, 'harmony_clusters_recode', 
#                                                 dystonia_genes, toupper(region), 
#                                                 vln_clr_recode)
# 
# 
# ### END NOTES ON PRELIMINARY CLUST LABELS  ------------------------------------------------
# 
# ## Plot gene expression violin plot  ------------------------------------------------------
# # Specific gene check
# vln_plot_specific <- create_stacked_vln_plot(seurat_sk_fcx, 'harmony_clusters', 
#                                              general_genes, toupper(region), 
#                                              fcx_vln_cols_recode)
# 
# seurat_sk_fcx$harmony_clusters_recode <- fcx_clusters_recode
# 
# # Gene check
# vln_plot_recode <- create_stacked_vln_plot(seurat_sk_fcx, 'harmony_clusters_recode',
#                                            general_genes, toupper(region), fcx_vln_cols_recode)
# 
# 
# 
# 
# 
# 
# hip_colours <- c('#76B5C5', '#B200ED', '#FAA0A0', '#EF0029', '#CEE5FD',  
#                  '#95D840FF', "#00BDD2", '#00FF00A5', "#DCBEFF", '#10A53DFF', 
#                  '#6F2DA8', '#ABDBE3', '#1E81B0', '#D2042D', '#006400',   
#                  '#FDE725FF', '#779CBA', '#F58231', '#9A6324')
# 
# str_umap_cols_recode <- c("Str-adult-InN-1" = '#3CBB75FF', "Str-adult-InN-2" = '#31C53F', "Str-adult-InN-3" = '#708238', 
#                           "Str-adult-InN-4" = '#B7FFB7', "Str-adult-InN-5" = '#006400', "Str-OPC-1" = '#FDE725FF', 
#                           "Str-adult-InN-6" = '#95D840FF', "Str-adult-InN-7" = '#2FF18B', "Str-adult-InN-8" = '#9DC183', 
#                           "Str-adult-Ast-1" = '#FF5959', "Str-adult-InN-9" = '#3CBB75FF', "Str-adult-InN-10" = '#31C53F', 
#                           "Str-adult-InN-11" = '#708238', "Str-OPC-2" = '#FDE725FF', "Str-adult-InN-12" = '#006400', 
#                           "Str-adult-InN-13" = '#95D840FF', "Str-adult-ExN-1" = '#00B6EB', "Str-adult-MG" = '#F58231', 
#                           "Str-adult-InN-14" = '#2FF18B', "Str-adult-InN-15" = '#9DC183', "Str-adult-ExN-2" = '#00B6EB', 
#                           "Str-adult-InN-16" = '#0CB702', "Str-adult-InN-17" = '#00BE67', "Str-adult-Ast-2" = '#FF5959', 
#                           "Str-adult-InN-18" = '#7CAE00', "Str-OPC-3" = '#FDE725FF')
# 
# str_vln_cols_recode <- c("Str-adult-InN-1" = '#3CBB75FF', "Str-adult-InN-2" = '#3CBB75FF', "Str-adult-InN-3" = '#3CBB75FF', 
#                          "Str-adult-InN-4" = '#3CBB75FF', "Str-adult-InN-5" = '#3CBB75FF', "Str-OPC-1" = '#FDE725FF', 
#                          "Str-adult-InN-6" = '#3CBB75FF', "Str-adult-InN-7" = '#3CBB75FF', "Str-adult-InN-8" = '#3CBB75FF', 
#                          "Str-adult-Ast-1" = '#FF5959', "Str-adult-InN-9" = '#3CBB75FF', "Str-adult-InN-10" = '#3CBB75FF', 
#                          "Str-adult-InN-11" = '#3CBB75FF', "Str-OPC-2" = '#FDE725FF', "Str-adult-InN-12" = '#3CBB75FF', 
#                          "Str-adult-InN-13" = '#3CBB75FF', "Str-adult-ExN-1" = '#00B6EB', "Str-adult-MG" = '#F58231', 
#                          "Str-adult-InN-14" = '#3CBB75FF', "Str-adult-InN-15" = '#3CBB75FF', "Str-adult-ExN-2" = '#00B6EB', 
#                          "Str-adult-InN-16" = '#3CBB75FF', "Str-adult-InN-17" = '#3CBB75FF', "Str-adult-Ast-2" = '#FF5959', 
#                          "Str-adult-InN-18" = '#3CBB75FF', "Str-OPC-3" = '#FDE725FF')
# 
# 
# 
# # Calculate aggr and aver expr for regional comparison of dystonia gene expression
# aver_exp_mat <- calculate_average_expression(seurat_object, paste0(region, '_adult'), dystonia_genes)
# aggr_exp_mat <- calculate_aggregated_expression(seurat_object, paste0(region, '_adult'), dystonia_genes)
# 
# saveRDS(object = aver_exp_mat, file = paste0(R_dir, "seurat_aver_exp_", region, ".Rds"))
# saveRDS(object = aggr_exp_mat, file = paste0(R_dir, "seurat_aggr_exp_", region, ".Rds"))
# 
# test <- readRDS(paste0(R_dir, "seurat_aggr_exp_", region, ".Rds"))
## Create markdown doc  ---------------------------------------------------------------
# rmarkdown::render(markdown_doc, output_file = markdown_html, output_dir = R_dir)

# Save Seurat object
saveRDS(seurat_object, paste0(R_dir, 'seurat_', region, '_basic.rds'))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

