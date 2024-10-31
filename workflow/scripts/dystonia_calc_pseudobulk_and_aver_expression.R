#--------------------------------------------------------------------------------------
#
#    Dystonia - Get aggregate (pseudobulk) and average expression
#
#--------------------------------------------------------------------------------------

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
seurat_object <- readRDS(paste0(R_dir, '03seurat_', region, '.rds'))
# fetal_region <- str_split_i(fetal_region, '_', 1)
# seurat_fetal <- readRDS(paste0(fetal_dir, 'seurat_', fetal_region, '_fetal.rds'))
# Idents(seurat_fetal) <- 'cellIDs'

# Switch to whole dataset
message('\nChanging to RNA object ...\n')
DefaultAssay(seurat_object) <- 'RNA'
seurat_object <- JoinLayers(seurat_object) # Do this for sketch and RNA independently
message('Any NAs in Idents: ', anyNA(Idents(seurat_object)))

# # Recode cluster IDs - sketch object 
# seurat_object[[paste0(region, '_clusters')]] <- recode_cluster_ids(seurat_object, 
#                                                                    region, 
#                                                                    'cluster_full')
# 
# # Set Idents 
# message("Any NAs after recode of cluster IDs ?",  seurat_object[[paste0(region, '_clusters')]] |> anyNA())
# #Idents(seurat_object) <- seurat_object[[paste0(region, '_clusters')]]

# Calculate aggregate and average expression - adult
ag_adult_mat <- calculate_aggregated_expression(seurat_object, dystonia_genes)
# av_adult_mat <- calculate_average_expression(seurat_object, dystonia_genes)
# 
# ag_fetal_mat <- calculate_aggregated_expression(seurat_fetal, dystonia_genes)
# av_fetal_mat <- calculate_average_expression(seurat_fetal, dystonia_genes)
# 
# colnames(ag_fetal_mat) <- str_replace(colnames(ag_fetal_mat), '-', '-fetal-')
# colnames(av_fetal_mat) <- str_replace(colnames(av_fetal_mat), '-', '-fetal-')

# Save - note for snakemake fetal ge is called fetal str
message('\nSaving objects ...\n')
saveRDS(ag_adult_mat, paste0(R_dir, 'seurat_aggr_adult_', region, '.rds'))
#saveRDS(av_adult_mat, paste0(R_dir, 'seurat_aver_adult_', region, '.rds'))
#saveRDS(ag_fetal_mat, paste0(R_dir, 'seurat_aggr_fetal_', region, '.rds'))
#saveRDS(av_fetal_mat, paste0(R_dir, 'seurat_aver_fetal_', region, '.rds'))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
