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

# Switch to whole dataset
message('\nChanging to RNA object ...\n')
DefaultAssay(seurat_object) <- 'RNA'
seurat_object <- JoinLayers(seurat_object) # Do this for sketch and RNA independently

# Recode cluster IDs - sketch object 
seurat_object[[paste0(region, '_clusters')]] <- recode_cluster_ids(seurat_object, 
                                                                   region, 
                                                                   'cluster_full')

# Set Idents 
message("Any NAs after recode of cluster IDs?",  seurat_object[[paste0(region, '_clusters')]] |> anyNA())
Idents(seurat_object) <- seurat_object[[paste0(region, '_clusters')]]

# Calculate aggregate and average expression
ag_mat <- calculate_aggregated_expression(seurat_obj, dystonia_genes)
av_mat <- calculate_average_expression(seurat_obj, dystonia_genes)

# Save 
saveRDS(ag_mat, paste0(R_dir, 'seurat_aggr_exp_', region, '.rds'))
saveRDS(av_mat, paste0(R_dir, 'seurat_aver_exp_', region, '.rds'))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------