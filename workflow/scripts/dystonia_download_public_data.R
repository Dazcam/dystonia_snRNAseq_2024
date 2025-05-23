#--------------------------------------------------------------------------------------
#
#    Download public data - Stiletti
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
# Download dissection data - run once - requires internet access
get_dissection_data(get(paste0(region, '_anns')), anns_table, R_dir, file_format = '.rds')

# Import Seurat object as BPCells object
seurat_object <- create_BPCell_seurat_object(get(paste0(region, '_anns')), '.rds', R_dir)

# Save Seurat object
saveRDS(seurat_object, paste0(R_dir, '01seurat_', region, '.rds'))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

                    
