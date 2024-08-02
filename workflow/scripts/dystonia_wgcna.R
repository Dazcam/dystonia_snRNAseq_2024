#--------------------------------------------------------------------------------------
#
#    Dystonia - hdWGCNA
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# Vingette: https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html

# Useful issues:

# https://github.com/smorabit/hdWGCNA/issues/275
# https://github.com/smorabit/hdWGCNA/issues/27


##  Load Packages, functions and variables  -------------------------------------------
message('Setting environment variables ...')
if (Sys.info()[["nodename"]] == "Darrens-iMac-2.local") {
  
  library(yaml)
  root_dir <- '~/Desktop/dystonia_snRNAseq_2024/'
  yaml_file <- yaml.load_file(paste0(root_dir, 'config/config.yaml'))
  region <- yaml.load(yaml_file$wgcna_region) # Note this is different from other scripts
  
  source(paste0(root_dir, 'workflow/scripts/dystonia_functions.R'))
  source(paste0(root_dir, 'workflow/scripts/dystonia_Renvs.R'))
  source(paste0(root_dir, 'workflow/scripts/dystonia_gene_lists.R'))
  
} else {
  
  source('scripts/dystonia_functions.R')
  source('scripts/dystonia_Renvs.R')
  source('scripts/dystonia_gene_lists.R')
  
}

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 20)

## Load Data --------------------------------------------------------------------------
# Seurat objects  ----
if (stringr::str_detect(region, 'fetal')) {
  seurat_obj <- readRDS(paste0(fetal_dir, 'seurat_', region,'.rds'))
} else { 
  seurat_obj <- readRDS(paste0(R_dir, '03seurat_', region, '.rds'))}

# Note hdWGCNA cannot handle BPCell or sketch objects
message("Converting BPCells counts and data matrix to in memory format for hdWGCNA ...")
DefaultAssay(seurat_obj) <- 'RNA'
seurat_obj <- JoinLayers(seurat_obj)
seurat_obj[["RNA"]]$counts <- as(object = seurat_obj[["RNA"]]$counts, Class = "dgCMatrix")
seurat_obj[["RNA"]]$data <- as(object = seurat_obj[["RNA"]]$data, Class = "dgCMatrix")

message("Recode cluster IDs ... ")
seurat_obj$cellIDs <- recode_cluster_ids(seurat_obj, region, 'cluster_full')
unique(seurat_obj$cellIDs)

if (subset_seurat_hdWGCNA == TRUE) {
  message("Subsetting seurat object for testing ...")
  seurat_obj <- subset(seurat_obj, downsample = 1000)}

## Aggregate cells? Set in dystonia_Renvs.R
#if (aggregate_cells == TRUE) 
#  seurat_obj <- recode_wgcna_clusters(seurat_obj, region)

counts_by_sample <- seurat_obj@meta.data %>%
  as_tibble() %>%
  select(cellIDs, sample_id) %>%
  group_by(sample_id) %>%
  dplyr::count(cellIDs)

# Set up object for WGCNA to get metacells
seurat_obj <- create_wgcna_metacells(seurat_obj, gene_select, paste0(region, '_wgcna'))



message("\nGet meta cell types ...\n")
meta_obj <- GetMetacellObject(seurat_obj)
meta_cell_types <- colnames(meta_obj) |>
  as_tibble() |>
  separate(value, c('value', NA), sep = '#') |>
  distinct() |>
  filter(str_detect(value, 'InN|ExN')) |>
  pull()

message("\nRunning the following meta cell types:\n")
message(paste0(capture.output(meta_cell_types), collapse = "\n"), '\n')

run_wgcna_orig(seurat_obj = seurat_obj, 
               cell_types = meta_cell_types, 
               cluster_column = 'cellIDs', 
               meta_column = 'sample_id',
               region = region, 
               outdir = wgcna_dir,
               wgcna_name = paste0(region, '_wgcna'))

### ------ Testing -------
# Aggregate all cell type specific WGCNA objects in single Seurat object
if (aggregate_misc == TRUE) {
  # Run again to set up metacell type specific objects in seurat_obj@misc
  seurat_obj <- create_wgcna_metacells(seurat_obj, 
                                       gene_select, 
                                       paste0(meta_cell_types[i], '_wgcna'),
                                       paste0(region, '_wgcna'),
                                       meta_cell_types)
  
  # Run WGCNA - Throwing errors
  run_wgcna(seurat_obj, meta_cell_types, region, wgcna_dir)}


# Issue from here is these function don't return the modified seurat object
# And currently there is a single object for each cell type per region

### ------

message("\nWriting metacell list to file ...\n")
write_tsv(as_tibble(meta_cell_types), paste0(wgcna_dir, region, '_metacells.tsv'), col_names = FALSE)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

