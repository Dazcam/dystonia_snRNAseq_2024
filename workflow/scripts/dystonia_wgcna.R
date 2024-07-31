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
message("Converting matrix from BPCells to in memory ...")
DefaultAssay(seurat_obj) <- 'RNA'
seurat_obj <- JoinLayers(seurat_obj)
seurat_obj[["RNA"]]$counts <- as(object = seurat_obj[["RNA"]]$counts, Class = "dgCMatrix")
seurat_obj[["RNA"]]$data <- as(object = seurat_obj[["RNA"]]$data, Class = "dgCMatrix")
seurat_obj[["RNA"]]$counts[1:10, 1:10]

message("Recode cluster IDs ... ")
seurat_obj$cellIDs <- recode_cluster_ids(seurat_obj, region, 'cluster_full')
unique(seurat_obj$cellIDs)

message("Subsetting seurat object for testing ...")
seurat_obj <- subset(seurat_obj, downsample = 1000)

# Recode clusters if adult data and set to 'cellIDs'
#if (!stringr::str_detect(region, 'fetal')) {
#  seurat_obj$cellIDs <- recode_cluster_ids(seurat_obj, region, 'harmony_clusters_0.1')
#  seurat_obj$Sample <- seurat_obj$sample_id
#  DefaultAssay(seurat_obj) <- 'sketch'
#  seurat_obj <- JoinLayers(seurat_obj) # Seurat 5 objects layers must be joined
#  }

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

head(seurat_obj@misc)

# Run basic stats - move to markdown?
clusters_rm <- catch_clusters_rm_warning()
meta_obj <- GetMetacellObject(seurat_obj)
cell_props <- get_wgcna_cell_props(seurat_obj, meta_obj)

message("\nGet meta cell types ...\n")
meta_cell_types <- colnames(meta_obj) %>%
  as_tibble() %>%
  separate(value, c('value', NA), sep = '#') %>%
  distinct() %>%
  pull()


meta_cell_types

message("\nSet dat exp ...\n")
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "Str-adult-InN-1", # the name of the group of interest in the group.by column
  group.by='cellIDs', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

# Select soft power threshold
message("\nGet soft power Thresh ...\n")
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed', # you can also use "unsigned" or "signed hybrid"
  wgcna_name = paste0(region, '_wgcna'))
    
message("\nSet power val ...\n")
power_val <- GetPowerTable(seurat_obj, paste0(region, '_wgcna')) %>%
  select(Power, SFT.R.sq) %>%
  filter(SFT.R.sq > 0.8) %>%
  pull(Power) %>%
  dplyr::first()
    
  message("\nSoft Power threshold set to: ", power_val, '\n')

message("\nConstruct co-expression network ...\n")
 seurat_obj <- ConstructNetwork(
      seurat_obj,
      tom_name = 'Str-adult-InN-1', # name of the topoligical overlap matrix written to disk
      soft_power = power_val,
      overwrite_tom = T,
      wgcna_name = paste0(region, '_wgcna'),
      tom_outdir = '../results/05wgcna/'
    )

# Required to avoid harmony error https://github.com/smorabit/hdWGCNA/issues/17
message("\nScaling Data ...\n")
seurat_obj <- ScaleData(seurat_obj)

    message("\nCompute Eigengenes ...\n")
    # Compute Eigengenes and Connectivity
    # Compute all MEs in the full single-cell dataset
    seurat_obj <- ModuleEigengenes(
      seurat_obj,
      modules = GetModules(seurat_obj, wgcna_name = paste0(region, '_wgcna')),
      group.by.vars = "sample_id",
      wgcna_name = paste0(region, '_wgcna')
    )

    message("\nModule Connect ...\n")
    # Compute module connectivity
    # compute eigengene-based connectivity (kME):
    seurat_obj <- ModuleConnectivity(
      seurat_obj,
      group.by = 'cellIDs', 
      group_name = 'Str-adult-InN-1',
      wgcna_name = paste0(region, '_wgcna')
    )
    
    message("\nRename modules ...\n")
    # rename the modules
    seurat_obj <- ResetModuleNames(
      seurat_obj,
      new_name = paste0("Str-adult-InN-1-M"),
      wgcna_name = paste0(region, '_wgcna')
    )

  message("\nSaving RDS file ...\n")
  saveRDS(seurat_obj, file = '../results/03wgcna/Str-adult-InN-1.rds')

#run_wgcna_orig(seurat_obj, meta_cell_types, region, wgcna_dir)

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

# Run overlaps between dystonia genes and WGCNA modules - atm 50 hub genes
message("Calculating overlap genes ...")
overlap_genes <- run_dyst_gene_overlap(seurat_obj, 'Str-adult-InN-1', region, 
                                       paste0(region, '_wgcna'), wgcna_dir)
overlap_genes

message("Getting stats tbl ...")
wgcna_stats_tbl <- get_wgcna_stats(seurat_obj, 'Str-adult-InN-1', region, 
                                   paste0(region, '_wgcna'), wgcna_dir)

  
## Create markdown doc  ---------------------------------------------------------------
# Modify html name for testing - this will only work locally 
#markdown_wgcna_html <- paste0(str_split(markdown_wgcna_html, '\\.')[[1]][1], '_',
#                              str_extract(gene_select, "^.{3}"), 
#                              if (aggregate_cells) '_agg.' else '.',
#                              str_split(markdown_wgcna_html, '\\.')[[1]][2])

rmarkdown::render(markdown_wgcna_doc, output_file = markdown_wgcna_html, output_dir = wgcna_dir)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

