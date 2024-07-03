#--------------------------------------------------------------------------------------
#
#    Dystonia - hdWGCNA
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# Vingette: https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html

# Useful issues:

# https://github.com/smorabit/hdWGCNA/issues/275


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

### Need to configure these for hawk ####

# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# Markdown
library(kableExtra)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 4)

# load the Zhou et al snRNA-seq dataset
#seurat_obj <- readRDS('Zhou_2020.rds')

### Need to configure these for hawk ####

## Load Data --------------------------------------------------------------------------
# Seurat objects  ----
if (stringr::str_detect(region, 'fetal'))
  seurat_obj <- readRDS(paste0(fetal_dir, 'seurat_', region,'.rds'))
if (stringr::str_detect(region, 'adult'))
  seurat_obj <- readRDS(paste0(R_dir, 'seurat.', region, '.rds'))
  
# Set up object for WGCNA to get metacells
gene_select <- 'fraction'
seurat_obj <- create_wgcna_metacells(seurat_obj, gene_select, paste0(region, '_wgcna'))

# Run basic stats - move to markdown?
clusters_rm <- catch_clusters_rm_warning()
meta_obj <- GetMetacellObject(seurat_obj)
cell_props <- get_wgcna_cell_props(seurat_obj, meta_obj)

meta_cell_types <- colnames(meta_obj) %>%
  as_tibble() %>%
  separate(value, c('value', NA), sep = '#') %>%
  distinct() %>%
  pull()

run_wgcna_orig(seurat_obj, meta_cell_types, region, wgcna_dir)

### ------ Testing -------

# Run again to set up metacell type specific objects in seurat_obj@misc
seurat_obj <- create_wgcna_metacells(seurat_obj, 
                                     gene_select, 
                                     paste0(meta_cell_types[i], '_wgcna'),
                                     paste0(region, '_wgcna'),
                                     meta_cell_types)

# Run WGCNA - Throwing errors
run_wgcna(seurat_obj, meta_cell_types, region, wgcna_dir)

### ------

# Run overlaps between dystonia genes and WGCNA modules - atm 50 hub genes
overlap_genes <- run_dyst_gene_overlap(seurat_obj, meta_cell_types, region, 
                                       paste0(region, '_wgcna'), wgcna_dir)

# Get 
wgcna_stats_tbl <- get_wgcna_stats(seurat_obj, meta_cell_types, region, 
                                   paste0(region, '_wgcna'), wgcna_dir)

  
## Create markdown doc  ---------------------------------------------------------------
rmarkdown::render(markdown_wgcna_doc, output_file = markdown_wgcna_html, output_dir = wgcna_dir)


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


