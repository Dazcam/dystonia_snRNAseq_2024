#--------------------------------------------------------------------------------------
#
#    R environments 
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# Load R environment - options are local or cluster

##  Load Packages  --------------------------------------------------------------------
message('Loading packages ...')
library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(SeuratWrappers)
library(Azimuth) 
library(dplyr)
library(scCustomize)
library(readxl)
library(cowplot)
library(scuttle)
library(scater)

## Set variables  ---------------------------------------------------------------------
if (exists("snakemake")) { 
  root_dir <- snakemake@params[['root_dir']]
  region <- snakemake@params[['region']]
  future::plan("multicore", workers = snakemake@threads) 
  plan()
  log_smk() 
  # Also try this: plan("cluster")
  #future::plan("cluster", workers = 19)
  #future::plan()
}

data_dir <- paste0(root_dir, 'resources/')
script_dir <- paste0(root_dir, 'workflow/scripts/')
results_dir <- paste0(root_dir, 'results/')
stiletti_dir <- paste0(data_dir, 'public_data/stiletti_2023/')
R_dir <- paste0(results_dir, '01R_objects/')
markdown_doc <- paste0(script_dir, 'dystonia_qc.Rmd')
markdown_html <- paste0('dystonia_qc_', toupper(region), '.html')
regions <- c('fcx', 'str', 'cer')
fcx_anns <- c('A13', 'A14', 'A25', 'A32', 'A44-A45', 'A46', 'FI', 'M1C')
str_anns <- c('CaB', 'Pu')
cer_anns <- c('CBL', 'CBV', 'CbDN')
all_anns <- c(fcx_anns, str_anns, cer_anns)
anns_table <- read_excel(paste0(data_dir, 'sheets/Stiletti_downloads_table.xlsx'))
dystonia_genes <- read_excel(paste0(data_dir, 'sheets/Dystonia_Genes_Clinical_5.0.xlsx'), range = 'D1:D26') %>%
  pull(GeneName) 
sample_split <- 'sample_id' # Meta_id col to split seurat object by for qc
resolution_set <- seq(0.1, 0.8, 0.1) # Set of res params to test
pc_thresh <- ifelse(region == 'fcx', 50, 30) # Set PC thresh

options(future.globals.maxSize = 3e+09, future.seed = T) # set this option when analyzing large datasets
options(digits = 1) # Set default decimal points
options(scipen = 999)
options(ggrepel.max.overlaps = Inf) # For DimPlots




#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
