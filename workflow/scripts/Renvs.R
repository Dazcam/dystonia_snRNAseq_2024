#--------------------------------------------------------------------------------------
#
#    R environments 
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# Load R environment - options are local of cluster

##  Load Packages  --------------------------------------------------------------------
library(BPCells)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(SeuratWrappers)
library(Azimuth) 
library(tidyverse)
library(scCustomize)
library(readxl)
library(cowplot)
library(scuttle)
library(scater)

if (locale == 'local') { root_dir <- '~/Desktop/dystonia_snRNAseq_2024/'}
if (locale == 'remote') { 
  root_dir <- snakemake@params[['root_dir']]
  region <- snakemake@params[['region']]
}

## Set variables  ---------------------------------------------------------------------
data_dir <- paste0(root_dir, 'resources/')
script_dir <- paste0(root_dir, 'workflow/scripts/')
results_dir <- paste0(root_dir, 'results/')
stiletti_dir <- paste0(data_dir, 'public_data/stiletti_2023/')
R_dir <- paste0(results_dir, '01R/')
markdown_doc <- paste0(script_dir, 'dystonia_qc.Rmd')
markdown_html <- paste0('dystonia_qc_', toupper(region), '.html')

regions <- c('fcx', 'str', 'cer')
fcx_anns <- c('A13', 'A14', 'A25', 'A32', 'A44-A45', 'A46', 'FI', 'M1C')
str_anns <- c('CaB', 'Pu')
cer_anns <- c('CBL', 'CBV', 'CbDN')
#glp_anns <- c('GPi', 'GPe')
all_anns <- c(fcx_anns, str_anns, cer_anns)
anns_table <- read_excel(paste0(data_dir, 'sheets/Stiletti_downloads_table.xlsx'))
dystonia_genes <- read_excel(paste0(data_dir, 'sheets/Dystonia_Genes_Clinical_5.0.xlsx'), range = 'D1:D26') %>%
  pull(GeneName) 

options(future.globals.maxSize = 3e+09) # set this option when analyzing large datasets
options(digits = 1) # Set default decimal points
options(scipen = 999)

# Set region specific variables
# Remove orig.ident after running this remotely
if (region == 'fcx') {
  
  pc_thresh <- 50
  sample_split <- 'orig.ident'
  
} else if (region == 'cer') {
  
  pc_thresh <- 30
  sample_split <- 'orig.ident'
  
} else {
  
  pc_thresh <- 30
  sample_split <- 'sample_id'
  
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------