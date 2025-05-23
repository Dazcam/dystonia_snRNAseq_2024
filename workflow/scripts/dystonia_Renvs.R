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
library(egg) 
#library(MatrixExtra) # This clashes with Seurat. Used here for hdWGCNA markdown, if needed, call it explicitly 
# Need to load within script just before markdown. See https://github.com/satijalab/seurat/issues/9221

library(tidyverse)
library(cowplot)
library(patchwork)
#library(WGCNA)
#library(hdWGCNA)

library(kableExtra)

#library(enrichR)
#library(GeneOverlap)
library(rmarkdown)

library(EWCE)
library(ewceData)

# Only need these package locally
if (Sys.info()[["nodename"]] == "Darrens-iMac-2.local") {
  library(pheatmap)
  library(grid)
  library(ComplexHeatmap)
  library(biomaRt)
}

## Set variables  ---------------------------------------------------------------------
if (exists("snakemake")) { 
  log_smk() 
  root_dir <- snakemake@params[['root_dir']]
  region <- snakemake@params[['region']]
  threads <- snakemake@threads

# Only works atm when snakemake@params[['wgcna']] is present in rule  
#  if (exists(snakemake@params[['wgcna']])) {
#    enableWGCNAThreads(snakemake@threads)
#    cat("WGCNA threads enabled. Set to: ", threads, '\n')
#  } else {
#    cat("Threads are set to: ", threads, '\n')
#    future::plan("multicore", workers = snakemake@threads)
#    plan()
#  }

}

data_dir <- paste0(root_dir, 'resources/')
script_dir <- paste0(root_dir, 'workflow/scripts/')
results_dir <- paste0(root_dir, 'results/')
stiletti_dir <- paste0(data_dir, 'public_data/stiletti_2023/')
fetal_dir <- paste0(data_dir, 'public_data/cameron_2023/')
R_dir <- paste0(results_dir, '01R_objects/')
bulk_dir <- paste0(results_dir, '02Bulk_data/')
table_dir <- paste0(results_dir, '05tables/')
fig_dir <- paste0(results_dir, '04figs/')
sheets_dir <- paste0(root_dir, 'resources/sheets/')
brain_span_dir <- paste0(data_dir, 'public_data/brain_span/genes_matrix_csv/')
wgcna_dir <- paste0(results_dir, '03wgcna/')
specificity_dir <- paste0(root_dir, 'resources/Source_Data_v241210/Specificity_per_cell_type/')
markdown_prep_doc <- paste0(script_dir, 'dystonia_qc.Rmd')
markdown_prep_html <- paste0('dystonia_qc_', region, '.html')
markdown_ann_doc <- paste0(script_dir, 'dystonia_ann.Rmd')
markdown_ann_html <- paste0('dystonia_ann_', region, '.html')
markdown_wgcna_doc <- paste0(script_dir, 'dystonia_wgcna.Rmd')
markdown_wgcna_html <- paste0('dystonia_wgcna_', region, '.html')
regions <- c('fcx', 'str', 'cer')
fcx_anns <- c('A13', 'A14', 'A25', 'A32', 'A44-A45', 'A46', 'FI', 'M1C')
str_anns <- c('CaB', 'Pu')
cer_anns <- c('CBL', 'CBV', 'CbDN')
all_anns <- c(fcx_anns, str_anns, cer_anns)
anns_table <- read_excel(paste0(data_dir, 'sheets/Stiletti_downloads_table.xlsx'))
sample_split <- 'sample_id' # Meta_id col to split seurat object by for qc
resolution_set <- seq(0.1, 0.8, 0.1) # Set of res params to test
pc_thresh <- ifelse(region == 'fcx', 50, 30) # Set PC thresh

options(future.globals.maxSize = 3e+09, future.seed = T) # set this option when analyzing large datasets
options(digits = 1) # Set default decimal points
options(scipen = 999) # Prevents wonky scientific notation
options(ggrepel.max.overlaps = Inf) # For DimPlots

dystonia_genes <- c("ACTB", "ADCY5", "ANO3", "AOPEP", "ATP1A3", "BCAP31", "CACNA1A", 
                    "COX20", "DDC", "DNAJC12", "EIF2AK2", "FITM2", "FOXG1", "GCH1", 
                    "GNAL", "GNAO1", "GNB1", "HPCA", "KCNA1", "KCNMA1", "KCTD17", 
                    "KMT2B", "MECR", "PNKD", "PRKRA", "RHOBTB2", "SCN8A", "SERAC1", 
                    "SGCE", "SLC2A1", "SLC6A3", "SPR", "SQSTM1", "TAF1", "TH", "THAP1", 
                    "TIMM8A", "TMEM151A", "TOR1A", "TSPOAP1", "TUBB4A", "VAC14", 
                    "VPS16", "YY1")
message(length(dystonia_genes), ' genes in dystonia hgnc_hg38 gene list.')

# Set num of cells per sample to create sketch object of ~60K
# Note FCX sample 10X418_6 only has 366 cells
sketch_num <- dplyr::case_when(
  region == "fcx" ~ 1100, # 56 samples (56 * 1100 = 61600 cells)
  region == "str" ~ 5000, # 12 samples (12 * 5000 = 60000 cells)
  region == "cer" ~ 2100) # 28 samples (28 * 2200 = 60200 cells)

### For final violin plots

# Set fetal region based on adult region
fetal_region <- dplyr::case_when(
  region == "fcx" ~ "fcx_fetal", 
  region == "str" ~ "ge_fetal", 
  region == "cer" ~ "cer_fetal")

# Region recode for cell IDs in seurat objects
region_recode <- dplyr::case_when(
  region == "fcx" ~ "FC", 
  region == "str" ~ "GE", 
  region == "cer" ~ "Cer")

adult_title <- dplyr::case_when(
  region == "fcx" ~ "Frontal Cortex", 
  region == "str" ~ "Striatum", 
  region == "cer" ~ "Cerebellum")

fetal_title <- dplyr::case_when(
  region == "fcx" ~ "Fetal Frontal Cortex", 
  region == "str" ~ "Fetal Ganglionic Eminences", 
  region == "cer" ~ "Fetal Cerebellum")


### For WGCNA
subset_seurat_hdWGCNA = FALSE
gene_select <- 'fraction'
set_k <- 25
aggregate_cells <- FALSE
aggregate_misc <- FALSE

### For Brain span
brain_levels <- c("PFC", "PMSC", "NPFC", "Str", "Cer", "Hip", "Tha", "Amy")
dev_levels <- c("EarlyFetal", "MidFetal", "LateFetal", "Infancy", "Childhood",
                "Adolescence", "Adulthood")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
