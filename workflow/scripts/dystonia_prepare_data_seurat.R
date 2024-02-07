#--------------------------------------------------------------------------------------
#
#    Stiletti analysis - Prepare data
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Prepare data
#  Stiletti data: https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443
#  Note that the Stiletti data has counts stored in gpi_obj[["RNA"]]$data

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
options(future.globals.maxSize = 3e+09) # set this option when analyzing large datasets


## Set variables  ---------------------------------------------------------------------
data_dir <- '~/Desktop/dystonia_snRNAseq_2024/resources/'
script_dir <- '~/Desktop/dystonia_snRNAseq_2024/workflow/scripts/'
results_dir <- '~/Desktop/dystonia_snRNAseq_2024/results/'
stiletti_dir <- paste0(data_dir, 'public_data/stiletti_2023/')
markdown_doc <- paste0(script_dir, 'dystonia_qc.Rmd')
markdown_dir <- paste0(results_dir, '01markdown/')
markdown_html <- 'dystonia_qc.html'
regions <- c('fcx', 'str', 'glp', 'cer')
fcx_anns <- c('A13', 'A14', 'A25', 'A32', 'A44-A45', 'A46', 'FI', 'M1C')
str_anns <- c('CaB', 'Pu')
cer_anns <- c('CBL', 'CBV', 'CbDN')
glp_anns <- c('GPi', 'GPe')
all_anns <- c(fcx_anns, str_anns, cer_anns, glp_anns)
anns_table <- read_excel(paste0(data_dir, 'sheets/Stiletti_downloads_table.xlsx'))


source('~/Desktop/dystonia_snRNAseq_2024/workflow/scripts/dystonia_functions.R')
  
## Load Data --------------------------------------------------------------------------
# Download dissection data
#get_dissection_data(glp_anns, anns_table, '~/Desktop/test/', file_format = '.rds')
#get_dissection_data(fcx_anns, anns_table, '~/Desktop/test/', file_format = '.rds')
#get_dissection_data(str_anns, anns_table, '~/Desktop/test/', file_format = '.rds')
#get_dissection_data(cer_anns, anns_table, '~/Desktop/test/', file_format = '.rds')

# Get the original cluster annotations
stiletti_cluster_anns <- get_cluster_anns(stiletti_dir)

# Write the counts layer to a directory
# write_matrix_dir(mat = gpi_obj[["RNA"]]$data, dir = '~/Desktop/test/gpi')
# counts.mat <- open_matrix_dir(dir = '~/Desktop/test/gpi')
seurat_str <- create_BPCell_seurat_object(str_anns, '.rds', '~/Desktop/test/') 
seurat_glp <- create_BPCell_seurat_object(glp_anns, '.rds', '~/Desktop/test/') 
seurat_cer <- create_BPCell_seurat_object(cer_anns, '.rds', '~/Desktop/test/')
seurat_fcx <- create_BPCell_seurat_object(fcx_anns, '.rds', '~/Desktop/test/') 

# Join and splitting object - not necessary here as objects are already split
# seurat_str[["RNA"]] <- split(seurat_str[["RNA"]], f = seurat_str$dataset)
# seurat_join <- JoinLayers(seurat_str[["RNA"]])
seurat_sk_str <- create_sketch_object(seurat_str, 30)
seurat_sk_glp <- create_sketch_object(seurat_glp, 30)
seurat_sk_cer <- create_sketch_object(seurat_cer, 30)
seurat_sk_fcx <- create_sketch_object(seurat_fcx, 30)

# Plot QCs
qc_plot_str <- create_qc_plot(seurat_sk_str, 30)
qc_plot_glp <- create_qc_plot(seurat_sk_glp, 30)
qc_plot_cer <- create_qc_plot(seurat_sk_cer, 30)
qc_plot_fcx <- create_qc_plot(seurat_sk_fcx, 30)

## Integration
seurat_sk_str <- run_integration_all(seurat_sk_str)
seurat_sk_glp <- run_integration_all(seurat_sk_glp)
seurat_sk_cer <- run_integration_all(seurat_sk_cer)
seurat_sk_fcx <- run_integration_all(seurat_sk_fcx)

int_plot_str <- create_integration_compare_plot(seurat_sk_str, c('harmony', 'cca', 'rpca', 'fastmnn'), 'dataset') 
int_plot_glp <- create_integration_compare_plot(seurat_sk_glp, c('harmony', 'cca', 'rpca', 'fastmnn'), 'dataset') 
int_plot_cer <- create_integration_compare_plot(seurat_sk_cer, c('harmony', 'cca', 'rpca', 'fastmnn'), 'dataset') 
int_plot_fcx <- create_integration_compare_plot(seurat_sk_fcx, c('harmony', 'cca', 'rpca', 'fastmnn'), 'dataset') 


## Create markdown doc  ---------------------------------------------------------------
rmarkdown::render(markdown_doc, output_file = markdown_html, output_dir = markdown_dir)


# Join the layers
seurat_sketch[["sketch"]] <- JoinLayers(seurat_sketch[["sketch"]])

# Add cluster annotations to seurat metadata
ann <- seurat_sketch@meta.data %>%
  as_tibble(rownames = 'cells') %>%
  mutate(cluster_id = as.double(cluster_id)) %>%
  left_join(stiletti_cluster_anns, by = join_by(cluster_id)) %>%
  pull(clusterAnn)
  

fc_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                 "MKI67", "C3", "ITM2A", "SST", "CALB2", 
                 "SCGN", "TLE3", "FEZF2", "CRYM", "LHX2")
ge_features <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                 "MKI67", "C3", "ITM2A", "LHX6", "SIX3", 
                 "PROX1", "TSHZ1", "DLX1", "SCGN")
hip_features <- c("NEUROD1", "GRIK4", "EOMES", "GLI3", "OLIG1", 
                  "MKI67", "C3", "ITM2A", "SLC17A6", "ADARB2",
                  "GAD2", "TNC", "PROX1", "RELN", "LHX6")
tha_features <- c("EOMES", "GLI3", "OLIG1", "MKI67", "C3", 
                  "ITM2A", "SLC1A6", "LHX9", "TNC", "GAD2", 
                  "PAX6", "SLC17A6")
cer_features <- c("GAD1", "EOMES", "GLI3", "OLIG1", "MKI67", 
                  "C3", "ITM2A", "CA8", "ITPR1", "RBFOX3", 
                  "RELN")

general_features <- c('VGLUT1', 'VGLUT2', 'SST', 'NPY', 'GAD2', 'C3',
                      'OLIG1', 'OLIG2')



seurat_merged <- SCTransform(seurat_str) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30)









for (i in 1:length(FILE_SET)) {
  path <- paste0(STILETTI_DIR, FILE_SET[i])
  data <- open_matrix_anndata_hdf5(path)
  write_matrix_dir(
    mat = data,
    dir = paste0(gsub(".h5ad", "", path), "_BP"),
    overwrite = TRUE
    
  )
  # Load in BP matrices
  mat <- open_matrix_dir(dir = paste0(gsub(".h5ad", "", path), "_BP"))
  mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
  
  # Get metadata
  metadata.list[[i]] <- LoadH5ADobs(path = path)
  data.list[[i]] <- mat
}

# Name layers
names(data.list) <- c("neurons", "non_neurons")

# Add Metadata
for (i in 1:length(metadata.list)) {
  metadata.list[[i]]$dataset <- names(data.list)[i]
}
metadata.list <- lapply(metadata.list, function(x) {
  x <- x[, c("dataset", "ROIGroup", "ROIGroupCoarse", "ROIGroupFine", "roi", "organism_ontology_term_id", 
             "disease_ontology_term_id", "self_reported_ethnicity_ontology_term_id", 
             "assay_ontology_term_id", "sex_ontology_term_id", "development_stage_ontology_term_id", 
             "donor_id", "suspension_type", "dissection", "fraction_mitochondrial", 
             "fraction_unspliced", "cell_cycle_score", "total_genes", "total_UMIs", 
             "sample_id", "supercluster_term", "cluster_id", "subcluster_id", 
             "cell_type_ontology_term_id", "tissue_ontology_term_id", "is_primary_data", 
             "cell_type", "assay", "disease", "organism", "sex", "tissue", 
             "self_reported_ethnicity", "development_stage")]
  return(x)
})

metadata <- Reduce(rbind, metadata.list)
merged.object <- CreateSeuratObject(counts = data.list, meta.data = metadata)
merged.object <- NormalizeData(merged.object, verbose = FALSE)
merged.object <- ScaleData(merged.object)

merged.object
