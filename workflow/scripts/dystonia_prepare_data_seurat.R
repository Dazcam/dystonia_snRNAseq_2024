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
library(scuttle)
library(scater)
options(future.globals.maxSize = 3e+09) # set this option when analyzing large datasets
options(digits = 1) # Set default decimal points
options(scipen = 999)

## Set variables  ---------------------------------------------------------------------
data_dir <- '~/Desktop/dystonia_snRNAseq_2024/resources/'
script_dir <- '~/Desktop/dystonia_snRNAseq_2024/workflow/scripts/'
results_dir <- '~/Desktop/dystonia_snRNAseq_2024/results/'
stiletti_dir <- paste0(data_dir, 'public_data/stiletti_2023/')
R_dir <- paste0(results_dir, '01R/')
markdown_doc <- paste0(script_dir, 'dystonia_qc.Rmd')
markdown_html <- 'dystonia_qc.html'
regions <- c('fcx', 'str', 'glp', 'cer')
fcx_anns <- c('A13', 'A14', 'A25', 'A32', 'A44-A45', 'A46', 'FI', 'M1C')
str_anns <- c('CaB', 'Pu')
cer_anns <- c('CBL', 'CBV', 'CbDN')
glp_anns <- c('GPi', 'GPe')
all_anns <- c(fcx_anns, str_anns, cer_anns, glp_anns)
anns_table <- read_excel(paste0(data_dir, 'sheets/Stiletti_downloads_table.xlsx'))
dystonia_genes <- read_excel(paste0(data_dir, 'sheets/Dystonia_Genes_Clinical_5.0.xlsx'), range = 'D1:D26') %>%
  pull(GeneName) 
  

source('~/Desktop/dystonia_snRNAseq_2024/workflow/scripts/dystonia_functions.R')
source('~/Desktop/dystonia_snRNAseq_2024/workflow/scripts/dystonia_gene_lists.R')
  
## Load Data --------------------------------------------------------------------------
# Download dissection data
get_dissection_data(glp_anns, anns_table, R_dir, file_format = '.rds')
get_dissection_data(fcx_anns, anns_table, R_dir, file_format = '.rds')
get_dissection_data(str_anns, anns_table, R_dir, file_format = '.rds')
get_dissection_data(cer_anns, anns_table, R_dir, file_format = '.rds')

# Get the original cluster annotations
stiletti_cluster_anns <- get_cluster_anns(stiletti_dir)

# Write the counts layer to a directory
# write_matrix_dir(mat = gpi_obj[["RNA"]]$data, dir = '~/Desktop/test/gpi')
# counts.mat <- open_matrix_dir(dir = '~/Desktop/test/gpi')
seurat_str <- create_BPCell_seurat_object(str_anns, '.rds', R_dir) 
seurat_glp <- create_BPCell_seurat_object(glp_anns, '.rds', R_dir) 
seurat_cer <- create_BPCell_seurat_object(cer_anns, '.rds', R_dir)
seurat_fcx <- create_BPCell_seurat_object(fcx_anns, '.rds', R_dir) 

# Basic qc plots
basic_qc_plot_str <- create_basic_qc_plots(seurat_str) 
basic_qc_plot_glp <- create_basic_qc_plots(seurat_glp) 
basic_qc_plot_cer <- create_basic_qc_plots(seurat_cer) 
basic_qc_plot_fcx <- create_basic_qc_plots(seurat_fcx) 

# Join layers to apply filters
seurat_str_join <- JoinLayers(seurat_str)


# Not finished
# Necessary?
# remove_undetected_genes <- function(
#     
#   seurat_obj = NULL
#     
#   ) {
#   
#   genes_per_cell <- colSums(seurat_obj[["RNA"]]$counts)
#   plot(density(genes_per_cell), main="", xlab="Genes per cell")
#   detected_genes <- rowSums(seurat_str_join[["RNA"]]$counts) > 0
#   seurat_obj <- seurat_obj[detected_genes, ]
#   
#   #return(seurat_obj)
#   
# }

###. NEW FILTERING SECTION. -----------------------------------------------------------

# Biol Psych thresholds set for fetal samples
# Cells exp < 1000 genes > 5000 genes
# > 5% mito
# > 10% ribo
# Genes from mito genome excluded
# Genes expressed in fewer than 3 cells excluded
# Doublets

#Note: Big difference in UMI / Genes captured may be worth downsampling as control / check

# Convert to single cell experiment
sce_str <- SingleCellExperiment(list(counts = as(seurat_str_join[["RNA"]]$counts, "dgCMatrix")),
                                colData = seurat_str_join@meta.data)

sce_str <- get_cell_outliers(sce_str, 3, 'both', 5, 5)

# Remove genes from the seurat object - not SCE!!
sce.cer <- sce.cer[!grepl("MALAT1", rownames(sce.cer)), ]
sce.cer <- sce.cer[!grepl("^MT-", rownames(sce.cer)), ]

# Discard cell outliers
discard <- !all_outliers
seurat_str_join$discard <- discard
subset(x = seurat_str_join, subset = discard == TRUE)

create_outlier_plots(sce_str)

            
# These are scCustomize functions - decided to go Bioconductor route            
# sum(duplicated(rownames(sce)))
# seurat_str_join <-  Add_Mito_Ribo_Seurat(seurat_object = seurat_str_join, species = "Human", overwrite = T)
# seurat_str_join <-  Add_Cell_Complexity_Seurat(seurat_object = seurat_str_join, species = "Human", overwrite = T)       
# 
# MAD_stats <- Median_Stats(seurat_object = seurat_str_join, group_by_var = "orig.ident")[1:3] %>%
#   mutate(mad_x_3_nCount = 3 * mad(seurat_str_join$nCount_RNA))
#   mutate()



get_meta_col_counts(seurat_str, 'dissection')

### END FILTERING SECTION. ------------------------------------------------------------


# Join and splitting object - not necessary here as objects are already split
# seurat_str[["RNA"]] <- split(seurat_str[["RNA"]], f = seurat_str$dataset)
# seurat_join <- JoinLayers(seurat_str[["RNA"]])
seurat_sk_str <- create_sketch_object(seurat_str, 30)
seurat_sk_glp <- create_sketch_object(seurat_glp, 30)
seurat_sk_cer <- create_sketch_object(seurat_cer, 30)
seurat_sk_fcx <- create_sketch_object(seurat_fcx, 30)

saveRDS(object = seurat_sk_str, file = paste0(R_dir, "seurat_sk_str.Rds"))
saveRDS(object = seurat_sk_glp, file = paste0(R_dir, "seurat_sk_glp.Rds"))
saveRDS(object = seurat_sk_cer, file = paste0(R_dir, "seurat_sk_cer.Rds"))
saveRDS(object = seurat_sk_fcx, file = paste0(R_dir, "seurat_sk_fcx.Rds"))

# Plot QCs
cluster_qc_plot_str <- create_cluster_qc_plot(seurat_sk_str, 30)
cluster_qc_plot_glp <- create_cluster_qc_plot(seurat_sk_glp, 30)
cluster_qc_plot_cer <- create_cluster_qc_plot(seurat_sk_cer, 30)
cluster_qc_plot_fcx <- create_cluster_qc_plot(seurat_sk_fcx, 30)

# Cluster counts - Need to change to full dataset first data first
seurat_sk_fcx$seurat_clusters

# Choose PC threshold
pca_plot_str <- DimHeatmap(seurat_sk_str, dims = 20:30, cells = 500, balanced = TRUE)
pca_plot_glp <- DimHeatmap(seurat_sk_glp, dims = 20:30, cells = 500, balanced = TRUE)
pca_plot_cer <- DimHeatmap(seurat_sk_cer, dims = 20:30, cells = 500, balanced = TRUE)
pca_plot_fcx <- DimHeatmap(seurat_sk_fcx, dims = 20:30, cells = 500, balanced = TRUE)

## Integration
seurat_sk_str <- run_integration(seurat_sk_str)
seurat_sk_glp <- run_integration(seurat_sk_glp)
seurat_sk_cer <- run_integration(seurat_sk_cer)
seurat_sk_fcx <- run_integration(seurat_sk_fcx)

int_plot_str <- create_integration_plot(seurat_sk_str) 
int_plot_glp <- create_integration_plot(seurat_sk_glp) 
int_plot_cer <- create_integration_plot(seurat_sk_cer) 
int_plot_fcx <- create_integration_plot(seurat_sk_fcx) 

# Join layers - Need to do this before stacked vln plotting
seurat_sk_str <- JoinLayers(seurat_sk_str)
seurat_sk_glp <- JoinLayers(seurat_sk_glp)
seurat_sk_cer <- JoinLayers(seurat_sk_cer)
seurat_sk_fcx <- JoinLayers(seurat_sk_fcx)

# Gene check
vln_plot_general_str <- create_stacked_vln_plot(seurat_sk_str, general_genes, 'Striatum')
vln_plot_general_glp <- create_stacked_vln_plot(seurat_sk_glp, general_genes, 'Globus Pallidus')
vln_plot_general_cer <- create_stacked_vln_plot(seurat_sk_cer, general_genes, 'Cerebellum')
vln_plot_general_fcx <- create_stacked_vln_plot(seurat_sk_fcx, general_genes, 'Frontal Cortex')

# Specific gene check
vln_plot_specific_cer <- create_stacked_vln_plot(seurat_sk_cer, cer_genes, 'Cerebellum')
vln_plot_specific_fcx <- create_stacked_vln_plot(seurat_sk_fcx, fcx_genes, 'Frontal Cortex')

# Dystonia gene check
vln_plot_dyst_str <- create_stacked_vln_plot(seurat_sk_str, dystonia_genes, 'Striatum')
vln_plot_dyst_glp <- create_stacked_vln_plot(seurat_sk_glp, dystonia_genes, 'Globus Pallidus')
vln_plot_dyst_cer <- create_stacked_vln_plot(seurat_sk_cer, dystonia_genes, 'Cerebellum')
vln_plot_dyst_fcx <- create_stacked_vln_plot(seurat_sk_fcx, dystonia_genes, 'Frontal Cortex')



create_stacked_vln_plot(seurat_sk_cer, cer_genes, 'Cerebellum')

# Run Nick's average expression analysis with the 25 dyst genes
# Extract av. exp. for each cluster across each region fetal / adult for the 25 dyst genes merge 
# everything into one matrix and then cluster to ID patterns (use PCA in the first instance)
# Note these analyses may have been run on sketch data
# Load fetal data
for (region in c('cer', 'hip', 'pfc', 'wge', 'tha')) {
  
  seurat_obj <- readRDS(paste0('~/Desktop/fetal_brain_snRNAseq_110122/resources/R_objects/seurat.', region, '.final.rds'))
  assign(paste0('seurat_fetal_', region), seurat_obj, envir = .GlobalEnv)
  
}

av_exp_mat_cer_ftl <- calculate_average_expression(seurat_fetal_cer, 'cer_fetal', dystonia_genes)
av_exp_mat_fcx_ftl <- calculate_average_expression(seurat_fetal_pfc, 'fcx_fetal', dystonia_genes)
av_exp_mat_hip_ftl <- calculate_average_expression(seurat_fetal_hip, 'hip_fetal', dystonia_genes)
av_exp_mat_tha_ftl <- calculate_average_expression(seurat_fetal_tha, 'tha_fetal', dystonia_genes)
av_exp_mat_gem_ftl <- calculate_average_expression(seurat_fetal_wge, 'gem_fetal', dystonia_genes)

av_exp_mat_str <- calculate_average_expression(seurat_sk_str, 'str_adult', dystonia_genes)
av_exp_mat_glp <- calculate_average_expression(seurat_sk_glp, 'glp_adult', dystonia_genes)
av_exp_mat_cer <- calculate_average_expression(seurat_sk_cer, 'cer_adult', dystonia_genes)
av_exp_mat_fcx <- calculate_average_expression(seurat_sk_fcx, 'fcx_adult', dystonia_genes)

# AggregateExpression(seurat_sk_cer, features = dystonia_genes)
av_exp_mat <- cbind(av_exp_mat_str, av_exp_mat_glp, av_exp_mat_cer, av_exp_mat_fcx)
av_exp_mat_scaled <- t(apply(av_exp_mat, 1, scale)) # centre and scale each col (Z-Score) then t
colnames(av_exp_mat_scaled) <- colnames(av_exp_mat) 

av_exp_pca <- prcomp(t(av_exp_mat_scaled))
ComplexHeatmap::Heatmap(t(av_exp_mat_scaled), column_labels = rownames(av_exp_mat_scaled), cluster_columns = F, cluster_rows = T)


### test code 
gc_vst <- read.table("~/Downloads/counts_vst.txt", header = T, row.names = 1, sep = "\t")
vst_pca <- prcomp(t(gc_vst))

# Calculate variance explained
frac_var <- function(x) x^2 / sum(x^2)

av_exp_pca$sdev %>% 
  as_tibble() %>% 
  frac_var() %>% 
  mutate(Comp = colnames(av_exp_pca$x)) %>% 
  slice(1:9) %>% 
  ggplot(aes(x=Comp, y = value)) + 
  geom_bar(stat = "identity", fill = "#4DC5F9") +
  geom_hline(yintercept = 0.03, linetype = 2) +
  xlab("Principal Components") +
  scale_y_continuous(name = "Variance Explained", breaks = seq(0,0.8,0.1), labels = scales::percent_format(accuracy = 5L)) +
  theme_classic(base_size = 14)

# Prep PCA plot
genes.selected <- av_exp_pca$rotation[c(which.max(av_exp_pca$rotation[,"PC1"]), 
                                  which.min(av_exp_pca$rotation[,"PC1"]), 
                                  which.max(av_exp_pca$rotation[,"PC2"]), 
                                  which.min(av_exp_pca$rotation[,"PC2"])),
                                c("PC1", "PC2")]
genes.selected <- genes.selected %>%
  as.data.frame() %>%
  rownames_to_column(var = "genes")
genes.selected

pca_all <- av_exp_pca$x %>%
  as_tibble(rownames = 'cell_type') %>%
  mutate(region = str_split_i(cell_type, pattern = '_', 1)) %>%
  ggplot(aes(x= PC1, y= PC2, color = region)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype=2) +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

ggplot(genes_tbl, aes(x=PC1, y=PC2)) +
  geom_point() +
  geom_segment(aes(xend=PC1, yend=PC2), x=0, y=0, color="Grey") +
  geom_label(aes(x=PC1, y=PC2, label=genes), size=2, vjust="outward") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  panel.grid.minor = element_blank()) 


# Bioconductor method - https://www.bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html
pca_test <- PCAtools::pca(av_exp_mat_scaled, scale = F, metadata = data.frame(region = str_split_i(colnames(av_exp_mat_scaled), 'ult', 1), 
                                                                              row.names = colnames(av_exp_mat_scaled)))
scree_plot <- PCAtools::screeplot(pca_test, axisLabSize = 18, titleLabSize = 22)
pca_plot <- PCAtools::biplot(pca_test, showLoadings = FALSE, ntopLoadings = 25, colby = 'region',
                            hline = 0, vline = 0)
pca_load_plot <- PCAtools::biplot(pca_test, showLoadings = TRUE, ntopLoadings = 25, colby = 'region',
                            hline = 0, vline = 0)
pairs_plot <- PCAtools::pairsplot(pca_test, colby = 'region')
plotloadings(pca_test,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1, PC2, PC3, PC4, PC5', 
             caption = 'Top 1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)

plot_grid(pca_all, pca_bio, align = 'hvlr')

## ------ Still to do
# mad_stats <- MAD_Stats(seurat_object = seurat_sk_fcx, group_by_var = "orig.ident", mad_num = 3)
# Number of undetected genes
detected_genes <- rowSums(counts(sce)) > 0
table(detected_genes)


## Create markdown doc  ---------------------------------------------------------------
rmarkdown::render(markdown_doc, output_file = markdown_html, output_dir = R_dir)


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


# Join the layers
seurat_fcx[["sketch"]] <- JoinLayers(seurat_fcx[["sketch"]])

# Add cluster annotations to seurat metadata
ann <- seurat_sketch@meta.data %>%
  as_tibble(rownames = 'cells') %>%
  mutate(cluster_id = as.double(cluster_id)) %>%
  left_join(stiletti_cluster_anns, by = join_by(cluster_id)) %>%
  pull(clusterAnn)
  

          


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
