#--------------------------------------------------------------------------------------
#
#    Dystonia - Annotate clusters
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# Needs to be checked interactively and coordinated with gene and colour lists in 
# dystonia_gene_lists.R and recode_cluster_ids() 

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
seurat_object <- readRDS(paste0(R_dir, '02seurat_', region, '.rds'))

##  Plots for sketch object  ----------
# Recode cluster IDs - sketch object 
seurat_object[[paste0(region, '_clusters')]] <- recode_cluster_ids(seurat_object, 
                                                                   region, 
                                                                   'harmony_clusters_0.1')

# Plot paired vln
vln_plots_sketch <- plot_paired_vlns(seurat_object, 
                                     paste0(region, '_clusters'), 
                                     general_genes,
                                     get(paste0(region, '_genes')), 
                                     get(paste0(region, '_vln_cols_recode')))

# Plot paired umap and vln
umap_vln_plots_sketch <- plot_paired_umap_vln(seurat_object, 
                                              paste0(region, '_clusters'), 
                                              get(paste0(region, '_final_genes')),
                                              get(paste0(region, '_umap_cols_recode')), 
                                              get(paste0(region, '_vln_cols_recode')))

# Plot dystonia genes
dystonia_plot_sketch <- create_stacked_vln_plot(seurat_object, 
                                                paste0(region, '_clusters'), 
                                                dystonia_genes,
                                                toupper(region), 
                                                get(paste0(region, '_vln_cols_recode')))

##  Plots for sketch object  ----------
# Switch to whole dataset
message('\nChanging to RNA object ...')
DefaultAssay(seurat_object) <- 'RNA'

# Recode cluster IDs - sketch object 
seurat_object[[paste0(region, '_clusters')]] <- recode_cluster_ids(seurat_object, 
                                                                   region, 
                                                                   'harmony_clusters_0.1')

# Plot paired vln
vln_plots_rna <- plot_paired_vlns(seurat_object, 
                                  paste0(region, '_clusters'), 
                                  general_genes,
                                  get(paste0(region, '_genes')), 
                                  get(paste0(region, '_vln_cols_recode')))

# Plot paired umap and vln
umap_vln_plots_rna <- plot_paired_umap_vln(seurat_object, 
                                           paste0(region, '_clusters'), 
                                           get(paste0(region, '_final_genes')),
                                           get(paste0(region, '_umap_cols_recode')), 
                                           get(paste0(region, '_vln_cols_recode')))

# Plot dystonia genes
dystonia_plot_rna <- create_stacked_vln_plot(seurat_object, 
                                             paste0(region, '_clusters'), 
                                             dystonia_genes,
                                             toupper(region),
                                             get(paste0(region, '_vln_cols_recode')))


## Find differential expressed marker genes in clusters  ------------------------------
message('\nCalculating diff exp markers ...')
Idents(seurat_object) <- paste0(region, '_clusters')
marker_genes <- FindAllMarkers(seurat_object)
readr::write_tsv(marker_genes, paste0(R_dir, region, '_marker_genes.tsv'))

## Back to sketch
message('\nChanging to sketch object ...')
DefaultAssay(seurat_object) <- "sketch"

## Create markdown doc  ---------------------------------------------------------------
rmarkdown::render(markdown_ann_doc, output_file = markdown_ann_html, output_dir = R_dir)

saveRDS(seurat_object, paste0(R_dir, '03seurat_', region, '.rds'))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------


# Code for working out cell types genes / cols
# vln_pair <- plot_paired_vlns(seurat_object, 'harmony_clusters_0.1', general_genes,
#                              get(paste0(region, '_genes')))
# stiletti <- DimPlot_scCustom(seurat_object, group.by = 'cell_type',
#                              alpha = 0.1, pt.size = 0.1,
#                              repel = T)
# clust <- DimPlot_scCustom(seurat_object, group.by = 'harmony_clusters_0.1',
#                           alpha = 0.1, pt.size = 0.1, label = T)
# 
# vln_pair | clust / stiletti


