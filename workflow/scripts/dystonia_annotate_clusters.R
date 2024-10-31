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
message('Root dir set to:', root_dir)
message('\nLoading R object ...\n')
seurat_object <- readRDS(paste0(R_dir, '02seurat_', region, '.rds'))

##  Plots for sketch object  ----------
# Recode cluster IDs - sketch object 
message('\nChanging to Sketch object ...\n')
DefaultAssay(seurat_object) <- 'sketch'

message("Recode cluster IDs ... ")
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
                                              'umap.harmony',
                                              paste0(region, '_clusters'), 
                                              get(paste0(region, '_final_genes')),
                                              adult_title,
                                              get(paste0(region, '_umap_cols_recode')), 
                                              get(paste0(region, '_vln_cols_recode')))



# create_stacked_vln_plot(seurat_object, 
#                         paste0(region, '_clusters'), 
#                         dystonia_genes,
#                         adult_title, 
#                         get(paste0(region, '_vln_cols_recode')))

#DimPlot(seurat_object, group.by = 'cell_type', reduction = 'umap.harmony')


# Plot dystonia genes
dystonia_plot_sketch <- create_stacked_vln_plot(seurat_object, 
                                                paste0(region, '_clusters'), 
                                                dystonia_genes,
                                                toupper(region), 
                                                get(paste0(region, '_vln_cols_recode')))

##  Plots for whole object  ----------
# Switch to whole dataset
message('\nChanging to RNA object ...\n')
DefaultAssay(seurat_object) <- 'RNA'
seurat_object <- JoinLayers(seurat_object) # Do this for sketch and RNA independently

# Recode cluster IDs - sketch object 
message("Recode cluster IDs ... ")
seurat_object$cellIDs <- recode_cluster_ids(seurat_object, region, 'cluster_full')
Idents(seurat_object) <- seurat_object$cellIDs
message('Number of NAs in Idents: ', anyNA(Idents(seurat_object)))

# Plot paired vln
vln_plots_rna <- plot_paired_vlns(seurat_object, 
                                  'cellIDs', 
                                  general_genes,
                                  get(paste0(region, '_genes')), 
                                  get(paste0(region, '_vln_cols_recode')))

# Plot paired umap and vln
umap_vln_plots_rna <- plot_paired_umap_vln(seurat_object, 
                                           'umap.full',
                                           'cellIDs', 
                                           get(paste0(region, '_final_genes')),
                                           adult_title,
                                           get(paste0(region, '_umap_cols_recode')), 
                                           get(paste0(region, '_vln_cols_recode')))

umap_vln_plots_sketch <- plot_paired_umap_vln(seurat_object, 
                                              'umap.harmony',
                                              paste0(region, '_clusters'), 
                                              get(paste0(region, '_final_genes')),
                                              adult_title,
                                              get(paste0(region, '_umap_cols_recode')), 
                                              get(paste0(region, '_vln_cols_recode')))

saveRDS(umap_vln_plots_rna, paste0(R_dir, '03seurat_umap_vln_plt_rna_', region, '.rds'))

# Plot dystonia genes
dystonia_plot_rna <- create_stacked_vln_plot(seurat_object, 
                                             'cellIDs', 
                                             dystonia_genes,
                                             toupper(region),
                                             get(paste0(region, '_vln_cols_recode')))

## Find differential expressed marker genes in clusters  ------------------------------
message('\nCalculating diff exp markers ...\n')
# Idents(seurat_object) <- 'cellIDs'
# marker_genes <- FindAllMarkers(seurat_object)
# readr::write_tsv(marker_genes, paste0(R_dir, region, '_marker_genes.tsv'))

## Back to sketch
message('\nChanging to sketch object ...\n')
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

# bergmann <- c('NPY', 'TNC', 'LINC01727', 'FST', 'MT2A', 'PIFO', 'RSPH1')
# kozareva <- c('PPP1R17', 'GABRA6', 'EOMES', 'LYPD6', 'PRKCD', 'SORC3', 
#               'PTPRK', 'PRKCD', 'NXPH1', 'CDH22', 'KLHL1', 'ALDH1A3', 'SLC6A5', 'HTR2A', 'EDIL3',
#               'DCN', 'KCNJ8', 'MRC1', 'FIT1', 'FOXJ1', 'SLC6A5', 'GRM2', 'SST', 'PTPRC')
# leuko <- c("PTPRC", "SKAP1", "ARHGAP15", "PRKCH", "IKZF1", "STAT4", "DOCK8", 
#            "CD247", "TC2N", "IQGAP2", "FYB1", "SAMD3", "BCL11B", "CARD11", 
#            "EMB", "ETS1", "HLA-E", "LCP1", "CD96", "THEMIS", "STK17B", "APBB1IP", 
#            "IKZF3", "TNFAIP8", "CLEC2D", "GNG2", "CCL5", "CD53", "FLI1", 
#            "ZC3HAV1")
# 
# dput(read_tsv('~/Desktop/dystonia_snRNAseq_2024/results/01R_objects/cer_marker_genes.tsv') %>%
#   filter(cluster == 'Cer-adult-Leuko?') %>%
#   slice_head(n = 30) %>%
#     pull(gene))
#   
# 
# seurat_object@meta.data |>
#   as_tibble() 
