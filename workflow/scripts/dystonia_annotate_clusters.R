#--------------------------------------------------------------------------------------
#
#    Dystonia - Annotate clusters
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# Needs to be checked interactively and coordinated with gene and colour lists in 
# dystonia_gene_lists.R

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
seurat_object <- readRDS(paste0(R_dir, 'seurat_', region, '.rds'))

# Join layers - Need to do this before stacked vln plotting of diff expression
seurat_object <- JoinLayers(seurat_object)

# seurat_object <- RunUMAP(seurat_object, reduction = "harmony.full", 
#                          dims = 1:pc_thresh, reduction.name = "umap.full",
#                          reduction.key = "UMAPfull_")

if (region %in% c('str', 'cer')) {

  # Recode cluster IDs - sketch object only for now
  seurat_object[[paste0(region, '_clusters')]] <- recode_cluster_ids(seurat_object, region, 'harmony_clusters_0.1')
  
  # Plot paired vln
  vln_plots <- plot_paired_vlns(seurat_object, paste0(region, '_clusters'), general_genes,
                                get(paste0(region, '_genes')), get(paste0(region, '_vln_cols_recode')))
  
  # Plot paired umap and vln
  umap_vln_plots <- plot_paired_umap_vln(seurat_object, paste0(region, '_clusters'), get(paste0(region, '_final_genes')), 
                                         get(paste0(region, '_umap_cols_recode')), get(paste0(region, '_vln_cols_recode')))
  
  dystonia_plot <- create_stacked_vln_plot(seurat_object, paste0(region, '_clusters'), dystonia_genes,
                                           toupper(region), get(paste0(region, '_vln_cols_recode')))
  
  # Find differential expressed marker genes in clusters
  message('\nCalculating diff exp markers ...')
  Idents(seurat_object) <- paste0(region, '_clusters')
  marker_genes <- FindAllMarkers(seurat_object)
  readr::write_tsv(marker_genes, paste0(R_dir, region, '_marker_genes.tsv'))
  
  # Project data to whole object - # Crashes locally with Cer and FCX
  seurat_object[["sketch"]] <- split(seurat_object[["sketch"]], f = seurat_object[[sample_split]])
  seurat_object <- project_sketch_data(seurat_object,
                                       pc_thresh,
                                       'harmony',
                                       'umap.harmony',
                                       paste0(region, '_clusters'))
  
  # Calculate aggr and aver expr for regional comparison of dystonia gene expression
  message('Set Idents to:', paste0(region, '_clusters'), ' ...\n')
  Idents(seurat_object) <- paste0(region, '_clusters')
  
  message('Set default assay to RNA ...\n')
  DefaultAssay(seurat_object) <- "RNA"
  
  message('Count NAs in ', paste0(region, '_clusters'), ' in RNA assay: ', sum(is.na(seurat_object$cer_clusters)), '\n')
  
  # Calcuate average and aggreagte expression: relys on home dir being correct in BPCell object
  aver_exp_mat <- calculate_average_expression(seurat_object, paste0(region, '_adult'), dystonia_genes)
  aggr_exp_mat <- calculate_aggregated_expression(seurat_object, paste0(region, '_adult'), dystonia_genes)
  
  # Save objects
  saveRDS(object = aver_exp_mat, file = paste0(R_dir, "seurat_aver_exp_", region, ".Rds"))
  saveRDS(object = aggr_exp_mat, file = paste0(R_dir, "seurat_aggr_exp_", region, ".Rds"))
  
  ## Create markdown doc  ---------------------------------------------------------------
  rmarkdown::render(markdown_ann_doc, output_file = markdown_ann_html, output_dir = R_dir)
  
}

message('Set default assat to RNA ...\n')
DefaultAssay(seurat_object) <- "sketch"

# Join layers - Need to do this before stacked vln plotting of diff expression
seurat_object <- JoinLayers(seurat_object)

saveRDS(seurat_object, paste0(R_dir, 'ann_seurat_', region, '.rds'))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------



# # Code for working out cell types genes / cols
# vln_gen <- create_stacked_vln_plot(seurat_object, 'harmony_clusters_0.3', general_genes,
#                                    toupper(region))
# clust <- DimPlot_scCustom(seurat_object, group.by = 'harmony_clusters_0.3', 
#                           alpha = 0.1, pt.size = 0.1, label = T) +
#   NoLegend()
# stiletti <- DimPlot_scCustom(seurat_object, group.by = 'cell_type', 
#                              alpha = 0.1, pt.size = 0.1,
#                              repel = T) 
# vln_gen | clust / stiletti
# 
# plot_paired_vlns(seurat_object, 'harmony_clusters_0.1', general_genes,
#                  get(paste0(region, '_genes')))
# 
# 
# 
# 
# # Cluster counts - Need to change to full dataset first data first
# seurat_object$seurat_clusters
# 
# # ### NOTES ON PRELIMINARY CLUST LABELS  ------------------------------------------------
# if (region == 'str') {
# 
#   # Striatum
#   # 13 and 25 experesses GAD 1 and the OLIGS
#   # 22 does not express GAD but exp InN transporters
#   clusters_id_recode <- str_clusters_recode
#   vln_clr_recode <- str_vln_cols_recode
#   umap_clr_recode <- str_umap_cols_recode
# 
# } else if (region == 'cer') {
# 
#   pc_thresh <- 30
#   sample_split <- 'orig.ident'
# 
# } else {
# 
#   # FCX
#   seurat_sk_fcx$harmony_clusters_recode <- fcx_clusters_recode
#   umap_recode <- DimPlot(seurat_sk_fcx, reduction = 'umap.harmony',
#                          group.by = 'harmony_clusters_recode',
#                          cols = fcx_umap_cols_recode)
#   vln_plot_recode <- create_stacked_vln_plot(seurat_sk_fcx, 'harmony_clusters_recode',
#                                              fcx_genes, toupper(region),
#                                              fcx_vln_cols_recode)
# 
#   # Dystonia gene check
#   vln_plot_recode_dyst <- create_stacked_vln_plot(seurat_sk_fcx, 'harmony_clusters_recode', dystonia_genes,
#                                                   toupper(region), fcx_vln_cols_recode)
# 
# }
# 
# # Plot umap and vlns with cell types and colours specified
# # Recode cell IDs
# seurat_object$harmony_clusters_recode <- clusters_id_recode
# 
# # UMAP with recode cell IDs
# umap_recode <- DimPlot(seurat_object, reduction = 'umap.harmony',
#                        group.by = 'harmony_clusters_recode',
#                        cols = umap_clr_recode)
# 
# # Region specific gene check vln plot
# vln_plot_recode <- create_stacked_vln_plot(seurat_object, 'harmony_clusters_recode',
#                                            general_genes, toupper(region),
#                                            vln_clr_recode)
# # Dystonia gene check vln plot
# vln_plot_recode_dyst <- create_stacked_vln_plot(seurat_object, 'harmony_clusters_recode',
#                                                 dystonia_genes, toupper(region),
#                                                 vln_clr_recode)
# 
# 
# 
# # Calculate aggr and aver expr for regional comparison of dystonia gene expression
# aver_exp_mat <- calculate_average_expression(seurat_object, paste0(region, '_adult'), dystonia_genes)
# aggr_exp_mat <- calculate_aggregated_expression(seurat_object, paste0(region, '_adult'), dystonia_genes)
# 
# saveRDS(object = aver_exp_mat, file = paste0(R_dir, "seurat_aver_exp_", region, ".Rds"))
# saveRDS(object = aggr_exp_mat, file = paste0(R_dir, "seurat_aggr_exp_", region, ".Rds"))
# 
# test <- readRDS(paste0(R_dir, "seurat_aggr_exp_", region, ".Rds"))