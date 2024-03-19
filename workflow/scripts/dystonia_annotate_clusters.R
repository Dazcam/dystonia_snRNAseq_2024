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

seurat_object <- RunUMAP(seurat_object, reduction = "harmony.full", 
                         dims = 1:pc_thresh, reduction.name = "umap.full",
                         reduction.key = "UMAPfull_")

# Join layers - Need to do this before stacked vln plotting of diff expression
seurat_object <- JoinLayers(seurat_object)

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
Idents(seurat_object) <- 'harmony_clusters_0.1'
marker_genes <- FindAllMarkers(seurat_object)




## Create markdown doc  ---------------------------------------------------------------
rmarkdown::render(markdown_ann_doc, output_file = markdown_ann_html, output_dir = R_dir)

vln_reg <- create_stacked_vln_plot(seurat_object, 'harmony_clusters_0.1', str_genes,
                        toupper(region), str_vln_cols_recode)
vln_gen <- create_stacked_vln_plot(seurat_object, 'harmony_clusters_0.1', general_genes,
                               toupper(region), str_vln_cols_recode)
vln_str <- create_stacked_vln_plot(seurat_object, 'str_clusters', final_genes,
                                   toupper(region), str_vln_cols_recode)
clust <- DimPlot_scCustom(seurat_object, group.by = 'harmony_clusters_0.1', 
                          alpha = 0.3, pt.size = 0.1, label = T) +
  NoLegend()
stiletti <- DimPlot_scCustom(seurat_object, group.by = 'cell_type', 
                             alpha = 0.3, pt.size = 0.1,
                             repel = T) 
vln_gen | clust / stiletti







# Cluster counts - Need to change to full dataset first data first
seurat_object$seurat_clusters

# ### NOTES ON PRELIMINARY CLUST LABELS  ------------------------------------------------
if (region == 'str') {

  # Striatum
  # 13 and 25 experesses GAD 1 and the OLIGS
  # 22 does not express GAD but exp InN transporters
  clusters_id_recode <- str_clusters_recode
  vln_clr_recode <- str_vln_cols_recode
  umap_clr_recode <- str_umap_cols_recode

} else if (region == 'cer') {

  pc_thresh <- 30
  sample_split <- 'orig.ident'

} else {

  # FCX
  seurat_sk_fcx$harmony_clusters_recode <- fcx_clusters_recode
  umap_recode <- DimPlot(seurat_sk_fcx, reduction = 'umap.harmony',
                         group.by = 'harmony_clusters_recode',
                         cols = fcx_umap_cols_recode)
  vln_plot_recode <- create_stacked_vln_plot(seurat_sk_fcx, 'harmony_clusters_recode',
                                             fcx_genes, toupper(region),
                                             fcx_vln_cols_recode)

  # Dystonia gene check
  vln_plot_recode_dyst <- create_stacked_vln_plot(seurat_sk_fcx, 'harmony_clusters_recode', dystonia_genes,
                                                  toupper(region), fcx_vln_cols_recode)

}

# Plot umap and vlns with cell types and colours specified
# Recode cell IDs
seurat_object$harmony_clusters_recode <- clusters_id_recode

# UMAP with recode cell IDs
umap_recode <- DimPlot(seurat_object, reduction = 'umap.harmony',
                       group.by = 'harmony_clusters_recode',
                       cols = umap_clr_recode)

# Region specific gene check vln plot
vln_plot_recode <- create_stacked_vln_plot(seurat_object, 'harmony_clusters_recode',
                                           general_genes, toupper(region),
                                           vln_clr_recode)
# Dystonia gene check vln plot
vln_plot_recode_dyst <- create_stacked_vln_plot(seurat_object, 'harmony_clusters_recode',
                                                dystonia_genes, toupper(region),
                                                vln_clr_recode)


### END NOTES ON PRELIMINARY CLUST LABELS  ------------------------------------------------

## Plot gene expression violin plot  ------------------------------------------------------
# Specific gene check
vln_plot_specific <- create_stacked_vln_plot(seurat_sk_fcx, 'harmony_clusters',
                                             general_genes, toupper(region),
                                             fcx_vln_cols_recode)

seurat_sk_fcx$harmony_clusters_recode <- fcx_clusters_recode

# Gene check
vln_plot_recode <- create_stacked_vln_plot(seurat_sk_fcx, 'harmony_clusters_recode',
                                           general_genes, toupper(region), fcx_vln_cols_recode)






hip_colours <- c('#76B5C5', '#B200ED', '#FAA0A0', '#EF0029', '#CEE5FD',
                 '#95D840FF', "#00BDD2", '#00FF00A5', "#DCBEFF", '#10A53DFF',
                 '#6F2DA8', '#ABDBE3', '#1E81B0', '#D2042D', '#006400',
                 '#FDE725FF', '#779CBA', '#F58231', '#9A6324')

str_umap_cols_recode <- c("Str-adult-InN-1" = '#3CBB75FF', "Str-adult-InN-2" = '#31C53F', "Str-adult-InN-3" = '#708238',
                          "Str-adult-InN-4" = '#B7FFB7', "Str-adult-InN-5" = '#006400', "Str-OPC-1" = '#FDE725FF',
                          "Str-adult-InN-6" = '#95D840FF', "Str-adult-InN-7" = '#2FF18B', "Str-adult-InN-8" = '#9DC183',
                          "Str-adult-Ast-1" = '#FF5959', "Str-adult-InN-9" = '#3CBB75FF', "Str-adult-InN-10" = '#31C53F',
                          "Str-adult-InN-11" = '#708238', "Str-OPC-2" = '#FDE725FF', "Str-adult-InN-12" = '#006400',
                          "Str-adult-InN-13" = '#95D840FF', "Str-adult-ExN-1" = '#00B6EB', "Str-adult-MG" = '#F58231',
                          "Str-adult-InN-14" = '#2FF18B', "Str-adult-InN-15" = '#9DC183', "Str-adult-ExN-2" = '#00B6EB',
                          "Str-adult-InN-16" = '#0CB702', "Str-adult-InN-17" = '#00BE67', "Str-adult-Ast-2" = '#FF5959',
                          "Str-adult-InN-18" = '#7CAE00', "Str-OPC-3" = '#FDE725FF')

str_vln_cols_recode <- c("Str-adult-InN-1" = '#3CBB75FF', "Str-adult-InN-2" = '#3CBB75FF', "Str-adult-InN-3" = '#3CBB75FF',
                         "Str-adult-InN-4" = '#3CBB75FF', "Str-adult-InN-5" = '#3CBB75FF', "Str-OPC-1" = '#FDE725FF',
                         "Str-adult-InN-6" = '#3CBB75FF', "Str-adult-InN-7" = '#3CBB75FF', "Str-adult-InN-8" = '#3CBB75FF',
                         "Str-adult-Ast-1" = '#FF5959', "Str-adult-InN-9" = '#3CBB75FF', "Str-adult-InN-10" = '#3CBB75FF',
                         "Str-adult-InN-11" = '#3CBB75FF', "Str-OPC-2" = '#FDE725FF', "Str-adult-InN-12" = '#3CBB75FF',
                         "Str-adult-InN-13" = '#3CBB75FF', "Str-adult-ExN-1" = '#00B6EB', "Str-adult-MG" = '#F58231',
                         "Str-adult-InN-14" = '#3CBB75FF', "Str-adult-InN-15" = '#3CBB75FF', "Str-adult-ExN-2" = '#00B6EB',
                         "Str-adult-InN-16" = '#3CBB75FF', "Str-adult-InN-17" = '#3CBB75FF', "Str-adult-Ast-2" = '#FF5959',
                         "Str-adult-InN-18" = '#3CBB75FF', "Str-OPC-3" = '#FDE725FF')



# Calculate aggr and aver expr for regional comparison of dystonia gene expression
aver_exp_mat <- calculate_average_expression(seurat_object, paste0(region, '_adult'), dystonia_genes)
aggr_exp_mat <- calculate_aggregated_expression(seurat_object, paste0(region, '_adult'), dystonia_genes)

saveRDS(object = aver_exp_mat, file = paste0(R_dir, "seurat_aver_exp_", region, ".Rds"))
saveRDS(object = aggr_exp_mat, file = paste0(R_dir, "seurat_aggr_exp_", region, ".Rds"))

test <- readRDS(paste0(R_dir, "seurat_aggr_exp_", region, ".Rds"))