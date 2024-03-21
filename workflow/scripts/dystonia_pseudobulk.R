#--------------------------------------------------------------------------------------
#
#    Dystonia - Pseudobulk - Compare average and aggregate expression
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Script to compare expression of dystonia genes across all regions / ages
#  Compare average and aggregate expression of 25 dystonia genes across all regions

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

# Fetal
# fetal_regions <- c('pfc', 'cer', 'wge')
# fetal_dir <- '~/Desktop/fetal_brain_snRNAseq_110122/resources/R_objects/'
# for (region in fetal_regions) {
#   
#   seurat_fetal <- readRDS(paste0(fetal_dir, 'seurat.', region, '.final.rds'))
#   Idents(seurat_fetal) <- 'cellIDs'
#   aggr <- calculate_aggregated_expression(seurat_fetal, paste0(region, '_ftl'), dystonia_genes)
#   aver <- calculate_average_expression(seurat_fetal, paste0(region, '_ftl'), dystonia_genes)
#   assign(paste0('aggr_ftl_', region), aggr)
#   assign(paste0('aver_ftl_', region), aver)
#   
# }

## Section for Hawk as can't get run agg or aver if BPCells object was created online. ----
# Calcuate average and aggreagte expression: relys on home dir being correct in BPCell object
message('Reading ', region)
seurat_object <- readRDS(paste0(R_dir, 'ann_seurat_', region, '.rds'))
message('Recoding ...')
seurat_object[[paste0(region, '_clusters')]] <- recode_cluster_ids(seurat_object, region, 'harmony_clusters_0.1')
message('Setting Idents ...')
Idents(seurat_object) <- paste0(region, '_clusters')
seurat_object
aver_exp_mat <- calculate_average_expression(seurat_object, paste0(region, '_adult'), dystonia_genes)
aggr_exp_mat <- calculate_aggregated_expression(seurat_object, paste0(region, '_adult'), dystonia_genes)
message('Saving ...')
saveRDS(object = aver_exp_mat, file = paste0(R_dir, "seurat_aver_exp_", region, ".rds"))
saveRDS(object = aggr_exp_mat, file = paste0(R_dir, "seurat_aggr_exp_", region, ".rds"))

## ------------------------------------------------------------------------------------

# Get gene lengths
# mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# lookup <- biomaRt::getBM(mart = mart, 
#                          attributes=c("ensembl_gene_id", 
#                                       "external_gene_name", 
#                                       "start_position", 
#                                       "end_position"))
# gene_lengths <- as_tibble(dystonia_genes) %>%
#   left_join(lookup, by = join_by(value == external_gene_name)) %>%
#   mutate(length = end_position - start_position)
  



# -------------------------------------------------------

# Adult data
# seurat_object <- readRDS(paste0(R_dir, 'seurat_', region, '.rds'))
# 
# # Join layers - Need to do this before stacked vln plotting of diff expression
# seurat_object <- JoinLayers(seurat_object)
# 
# # Adult data
# for (region in c('fcx', 'str')) {
#   
#   aggr <- readRDS(paste0(R_dir, "seurat_aggr_exp_", region, ".Rds"))
#   aver <- readRDS(paste0(R_dir, "seurat_aver_exp_", region, ".Rds"))
#   
#   assign(paste0('aggr_adult_', region), aggr, envir = .GlobalEnv)
#   assign(paste0('aver_adult_', region), aver, envir = .GlobalEnv)
#   
# }
# 
# agg_exp_mat <- AggregateExpression(seurat_fetal, features = dystonia_genes)
# agg_exp_mat <- agg_exp_mat$RNA
# colnames(agg_exp_mat) <- paste0(prefix, colnames(agg_exp_mat))
# 
# ##  Join data 
# ag_exp_mat <- cbind(aggr_ftl_cer, aggr_ftl_wge, aggr_ftl_pfc, aggr_adl_fcx, aggr_adl_str)
# av_exp_mat <- cbind(aver_ftl_cer, aver_ftl_wge, aver_ftl_pfc, aver_adl_fcx, aver_adl_str)
# 
# # Scale data
# av_exp_mat <- normalizeCounts(av_exp_mat)
# av_exp_mat_scaled <- t(apply(av_exp_mat, 1, scale)) # centre and scale each col (Z-Score) then t
# colnames(av_exp_mat_scaled) <- colnames(av_exp_mat) 
# 
# # Normalise??
# ag_exp_mat <- normalizeCounts(ag_exp_mat)
# ag_exp_mat_scaled <- t(apply(ag_exp_mat, 1, scale)) # centre and scale each col (Z-Score) then t
# colnames(ag_exp_mat_scaled) <- colnames(ag_exp_mat) 
# 
# # Generate Heatmaps - huge in size
# jpeg(paste0('~/Desktop/dystonia_av_exp.jpg'), width = 960, height = 2500, 
#      units = "px", pointsize = 12, quality = 150)
# ComplexHeatmap::Heatmap(t(av_exp_mat_scaled), 
#                         column_labels = rownames(av_exp_mat_scaled), 
#                         cluster_columns = F, 
#                         cluster_rows = T,
#                         width = ncol(t(av_exp_mat_scaled))*unit(5, "mm"), 
#                         height = nrow(t(av_exp_mat_scaled))*unit(5, "mm"))
# 
# dev.off()
# 
# jpeg(paste0('~/Desktop/dystonia_ag_exp.jpg'), width = 960, height = 2500, 
#      units = "px", pointsize = 12, quality = 150)
# ComplexHeatmap::Heatmap(t(ag_exp_mat_scaled), 
#                         column_labels = rownames(ag_exp_mat_scaled), 
#                         cluster_columns = F, 
#                         cluster_rows = T,
#                         width = ncol(t(av_exp_mat_scaled))*unit(5, "mm"), 
#                         height = nrow(t(av_exp_mat_scaled))*unit(5, "mm"))
# dev.off()


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

k
