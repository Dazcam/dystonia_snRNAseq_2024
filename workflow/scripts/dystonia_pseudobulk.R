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
for (region in fetal_regions) {

  seurat_fetal <- readRDS(paste0(fetal_dir, 'seurat.', region, '.final.rds'))
  Idents(seurat_fetal) <- 'cellIDs'
  aggr <- calculate_aggregated_expression(seurat_fetal, dystonia_genes)
  aver <- calculate_average_expression(seurat_fetal, dystonia_genes)
  assign(paste0('aggr_ftl_', region), aggr)
  assign(paste0('aver_ftl_', region), aver)
  assign(paste0('seurat_ftl_', region), seurat_fetal)
  
}

## Section for Hawk as can't get run agg or aver if BPCells object was created online. ----
# Calcuate average and aggreagte expression: relys on home dir being correct in BPCell object
# message('Reading ', region)
# seurat_object <- readRDS(paste0(R_dir, 'ann_seurat_', region, '.rds'))
# message('Recoding ...')
# seurat_object[[paste0(region, '_clusters')]] <- recode_cluster_ids(seurat_object, region, 'harmony_clusters_0.1')
# message('Setting Idents ...')
# Idents(seurat_object) <- paste0(region, '_clusters')
# seurat_object
# aver_exp_mat <- calculate_average_expression(seurat_object, dystonia_genes)
# aggr_exp_mat <- calculate_aggregated_expression(seurat_object, dystonia_genes)
# message('Saving ...')
# saveRDS(object = aver_exp_mat, file = paste0(R_dir, "aver_exp_", region, ".rds"))
# saveRDS(object = aggr_exp_mat, file = paste0(R_dir, "aggr_exp_", region, ".rds"))

## ------------------------------------------------------------------------------------


  



# -------------------------------------------------------
tbl_list <- list()
for (region in regions) {
  
  for (age in c('adult', 'fetal')) {
    
    for (test in c('aggr', 'aver')) {
      
      if (age == 'adult' & region == 'cer') {
        
        # Add empty row of TH genes for Cer
        gene_mat <- readRDS(paste0(R_dir, "seurat_", test, "_", age, "_", region, ".rds")) 
        th_mat <- as(MatrixExtra::emptySparse(nrow = 1, ncol = 15), "dgCMatrix")
        rownames(th_mat) <- 'TH'
        
        gene_tbl <- rbind(gene_mat, th_mat) |>
          as_tibble(rownames = 'gene') |>
          arrange(gene)
        
      } else {
    
      gene_tbl <- readRDS(paste0(R_dir, "seurat_", test, "_", age, "_", region, ".rds")) |>
        as_tibble(rownames = 'gene') |>
        arrange(gene)
      
      if (age == 'fetal') {
        
        gene_tbl |>
          rename_with(~str_replace(., "_", "-fetal-"), contains('_'))
        
      }
      
      }
      
      tbl_list[[paste0(age, '_', test, '_', region)]] <- gene_tbl
    
    }
    
  }
  
}

writexl::write_xlsx(tbl_list, "~/Desktop/dystonia_snRNAseq_2024/results/01R_objects/aggr_exp_all_regions.xlsx")


##  Join data
ag_exp_mat <- cbind(aggr_ftl_cer, aggr_ftl_wge, aggr_ftl_pfc, aggr_adult_cer, aggr_adult_str)
av_exp_mat <- cbind(aver_ftl_cer, aver_ftl_wge, aver_ftl_pfc, aver_adult_cer, aver_adult_str)

# Scale data
av_exp_mat <- normalizeCounts(av_exp_mat)
av_exp_mat_scaled <- t(apply(av_exp_mat, 1, scale)) # centre and scale each col (Z-Score) then t
colnames(av_exp_mat_scaled) <- colnames(av_exp_mat)

# Normalise??
ag_exp_mat <- normalizeCounts(ag_exp_mat)
ag_exp_mat_scaled <- t(apply(ag_exp_mat, 1, scale)) # centre and scale each col (Z-Score) then t
colnames(ag_exp_mat_scaled) <- colnames(ag_exp_mat)

# Generate Heatmaps - huge in size
jpeg(paste0('~/Desktop/dystonia_av_exp.jpg'), width = 960, height = 2500,
     units = "px", pointsize = 12, quality = 150)
ComplexHeatmap::Heatmap(t(av_exp_mat_scaled),
                        column_labels = rownames(av_exp_mat_scaled),
                        cluster_columns = F,
                        cluster_rows = T,
                        width = ncol(t(av_exp_mat_scaled))*unit(5, "mm"),
                        height = nrow(t(av_exp_mat_scaled))*unit(5, "mm"))

dev.off()

jpeg(paste0('~/Desktop/dystonia_ag_exp_scaled.jpg'), width = 960, height = 2500,
     units = "px", pointsize = 12, quality = 150)
ComplexHeatmap::Heatmap(t(ag_exp_mat_scaled),
                        column_labels = rownames(ag_exp_mat_scaled),
                        cluster_columns = F,
                        cluster_rows = T,
                        width = ncol(t(ag_exp_mat_scaled))*unit(5, "mm"),
                        height = nrow(t(ag_exp_mat_scaled))*unit(5, "mm"))
dev.off()


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

