#--------------------------------------------------------------------------------------
#
#    Dystonia - Create vln plots - adult fetal
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
# R_dir <- '../results/01R_objects/'
# region <- 'cer'
# fetal_dir <- '../resources/public_data/cameron_2023/'

# Run fetal plots locally
if (Sys.info()[["nodename"]] == "Darrens-iMac-2.local") {

  # Fetal -----
  for (brain_region in c('ge_fetal', 'cer_fetal', 'fcx_fetal')) {
    
    fetal_object <- readRDS(paste0(fetal_dir, 'seurat_', brain_region, '.rds'))
  
    fetal_object$cellIDs <- stringr::str_replace_all(fetal_object$cellIDs,
                                                     paste0(region_recode, "-"), 
                                                     paste0(region_recode, "-fetal-"))
    assign(paste0(brain_region, '_obj'), fetal_object)
    rm(fetal_object)
    
  }
  
  # Set factors and colours
  for (brain_region in c('ge_fetal', 'cer_fetal', 'fcx_fetal')) {
    
    local_region <- str_replace(brain_region, '_fetal', '')
    fetal_object <- get(paste0(brain_region, '_obj'))
    
  
    if (brain_region == "fcx_fetal" || brain_region == "cer_fetal") {
      
      Idents(fetal_object) <- factor(names(get(paste0('fetal_', local_region, '_vln_recode'))), 
                                     levels = names(get(paste0('fetal_', local_region, '_vln_recode'))))
      plot_cols <- get(paste0('fetal_', local_region, '_vln_recode'))
      
      assign(paste0(brain_region, '_obj'), fetal_object)
      assign(paste0(brain_region, '_cols'), plot_cols)
      
      } else {
        
        Idents(fetal_object) <- factor(names(fetal_ge_vln_recode), 
                                       levels = names(fetal_ge_vln_recode))
        plot_cols <- fetal_ge_vln_recode
        
        assign(paste0(brain_region, '_obj'), fetal_object)
        assign(paste0(brain_region, '_cols'), plot_cols)
        
      }
    
    rm(fetal_object, plot_cols)
    
  }
  
  lhs_plot <- VlnPlot(fcx_fetal_obj, dystonia_genes, stack = TRUE, flip = TRUE,  
                      cols = fcx_fetal_cols,
                      same.y.lims = TRUE, fill.by = 'ident') +
    theme(legend.position = "none",
          plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16),
          text = element_text(size = 12),
          axis.text.x  = element_text(colour = "#000000", size = 12),
          axis.text.y  = element_blank(),
          axis.line.y.right = element_blank(),  # Removes the y-axis line on the right
          axis.ticks.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text.y.left = element_text(angle = 0, size = 16)) +
    ggtitle('Fetal Frontal Cortex') +
    facet_wrap(~feature,  ncol = 1, strip.position = "left") +
    scale_y_continuous(position = "right", limits=c(-0.00004, 5), breaks = 4)
  
  mid_plot <- VlnPlot(cer_fetal_obj, dystonia_genes, stack = TRUE, flip = TRUE,  
                      cols = cer_fetal_cols,
                      same.y.lims = TRUE, fill.by = 'ident') + 
    theme(legend.position = "none",
          plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16),
          axis.text.x  = element_text(colour = "#000000", size = 12),
          axis.text.y  = element_blank(),  # Hides y-axis text
          axis.ticks.y = element_blank(),  # Hides y-axis ticks
          axis.line.y.right = element_blank(),  # Removes the y-axis line on the right
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          strip.text = element_blank(),
          strip.background = element_blank()) +  
    ggtitle('Fetal Cerebellum') +
    facet_wrap(~feature, ncol = 1, strip.position = "left") +
    scale_y_continuous(position = "right", limits = c(-0.00004, 5), breaks = 4)
  
  rhs_plot <- VlnPlot(ge_fetal_obj, dystonia_genes, stack = TRUE, flip = TRUE,  
                      cols = ge_fetal_cols,
                      same.y.lims = TRUE, fill.by = 'ident') +
    theme(legend.position = "none",
          plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16),
          text = element_text(size = 12),
          axis.text.x  = element_text(colour = "#000000", size = 12),
          axis.text.y  = element_text(colour = "#000000", size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          strip.text.y.left = element_blank()) +
    ggtitle('Fetal Ganglionic Eminences') +
    facet_wrap(~feature,  ncol = 1, strip.position = "left", drop = F) +
    scale_y_continuous(position = "right", limits=c(-0.00004, 5), breaks = 4)
    
  # Join y-axes
  egg::ggarrange(lhs_plot, mid_plot, rhs_plot, nrow = 1)

  # Create adult plots on Hawk
  } else {

  ## Adult -------
  adult_object <- readRDS(paste0(R_dir, '02seurat_', region, '.rds'))
  
  # Recode cluster IDs - sketch object 
  DefaultAssay(adult_object) <- 'RNA'
  adult_object <- JoinLayers(adult_object)
  
  message("Recode cluster IDs ... ")
  adult_object$cellIDs <- recode_cluster_ids(adult_object, region, 'cluster_full')
  unique(adult_object$cellIDs)
  
  # Set plot colours
  plot_cols <- get(paste0(region, '_vln_cols_recode'))
  
  if (region == 'fcx') {
    
    message('Plotting ', region, ' plot ...')
    lhs_plot <- VlnPlot(adult_object, dystonia_genes, stack = TRUE, flip = TRUE,  
                        cols = plot_cols,
                        same.y.lims = TRUE, fill.by = 'ident') +
      theme(legend.position = "none",
            plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 16),
            text = element_text(size = 12),
            axis.text.x  = element_text(colour = "#000000", size = 12),
            axis.text.y  = element_blank(),
            axis.line.y.right = element_blank(),  # Removes the y-axis line on the right
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.text.y.left = element_text(angle = 0, size = 16)) +
      ggtitle('Fetal Frontal Cortex') +
      facet_wrap(~feature,  ncol = 1, strip.position = "left") +
      scale_y_continuous(position = "right", limits=c(-0.00004, 5), breaks = 4)
    
    saveRDS(lhs_plot, paste0(R_dir, region, '_vln_plot.rds'))}
  
  if (region == 'cer') {
    
    message("Add blank row to cer Seurat object for TH gene ... ")
    message("Converting BPCells counts and data matrix to in memory format first ...")
    
    message("Count matrix ...")
    adult_object[["RNA"]]$counts <- as(object = adult_object[["RNA"]]$counts, Class = "dgCMatrix")
    
    message("Data matrix ...")
    adult_object[["RNA"]]$data <- as(object = adult_object[["RNA"]]$data, Class = "dgCMatrix")
    
    th_mat <- as(MatrixExtra::emptySparse(nrow = 1, ncol = ncol(adult_object)), "dgCMatrix")
    rownames(th_mat) <- 'TH'
    
    new_counts <- rbind(adult_object[["RNA"]]$counts, th_mat)
    new_data <- rbind(adult_object[["RNA"]]$data, th_mat)
    obj_meta <- adult_object@meta.data
    
    message("Counts old: ", nrow(adult_object[["RNA"]]$counts), " ; Counts new: ", nrow(new_counts))
    message("Data old: ", nrow(adult_object[["RNA"]]$data), " ; Data new: ", nrow(new_data))
    
    new_cer_obj <- CreateSeuratObject(counts = new_counts, meta.data = obj_meta)
    new_cer_obj[["RNA"]]$data <- new_data
    
    message('Plotting ', region, ' plot ...')
    mid_plot <- VlnPlot(cer_fetal_obj, dystonia_genes, stack = TRUE, flip = TRUE,  
                        cols = cer_fetal_cols,
                        same.y.lims = TRUE, fill.by = 'ident') + 
      theme(legend.position = "none",
            plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 16),
            axis.text.x  = element_text(colour = "#000000", size = 12),
            axis.text.y  = element_blank(),  # Hides y-axis text
            axis.ticks.y = element_blank(),  # Hides y-axis ticks
            axis.line.y.right = element_blank(),  # Removes the y-axis line on the right
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.text = element_blank(),
            strip.background = element_blank()) +  
      ggtitle('Cerebellum') +
      facet_wrap(~feature, ncol = 1, strip.position = "left") +
      scale_y_continuous(position = "right", limits = c(-0.00004, 5), breaks = 4)
    
    saveRDS(mid_plot, paste0(R_dir, region, '_vln_plot.rds'))}
  
  if (region == 'str') {
    
    message('Plotting ', region, ' plot ...')
  
    rhs_plot <- VlnPlot(ge_fetal_obj, dystonia_genes, stack = TRUE, flip = TRUE,  
                        cols = ge_fetal_cols,
                        same.y.lims = TRUE, fill.by = 'ident') +
      theme(legend.position = "none",
            plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 16),
            text = element_text(size = 12),
            axis.text.x  = element_text(colour = "#000000", size = 12),
            axis.text.y  = element_text(colour = "#000000", size = 12),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 14),
            strip.text.y.left = element_blank()) +
      ggtitle('Fetal Ganglionic Eminences') +
      facet_wrap(~feature,  ncol = 1, strip.position = "left", drop = F) +
      scale_y_continuous(position = "right", limits=c(-0.00004, 5), breaks = 4)
    
    saveRDS(rhs_plot, paste0(R_dir, region, '_vln_plot.rds'))}
  
}

# Join y-axes
#egg::ggarrange(fetal_plot, adult_plot, nrow = 1)

# # Add summary stats
# for (sum_stat in c('mean', 'median')) {
#   
#   fetal_stat_plot <- fetal_plot +
#     stat_summary(fun = sum_stat, colour = "red", size = 2.5, # Adds text statistic
#                  geom = "text", aes(label = round(after_stat(y), 2)),
#                  position = position_nudge(x = 0.25, y = 2)) 
#   #  stat_summary(fun = median, geom = "point",  # Adds median line
#   #               size = 5, colour = "white", shape = 95) +
#   
#   adult_stat_plot <- adult_plot +
#     stat_summary(fun = sum_stat, colour = "red", size = 2.5,
#                  geom = "text", aes(label = round(after_stat(y), 2)),
#                  position = position_nudge(x = 0.25, y = 2)) 
#   #  stat_summary(fun = median, geom = "point", 
#   #               size = 5, colour = "white", shape = 95) +
#   
#   # Join y-axes
#   vln_stat_plt <- egg::ggarrange(fetal_stat_plot, adult_stat_plot, nrow = 1)
#   
#   ggsave(filename = paste0(R_dir, region, "_vln_plt_with_", sum_stat, ".svg"),
#          plot = vln_stat_plt, 
#          width = 40, 
#          height = 30, 
#          dpi = 300, 
#          units = "cm")
# } 



#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
