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
#R_dir <- '../results/01R_objects/'
#region <- 'cer'
#fetal_dir <- '../resources/public_data/cameron_2023/'
adult_object <- readRDS(paste0(R_dir, '02seurat_', region, '.rds'))
fetal_object <- readRDS(paste0(fetal_dir, 'seurat_', fetal_region, '.rds'))

# Fetal -----
fetal_object$cellIDs <- stringr::str_replace_all(fetal_object$cellIDs,
                                                 paste0(region_recode, "-"), 
                                                 paste0(region_recode, "-fetal-"))

if (region == "fcx" || region == "cer") {
  
  Idents(fetal_object) <- factor(names(get(paste0('fetal_', region, '_vln_recode'))), 
                                 levels = names(get(paste0('fetal_', region, '_vln_recode'))))
  plot_cols <- get(paste0('fetal_', region, '_vln_recode'))

  } else {
    
    Idents(fetal_object) <- factor(names(fetal_ge_vln_recode), 
                                   levels = names(fetal_ge_vln_recode))
    plot_cols <- fetal_ge_vln_recode
}

fetal_plot <- VlnPlot(fetal_object, dystonia_genes, stack = TRUE, flip = TRUE,  
                      cols = plot_cols,
                      same.y.lims = TRUE, fill.by = 'ident') +
  theme(legend.position = "none",
        plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        axis.text.x  = element_text(colour = "#000000", size = 12),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 16)) +
  ggtitle(fetal_title) +
  facet_wrap(~feature,  ncol = 1, strip.position = "left") +
  scale_y_continuous(position = "right", limits=c(-0.00004, 5), breaks = 4)

#layer_data(fetal_plot, 1) # Pull info from plot
# Adult -----
# Subset seurat object and add empty row for gene TH
# VlnPlot(adult_object, dystonia_genes, stack = TRUE, flip = TRUE,  
#         cols = get(paste0(region, '_vln_cols_recode')),
#         same.y.lims = TRUE, fill.by = 'ident') +
#   facet_wrap(feature ~ ., drop = F) +
#   theme(strip.text.y.left = element_text(angle = 0, size = 16)) +
#   scale_y_continuous(position = "right", limits=c(-0.00004, 5), breaks = 4)
# 
# th_mat <- as(MatrixExtra::emptySparse(nrow = 1, ncol = 15), "dgCMatrix")
# rownames(th_mat) <- 'TH'
# aggr_adult_cer <- rbind(aggr_adult_cer, th_mat)
# 
# str(dystonia_adult_plot_sketch$facet$params$drop)
# adult_plot$facet$params$drop = FALSE
# adult_2 <- subset(adult_object, features = dystonia_genes)


# Recode cluster IDs - sketch object 
message("Converting BPCells counts and data matrix to in memory format for hdWGCNA ...")
DefaultAssay(adult_object) <- 'RNA'
adult_object <- JoinLayers(adult_object)
adult_object[["RNA"]]$counts <- as(object = adult_object[["RNA"]]$counts, Class = "dgCMatrix")
adult_object[["RNA"]]$data <- as(object = adult_object[["RNA"]]$data, Class = "dgCMatrix")

message("Recode cluster IDs ... ")
adult_object$cellIDs <- recode_cluster_ids(adult_object, region, 'cluster_full')
unique(adult_object$cellIDs)

message("Add blank row  cluster IDs ... ")
if (region == 'cer') {

  th_mat <- as(MatrixExtra::emptySparse(nrow = 1, ncol = ncol(adult_object)), "dgCMatrix")
  rownames(th_mat) <- 'TH'
  adult_object[["RNA"]]$data <- rbind(adult_object[["RNA"]]$data, th_mat)
}

adult_object[[paste0(region, '_clusters')]] <- recode_cluster_ids(adult_object, 
                                                                  region, 
                                                                  'harmony_clusters_0.1')
dystonia_adult_plot_RNA <- create_stacked_vln_plot(adult_object,
                                                paste0(region, '_clusters'), 
                                                dystonia_genes,
                                                adult_title, 
                                                get(paste0(region, '_vln_cols_recode')))

adult_plot <- dystonia_adult_plot_sketch +
  theme(legend.position = "none",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 12),
        axis.text.x  = element_text(colour = "#000000", size = 12),
        axis.text.y  = element_text(colour = "#000000", size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        strip.text.y.left = element_blank()) +
  scale_y_continuous(position = "right", limits=c(-0.00004, 5), breaks = 4) +
  facet_wrap(~feature,  ncol = 1, strip.position = "left", drop = F)


# Join y-axes
egg::ggarrange(fetal_plot, adult_plot, nrow = 1)

# Add summary stats
for (sum_stat in c('mean', 'median')) {
  
  fetal_stat_plot <- fetal_plot +
    stat_summary(fun = sum_stat, colour = "red", size = 2.5, # Adds text statistic
                 geom = "text", aes(label = round(after_stat(y), 2)),
                 position = position_nudge(x = 0.25, y = 2)) 
  #  stat_summary(fun = median, geom = "point",  # Adds median line
  #               size = 5, colour = "white", shape = 95) +
  
  adult_stat_plot <- adult_plot +
    stat_summary(fun = sum_stat, colour = "red", size = 2.5,
                 geom = "text", aes(label = round(after_stat(y), 2)),
                 position = position_nudge(x = 0.25, y = 2)) 
  #  stat_summary(fun = median, geom = "point", 
  #               size = 5, colour = "white", shape = 95) +
  
  # Join y-axes
  vln_stat_plt <- egg::ggarrange(fetal_stat_plot, adult_stat_plot, nrow = 1)
  
  ggsave(filename = paste0(R_dir, region, "_vln_plt_with_", sum_stat, ".svg"),
         plot = vln_stat_plt, 
         width = 40, 
         height = 30, 
         dpi = 300, 
         units = "cm")
} 


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
