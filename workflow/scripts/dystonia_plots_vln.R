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
adult_object <- readRDS(paste0(R_dir, '02seurat_', region, '.rds'))
fetal_object <- readRDS(paste0(fetal_dir, 'seurat.', fetal_region, '.final.rds'))

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
  ggtitle("Fetal Ganglionic Eminences") +
  facet_wrap(~feature,  ncol = 1, strip.position = "left") +
  scale_y_continuous(position = "right", limits=c(-0.00004, 5), breaks = 4)


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
adult_object[[paste0(region, '_clusters')]] <- recode_cluster_ids(adult_object, 
                                                                  region, 
                                                                  'harmony_clusters_0.1')

dystonia_adult_plot_sketch <- create_stacked_vln_plot(adult_object,
                                                paste0(region, '_clusters'), 
                                                dystonia_genes,
                                                "Adult Striatum", 
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

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
