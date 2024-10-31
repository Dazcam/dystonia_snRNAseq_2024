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
  
  lhs_plot <- VlnPlot(ge_fetal_obj, dystonia_genes, stack = TRUE, flip = TRUE,  
                      cols = ge_fetal_cols,
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
    ggtitle('Fetal Ganglionic Eminences') +
    facet_wrap(~feature,  ncol = 1, strip.position = "left") +
    scale_y_continuous(position = "right", limits=c(-0.00004, 5), breaks = 4)
  
  mid_plot <- VlnPlot(fcx_fetal_obj, dystonia_genes, stack = TRUE, flip = TRUE,  
                      cols = fcx_fetal_cols,
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
    ggtitle('Fetal Frontal Cortex') +
    facet_wrap(~feature, ncol = 1, strip.position = "left") +
    scale_y_continuous(position = "right", limits = c(-0.00004, 5), breaks = 4)
  
  rhs_plot <- VlnPlot(cer_fetal_obj, dystonia_genes, stack = TRUE, flip = TRUE,  
                      cols = cer_fetal_cols,
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
    ggtitle('Fetal Cerebellum') +
    facet_wrap(~feature,  ncol = 1, strip.position = "left", drop = F) +
    scale_y_continuous(position = "right", limits=c(-0.00004, 5), breaks = 4)
    
  # Join y-axes
  egg::ggarrange(lhs_plot, mid_plot, rhs_plot, nrow = 1)

  # Create adult plots on Hawk
  } else {

  ## Adult -------
  adult_object <- readRDS(paste0(R_dir, '03seurat_', region, '.rds'))
  
  # Recode cluster IDs - sketch object 
  message('\nChanging to RNA object ...\n')
  DefaultAssay(adult_object) <- 'RNA'
  adult_object <- JoinLayers(adult_object)
  message('\nSetting idents ...\n')
  Idents(adult_object) <- adult_object$cellIDs
  message('Number of NAs in Idents: ', anyNA(Idents(adult_object)))
  
  # message("Recode cluster IDs ... ")
  # adult_object$cellIDs <- recode_cluster_ids(adult_object, region, 'cluster_full')
  # Idents(adult_object) <- adult_object$cellIDs
  # message('Number of NAs in Idents: ', anyNA(Idents(adult_object)))
  
  # Set plot colours
  plot_cols <- get(paste0(region, '_vln_cols_recode'))
  
  if (region == 'fcx') {
    
    message('Plotting ', region, ' plot ...')
    # lhs_plot <- VlnPlot(adult_object, dystonia_genes, stack = TRUE, flip = TRUE,  
    #                     cols = plot_cols, same.y.lims = TRUE, fill.by = 'ident') +
    #   theme(legend.position = "none",
    #         plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
    #         panel.grid.major = element_blank(), 
    #         panel.grid.minor = element_blank(),
    #         plot.title = element_text(hjust = 0.5, size = 16),
    #         text = element_text(size = 12),
    #         axis.text.x  = element_text(colour = "#000000", size = 12),
    #         axis.text.y  = element_blank(),
    #         axis.line.y.right = element_blank(),  # Removes the y-axis line on the right
    #         axis.ticks.y = element_blank(),
    #         axis.title.x = element_blank(),
    #         axis.title.y = element_blank(),
    #         strip.text.y.left = element_text(angle = 0, size = 16)) +
    #   ggtitle('Frontal Cortex') +
    #   facet_wrap(~feature,  ncol = 1, strip.position = "left") +
    #   scale_y_continuous(position = "right", limits=c(-0.00004, 5), breaks = 4)
    
    lhs_plot <- create_stacked_vln_plot(adult_object, 
                                        'cellIDs', 
                                        dystonia_genes,
                                        adult_title, 
                                        get(paste0(region, '_vln_cols_recode')))
    
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
    new_meta <- adult_object@meta.data
    
    message("Counts old: ", nrow(adult_object[["RNA"]]$counts), " ; Counts new: ", nrow(new_counts))
    message("Data old: ", nrow(adult_object[["RNA"]]$data), " ; Data new: ", nrow(new_data))
    
    new_cer_obj <- CreateSeuratObject(counts = new_counts, meta.data = new_meta)
    new_cer_obj[["RNA"]]$data <- new_data
    
    Idents(new_cer_obj) <- new_cer_obj$cellIDs
    message('Number of NAs in Idents: ', anyNA(Idents(new_cer_obj)))
    
    message('Plotting ', region, ' plot ...')
    # mid_plot <- VlnPlot(new_cer_obj, dystonia_genes, stack = TRUE, flip = TRUE,  
    #                     cols = plot_cols, same.y.lims = TRUE, fill.by = 'ident') +
    #   theme(legend.position = "none",
    #         plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
    #         panel.grid.major = element_blank(), 
    #         panel.grid.minor = element_blank(),
    #         plot.title = element_text(hjust = 0.5, size = 16),
    #         text = element_text(size = 12),
    #         axis.text.x  = element_text(colour = "#000000", size = 12),
    #         axis.text.y  = element_blank(),
    #         axis.line.y.right = element_blank(),  # Removes the y-axis line on the right
    #         axis.ticks.y = element_blank(),
    #         axis.title.x = element_blank(),
    #         axis.title.y = element_blank(),
    #         strip.text.y.left = element_text(angle = 0, size = 16)) +
    #   ggtitle('Cerebellum') +
    #   facet_wrap(~feature,  ncol = 1, strip.position = "left") +
    #   scale_y_continuous(position = "right", limits=c(-0.00004, 5), breaks = 4)
    
    mid_plot <- create_stacked_vln_plot(adult_object, 
                                        'cellIDs', 
                                        dystonia_genes,
                                        adult_title, 
                                        get(paste0(region, '_vln_cols_recode')))
    
    saveRDS(mid_plot, paste0(R_dir, region, '_vln_plot.rds'))}
  
  if (region == 'str') {
    
    message('Plotting ', region, ' plot ...')
  
    # rhs_plot <- VlnPlot(adult_object, dystonia_genes, stack = TRUE, flip = TRUE,  
    #                     cols = plot_cols, same.y.lims = TRUE, fill.by = 'ident') +
    #   theme(legend.position = "none",
    #         plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
    #         panel.grid.major = element_blank(), 
    #         panel.grid.minor = element_blank(),
    #         plot.title = element_text(hjust = 0.5, size = 16),
    #         text = element_text(size = 12),
    #         axis.text.x  = element_text(colour = "#000000", size = 12),
    #         axis.text.y  = element_blank(),
    #         axis.line.y.right = element_blank(),  # Removes the y-axis line on the right
    #         axis.ticks.y = element_blank(),
    #         axis.title.x = element_blank(),
    #         axis.title.y = element_blank(),
    #         strip.text.y.left = element_text(angle = 0, size = 16)) +
    #   ggtitle('Striatum') +
    #   facet_wrap(~feature,  ncol = 1, strip.position = "left") +
    #   scale_y_continuous(position = "right", limits=c(-0.00004, 5), breaks = 4)
    
    rhs_plot <- create_stacked_vln_plot(adult_object, 
                                        'cellIDs', 
                                        dystonia_genes,
                                        adult_title, 
                                        get(paste0(region, '_vln_cols_recode')))
    
    saveRDS(rhs_plot, paste0(R_dir, region, '_vln_plot.rds'))}
  
}

if (Sys.info()[["nodename"]] == "Darrens-iMac-2.local") {
  
  # manually created this 

  agg_tbl <- read_excel(paste0(R_dir, 'aggr_exp_all_regions_norm.xlsx'), sheet = 'adult_aggr_fcx')
  sum_tbl <- agg_tbl |>
    rowwise() |>
    mutate(ExN_log10 = log10(sum(c_across(contains('ExN'))))) |>
    mutate(InN_log10 = log10(sum(c_across(contains('InN'))))) |>
    mutate(Oligo_log10 = log10(sum(c_across(contains('Olig'))))) |>
    mutate(Astro_log10 = log10(sum(c_across(contains('Ast'))))) |>
    mutate(Micro_log10 = log10(sum(c_across(contains('MG'))))) |>
    ungroup() |>
    select(gene, ExN_log10, InN_log10, Oligo_log10, Astro_log10, Micro_log10) |>
    mutate(across(everything(), ~ ifelse(. == -Inf, 0, .)))
    
  
  
  # Generate bar charts with facet wrap for each cell type
  ggplot(sum_tbl, aes(x = ExN_log10, y = gene, fill = '#00B6EB')) +
    geom_segment(aes(x = ExN_log10, xend = ExN_log10, y = as.numeric(gene) - 0.2, yend = as.numeric(gene) + 0.2), 
                 color = "#1f77b4", size = 1) +
    facet_wrap(~cell_type, ncol = 1, strip.position = "left") +
    scale_fill_manual(values = cell_type_colors) +
    labs(title = "Fetal Cerebellum", x = "Log10 Aggregated Expression", y = "") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 16),
      text = element_text(size = 12),
      axis.text.x = element_text(colour = "#000000", size = 12),
      axis.text.y = element_text(colour = "#000000", size = 12),
      axis.title.x = element_text(size = 14),
      strip.text.y.left = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0), "cm")
    ) +
    scale_x_continuous(position = "top", limits = c(-0.00004, 20), breaks = seq(0, 20, by = 5))
  
  
  for (cell_abbr in c('ExN', 'InN', 'Oligo', 'Astro', 'Micro')) {
    
    x_var <- paste0(cell_abbr, '_log10')
    title <- case_when(
      cell_abbr == 'ExN' ~ 'Excitatory Neurons',
      cell_abbr == 'InN' ~ 'Inhibitory Neurons',
      cell_abbr == 'Oligo' ~ 'Oligodendrocytes',
      cell_abbr == 'Astro' ~ 'Astrocytes',
      cell_abbr == 'Micro' ~ 'Microglia',
      .default = as.character(cell_abbr)
    )
    color <- case_when(
        cell_abbr == 'ExN' ~ "#1f77b4",
        cell_abbr == 'InN' ~ '#3CBB75FF',
        cell_abbr == 'Oligo' ~ '#FDE725FF',
        cell_abbr == 'Astro' ~ '#FF5959',
        cell_abbr == 'Micro' ~ '#F58231',
        .default = as.character(cell_abbr)
      )
    
    if (cell_abbr %in% c('ExN', 'InN')) y_size <- 11 else y_size <- 10
    
    plt <- ggplot(sum_tbl, aes(x = .data[[x_var]], y = fct_reorder(gene, .data[[x_var]]))) +
      geom_point(fill = color, size = 2, shape = 23) +
      labs(x = expression("Log"[10]*" aggregated expression"), 
           y = "") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 16),
            axis.text.y = element_text(size = y_size),
            axis.text.x = element_text(colour = "#000000", size = 11),
            axis.title.x = element_text(size = 12, vjust = -1)) +
      scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 2)) +
      ggtitle(title)
    
    assign(paste0(cell_abbr, '_plt'), plt)
  
  }
  
  nrns_plt <- cowplot::plot_grid(ExN_plt, InN_plt, ncol = 1, labels = c('B', 'C'), label_size = 20, scale = c(0.9, 0.9))
  glia_plt <- cowplot::plot_grid(Oligo_plt, Astro_plt, Micro_plt, ncol = 1, 
                                 labels = c('D', 'E', 'F'), label_size = 20, scale = c(0.95, 0.95, 0.95))
  cell_plt <- cowplot::plot_grid(nrns_plt, glia_plt, ncol = 2)
  comb_plt <- cowplot::plot_grid(lhs_plot, cell_plt, ncol = 2, 
                                 labels = c('A', ''), label_size = 20, rel_widths = c(2,3))

}
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
