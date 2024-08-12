#--------------------------------------------------------------------------------------
#
#    Dystonia - hdWGCNA reporting
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# Vingette: https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html

# Useful issues:

# https://github.com/smorabit/hdWGCNA/issues/275
# https://github.com/smorabit/hdWGCNA/issues/27


##  Load Packages, functions and variables  -------------------------------------------
message('Setting environment variables ...')
if (Sys.info()[["nodename"]] == "Darrens-iMac-2.local") {
  
  library(yaml)
  root_dir <- '~/Desktop/dystonia_snRNAseq_2024/'
  yaml_file <- yaml.load_file(paste0(root_dir, 'config/config.yaml'))
  region <- yaml.load(yaml_file$region) # Note this is different from other scripts
  
  source(paste0(root_dir, 'workflow/scripts/dystonia_functions.R'))
  source(paste0(root_dir, 'workflow/scripts/dystonia_Renvs.R'))
  source(paste0(root_dir, 'workflow/scripts/dystonia_gene_lists.R'))
  
} else {
  
  source('scripts/dystonia_functions.R')
  source('scripts/dystonia_Renvs.R')
  source('scripts/dystonia_gene_lists.R')
  
}

# Read in cell types that hdWGCNA was run on
cell_types <- read_tsv(paste0(wgcna_dir, region, '_metacells.tsv'), col_names = 'cell_id') 

for (cell_type in cell_types$cell_id) {
  
  message('Generating hdWGCNA plot list for ', cell_type, '...\n')
  
  message('Read in Seurat object ...\n')
  seurat_obj <- readRDS(paste0(wgcna_dir, cell_type, '_hdWGCNA.rds'))
  
  # Moved to Rmd
  # message('Calculate metacell proprtions ...\n')
  # meta_obj <- GetMetacellObject(seurat_obj)
  # cell_props <- get_wgcna_cell_props(seurat_obj, meta_obj)
  
  message('Setting TOM file location ...')
  seurat_obj@misc[["str_wgcna"]][["wgcna_net"]][["consensusTOMInfo"]][["saveConsensusTOMs"]] <- paste0(wgcna_dir, cell_type, '_TOM.rda')
  seurat_obj@misc[["str_wgcna"]][["wgcna_net"]][["TOMFiles"]] <- paste0(wgcna_dir, cell_type, '_TOM.rda')
  
  # Run overlaps between dystonia genes and WGCNA modules - atm 50 hub genes
  message("Calculating overlap genes ...")
  overlap_genes_lst <- run_dyst_gene_overlap(seurat_obj, cell_type, region,
                                         paste0(region, '_wgcna'), wgcna_dir, 50)
  overlap_genes <- bind_rows(overlap_genes_lst)
  
  # Skip if no module gene overlaps are found
  if (nrow(overlap_genes) == 0) {
    message('WARNING! No dystonia gene overlaps for ', cell_type, ' \n')
    next
  }
  
  message("Getting stats tbl ...")
  wgcna_stats_tbl <- get_wgcna_stats(seurat_obj, cell_type, region,
                                     paste0(region, '_wgcna'), wgcna_dir,
                                     nrow(bind_rows(overlap_genes)))

  # Plot Dendogram
  message("Plotting Dendogram ...")
  PlotDendrogram(seurat_obj, main = paste0(cell_type, ' hdWGCNA Dendrogram'))
  
  dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')
  
  message("Running EnrichR for GO terms for each module ...")
  seurat_obj <- RunEnrichr(
    seurat_obj,
    dbs = dbs, # character vector of enrichr databases to test
    max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
  )
  
  # Retrieve the output table
  enrichr_tbl <- GetEnrichrTable(seurat_obj) |> 
    as_tibble() |>
    filter(Adjusted.P.value < 0.05) |>
    group_by(module)
  
  enrich_top20_lst <- enrichr_tbl |>
    slice_head(n = 20) |> 
    group_split(module) 
  names(enrich_top20_lst) <- group_keys(enrichr_tbl) |> pull()
  
  enrich_lst <- enrichr_tbl |>
    group_split(module) 
  names(enrich_lst) <- group_keys(enrichr_tbl) |> pull()
  
  # Calculate how many GO terms at P < 0.05 each module has
  go_terms <- vector()
  for (mod in unique(overlap_genes$module)) {
    
    tryCatch(
      {
        num_go_terms <- enrich_top20_lst[[mod]] |>
          filter(module == mod & Adjusted.P.value < 0.05) |>
          nrow()
        go_terms <- c(go_terms, num_go_terms)
      },
      error = function(err){
        message('WARNING! No GO terms for ', mod, ' at P < 0.05 \n')
        go_terms <<- c(go_terms, 0)
      }
    )
    
  }
  
  go_module_tbl <- overlap_genes |> 
    distinct(module) |>
    mutate('No. GO terms' = go_terms)
  
  # Plot only modules with dystonia gene overlaps
  module_overlaps <- overlap_genes |> 
    select(module) |> 
    filter(module %in% names(enrich_lst)) |> # Filt mods with overlaps P > 0.05
    unique() |> 
    pull()
  
  plot_list <- list()
  
  for (module in module_overlaps) {
    
    message("Number of modules with overalpping dystonia genes: ", length(module_overlaps), "...")
    message("Plotting GO plot for: ", module, "...")
    go_plot <- ggplot(data = enrich_top20_lst[[module]], 
                      aes(y = Term, 
                          x = module, 
                          color = -log10(Adjusted.P.value), 
                          size = Odds.Ratio)) +
      geom_point() +
      theme_bw() +
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            #      panel.grid.major = element_blank(), 
            #    panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", size = 1),
            plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(colour = "#000000", size = 17, vjust = -0.5),
            axis.title.y = element_text(colour = "#000000", size = 17),
            axis.text.x  = element_text(colour = "#000000", size = 15, vjust = 0.5),
            axis.text.y  = element_text(colour = "#000000", size = 15),
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 15)) +
      scale_color_gradient(low = "blue", high = "red") +
      theme(axis.text.x = element_text(colour = "#000000")) +
      ylab("Term") 
    #scale_radius(limits = c(1, 6), range = c(1,10))
    
    plot_list[[module]] <- go_plot
  
  }
  
  rmarkdown::render(markdown_wgcna_doc, output_file = paste0("dystonia_wgcna_", 
                                                             cell_type, ".html"), 
                    output_dir = wgcna_dir)
  rm(seurat_obj)
  
}

file.create(paste0(wgcna_dir, region, '_wgcna_plots.done'))

## Other plots  -----------------------------

# Plot average expression
# get hMEs from seurat object
# NAs in Idents after using sketch object throws error: https://github.com/satijalab/seurat/issues/8772
# Idents(seurat_obj) <- seurat_obj$cellIDs
# MEs <- GetMEs(seurat_obj, harmonized=TRUE)
# modules <- GetModules(seurat_obj)
# mods <- levels(modules$module); mods <- mods[mods != 'grey']
# seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)
# 
# module_av_exp_plt <- DotPlot(seurat_obj, features=mods, group.by = 'cellIDs') +
#   RotatedAxis() +
#   scale_color_gradient2(high='red', mid='grey95', low='blue')
# 
# # Module Network plot
# ModuleNetworkPlot(
#   seurat_obj,
#   outdir = paste0(wgcna_dir, 'ModuleNetworks')
# )
  
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------





