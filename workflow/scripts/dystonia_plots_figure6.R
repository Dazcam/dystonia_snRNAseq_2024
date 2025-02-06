#--------------------------------------------------------------------------------------
#
#    Dystonia - Create figure 6
#
#--------------------------------------------------------------------------------------

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
go_tbl_a <- read_delim(paste0(bulk_dir, 'FUMA_gene2func555089_111124/GS.txt')) |>
  mutate(negLog10 = -log10(adjP),
         Proportion = N_overlap / N_genes,
         GeneSet = str_replace_all(GeneSet, "_", " "),
         GeneSet = str_to_title(GeneSet),
         GeneSet = str_remove(GeneSet, "Gobp "),
         GeneSet = str_remove(GeneSet, "Gocc "),
         Category = str_replace_all(Category, "_", ":"),
         Category = str_to_upper(Category),
         Category = str_replace_all(Category, "CELL:TYPE:SIGNATURE", "CT-SIG"),
         Category = str_remove(Category, "GO:")) |>
  filter(GeneSet != 'Kegg Folate Biosynthesis')

go_tbl_b <- read_delim(paste0(bulk_dir, 'FUMA_gene2func555124_no_SPR_GCH1_TH_111124/GS.txt')) |>
  mutate(negLog10 = -log10(adjP),
         Proportion = N_overlap / N_genes,
         GeneSet = str_replace_all(GeneSet, "_", " "),
         GeneSet = str_to_title(GeneSet),
         GeneSet = str_remove(GeneSet, "Gocc "),
         Category = str_replace_all(Category, "_", ":"),
         Category = str_to_upper(Category),
         Category = str_remove(Category, "GO:"),
         GeneSet = str_replace(GeneSet, "Plasma Membrane Protein Complex", # Padding to help with plot alignment below
                               str_pad("Plasma Membrane Protein Complex", 
                                       width = 43, side = "left"))) |>
  filter(GeneSet != 'Kegg Folate Biosynthesis')

# Calculate the range of Proportion and -log[10]*P[adj] across both datasets
# common_proportion_range <- range(c(go_tbl_a$Proportion, go_tbl_b$Proportion), na.rm = TRUE)
# common_neglog_range <- range(c(go_tbl_a$negLog10, go_tbl_b$negLog10), na.rm = TRUE)


## Plot ------
# Prepare bar plots
for (i in c('a', 'b')) {
  
  go_tbl <- get(paste0('go_tbl_', i))
  
  data_long <- go_tbl %>%
    select(GeneSet, Proportion, negLog10, Category) %>%
    pivot_longer(cols = c(Proportion, negLog10), names_to = "Metric", values_to = "Value") %>%
    mutate(Value = ifelse(Metric == "Proportion", -Value, Value),
           Metric = ifelse(Metric == "negLog10", "-log[10]*P[adj]","Proportion"),
           Metric = factor(Metric, levels = c("Proportion", "-log[10]*P[adj]")),
           GeneSet = factor(GeneSet, levels = sort(unique(GeneSet), decreasing = TRUE))) 
  
  if (i == 'b') {

    # Need this as padding for b affects aplhabetical ordering
    data_long <- data_long %>%
      mutate(GeneSet = fct_relevel(GeneSet, levels(GeneSet)[11], after = 3))
  }
  
  # Now plot with the new label column
  go_plot <- ggplot(data_long, aes(x = GeneSet, y = Value, fill = Metric)) + 
    facet_grid(
      rows = vars(Category),  
      cols = vars(Metric),  
      scales = 'free',  # Allow separate scales for each facet
      space = 'free_y',  
      switch = "x",  
      labeller = label_parsed
    ) +  
    geom_col(width = 0.80, color = 'black', position = "identity") +  
    coord_flip() +  
    theme(
      plot.margin = unit(c(0.5, 0.11, 0.5, 0.5), "cm"),
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", size = 1, fill = NA),
      plot.title = element_text(hjust = 0.5, face = 'bold'),
      axis.title.x = element_text(colour = "#000000", size = 14),
      axis.title.y = element_text(colour = "#000000", size = 14),
      axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5),
      axis.text.y  = element_text(colour = "#000000", size = 12),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      strip.text.y = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text = element_text(size = 13),
      panel.spacing.y = unit(.38, "lines") # Different from b plot due to coord_flip
    ) + 
    ggh4x::facetted_pos_scales(
      y = list(
        `Proportion` = scale_y_continuous(
          limits = c(-0.2, 0),  # Set limits to ensure proper space for the bars
          breaks = c(0, -0.05, -0.10, -0.15),  # Define specific breaks for Proportion
          expand = expansion(mult = c(0.04, 0.04)), 
          labels = function(x) signif(abs(x), 3),  # Flip the axis values to show positive values
          
        ),
        `-log[10]*P[adj]` = scale_y_continuous(
          limits = c(0, 3.4),  
          breaks = c(0, 1, 2, 3), 
          expand = expansion(mult = c(0.04, 0.04)), 
          labels = function(x) signif(x, 3)
        )
      )
    ) + 
    NoLegend() + 
    xlab("") + 
    ylab("")  

  assign(paste0('go_', i, '_plt'), go_plot)
  
}

# Prepare geom tile plots
for (i in c('a', 'b')) {
  
  go_tbl <- get(paste0('go_tbl_', i))
  
  # Expand genes into separate rows
  expanded_genes <- go_tbl %>%
    separate_rows(genes, sep = ":")
  
  # Select only relevant cols to avoid conflicts
  # Keep GeneSet, Category, and genes only
  gene_presence_tbl <- expanded_genes %>%
    select(GeneSet, Category, genes) %>%  # Only keep necessary columns
    mutate(present = 1) %>%               # Create a presence indicator
    pivot_wider(names_from = genes, values_from = present, values_fill = 0)  # Fill absence with 0
  
  if (i == 'b') {
    
    # Add genes in A that are not in B for continuity
    gene_presence_tbl <- gene_presence_tbl |>
      mutate("ADCY5" = 0,
             "GCH1" = 0 ,
             "TH" = 0,
             "DNAJC12" = 0)
  }
  
  # Step 3: Pivot to long format for plotting, keeping GeneSet and Category
  heatmap_data <- gene_presence_tbl %>%
    pivot_longer(cols = -c(GeneSet, Category), names_to = "Gene", values_to = "Presence") |>
    mutate(GeneSet = factor(GeneSet, levels = sort(unique(GeneSet), decreasing = TRUE)))
  
  # Plot with ggplot
  tile_plt <- ggplot(heatmap_data, aes(x = Gene, y = GeneSet, fill = factor(Presence))) +
    geom_tile(color = "black") +
    scale_fill_manual(values = c("0" = "white", "1" = "#FFB300")) +  # Set colors for presence/absence
    facet_grid(rows = vars(Category), scales = "free_y", space = "free") +
    theme_minimal() +
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.11), "cm"),
      panel.border = element_rect(colour = "black", size = 1, fill = NA),
      axis.text.x = element_text(colour = "#000000", size = 12, angle = 45, hjust = 1),
      axis.text.y  = element_blank(), 
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.text.y = element_text(size = 13),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text = element_text(size = 13),
      axis.ticks.x = element_line(color = "black"),
      panel.spacing.x = unit(.6, "lines")
    ) +
    labs(fill = "Presence") + 
    scale_x_discrete(expand = c(0, 0)) +  # Ensure no space around x-axis
    scale_y_discrete(expand = c(0, 0)) +  # Ensure no space around y-axis
    NoLegend()
  
  assign(paste0('tile_', i, '_plt'), tile_plt)

}

join_a_plt <- plot_grid(go_a_plt, tile_a_plt, axis = 'tblr', align = 'h',
                        ncol = 2, labels = c('A', ''), label_size = 20)
join_b_plt <- plot_grid(go_b_plt, tile_b_plt, axis = 'tblr', align = 'h',
                        ncol = 2, labels = c('B', ''), label_size = 20)

plot_grid(join_a_plt, join_b_plt, rel_heights = c(1.8,1), axis = 'lr', ncol = 1, align = 'v',
          labels = 'AUTO', label_size = 20)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

tile_a_plt+
  scale_x_discrete(expand = expansion(mult = c(0.1, 0.1))) +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.1)))
