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
go_tbl_a <- read_delim(paste0(bulk_dir, 'FUMA_A/GS.txt')) |>
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

go_tbl_b <- read_delim(paste0(bulk_dir, 'FUMA_B/GS.txt')) |>
  mutate(negLog10 = -log10(adjP),
         Proportion = N_overlap / N_genes,
         GeneSet = str_replace_all(GeneSet, "_", " "),
         GeneSet = str_to_title(GeneSet),
         GeneSet = str_remove(GeneSet, "Gocc "),
         Category = str_replace_all(Category, "_", ":"),
         Category = str_to_upper(Category),
         Category = str_remove(Category, "GO:")) |>
  filter(GeneSet != 'Kegg Folate Biosynthesis')


## Plot ------
for (i in c('a', 'b')) {
  
  go_tbl <- get(paste0('go_tbl_', i))
  
  data_long <- go_tbl %>%
    select(GeneSet, Proportion, negLog10, Category) %>%
    pivot_longer(cols = c(Proportion, negLog10), names_to = "Metric", values_to = "Value") %>%
    mutate(Value = ifelse(Metric == "Proportion", -Value, Value),
           Metric = ifelse(Metric == "negLog10", "-log[10]*P[adj]","Proportion"),
           Metric = factor(Metric, levels = c("Proportion", "-log[10]*P[adj]")),
           GeneSet = factor(GeneSet, levels = sort(unique(GeneSet), decreasing = TRUE))) 
  
  # # Create a new column with text labels for facets
  # data_long$Metric <- ifelse(data_long$Metric == "negLog10", 
  #                                  "-log[10]*P[adj]", 
  #                                  "Proportion")
  
  # Now plot with the new label column
  go_plot <- ggplot(data_long, aes(x = GeneSet, y = Value, fill = Metric)) + 
    facet_grid(
      rows = vars(Category), 
      cols = vars(Metric), 
      scales = 'free', 
      space = 'free_y', 
      switch = "x",
      labeller = label_parsed) + 
    geom_col(width = 0.85, color = 'black', position = "identity") + 
    coord_flip() +  # Consistent limits across all plots
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
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
      strip.text.y = element_text(size = 13),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text = element_text(size = 13),
      panel.spacing.x = unit(.6, "lines")
    ) + 
    scale_y_continuous(
      expand = expansion(mult = c(0.03, 0.03)),
      labels = function(x) signif(abs(x), 3)
    ) + 
    NoLegend() + 
    xlab("") + 
    ylab("") 

  assign(paste0('go_', i, '_plt'), go_plot)
  
}

plot_grid(go_a_plt, go_b_plt, rel_heights = c(1.8,1), ncol = 1, align = 'v',
          labels = 'AUTO', label_size = 20)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

# Step 1: Expand genes into separate rows
expanded_genes <- go_tbl_a %>%
  separate_rows(genes, sep = ":")

# Step 2: Select only relevant columns to avoid conflicts
# Keep GeneSet, Category, and genes only
gene_presence_tbl <- expanded_genes %>%
  select(GeneSet, Category, genes) %>%  # Only keep necessary columns
  mutate(present = 1) %>%               # Create a presence indicator
  pivot_wider(names_from = genes, values_from = present, values_fill = 0)  # Fill absence with 0

# Step 3: Pivot to long format for plotting, keeping GeneSet and Category
heatmap_data <- gene_presence_tbl %>%
  pivot_longer(cols = -c(GeneSet, Category), names_to = "Gene", values_to = "Presence") |>
  mutate(GeneSet = factor(GeneSet, levels = sort(unique(GeneSet), decreasing = TRUE)))

# Plot with ggplot
tile_a_plt <- ggplot(heatmap_data, aes(x = Gene, y = GeneSet, fill = factor(Presence))) +
  geom_tile(color = "black") +
  scale_fill_manual(values = c("0" = "white", "1" = "#FFB300")) +  # Set colors for presence/absence
  facet_grid(rows = vars(Category), scales = "free_y", space = "free") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.border = element_rect(colour = "black", size = 1, fill = NA),
    axis.text.x = element_text(colour = "#000000", size = 12, angle = 45, hjust = 1),
    axis.text.y  = element_text(colour = "#000000", size = 12), 
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    strip.text.y = element_text(size = 13),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text = element_text(size = 13)
  ) +
  labs(fill = "Presence") + 
  NoLegend()

plot_grid(go_a_plt, go_b_plt, rel_heights = c(1.8,1), ncol = 1, align = 'v',
          labels = 'AUTO', label_size = 20)

plot_grid(go_a_plt, tile_a_plt), rel_heights = c(1.8,1), ncol = 1, align = 'v',
          labels = 'AUTO', label_size = 20)