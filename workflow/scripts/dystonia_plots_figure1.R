#--------------------------------------------------------------------------------------
#
#    Dystonia - Create figure 2 - Heatmap
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
source(paste0(script_dir, 'dystonia_brain_span_analysis.R'))
# Note we used P values here and not adjusted Ps and manually applied BF correction for 54 tests
deg_a_tbl <- read_delim(paste0(bulk_dir, 'FUMA_gene2func555089_111124/gtex_v8_ts_DEG.txt')) |>
  mutate(negLog10 = -log10(p),
         GeneSet = str_replace_all(GeneSet, "_", " ")) |>
  group_by(Category) |>
  arrange(p) |>
  ungroup() 
  

deg_gen_tbl <- read_delim(paste0(bulk_dir, 'FUMA_gene2func555089_111124/gtex_v8_ts_general_DEG.txt')) |>
  mutate(negLog10 = -log10(p),
         GeneSet = str_replace_all(GeneSet, "_", " ")) |>
  group_by(Category) |>
  arrange(p) |>
  ungroup() 


for (i in c("a", 'gen')) {
  
  deg_tbl <- get(paste0("deg_", i,"_tbl"))
  
  gene_order <- deg_tbl %>%
    filter(Category == "DEG.up") %>%
    arrange(desc(negLog10)) %>%
    pull(GeneSet)

  data_long <- deg_tbl %>%
    select(GeneSet, negLog10, Category) %>%
    pivot_longer(cols = c(negLog10), names_to = "Metric", values_to = "negLog10") |>
    mutate(GeneSet = factor(GeneSet, levels = gene_order),
           Category = str_replace_all(Category, "\\.", " "),
           Category = factor(Category, levels = c("DEG up", "DEG twoside", "DEG down"),),
           sig = ifelse(negLog10 > (-log10(0.05 / 54)), TRUE, FALSE))

  deg_plot <-  ggplot(data_long, aes(x = GeneSet, y = negLog10, fill = sig)) + 
    facet_grid(
      rows = vars(Category), 
      cols = vars(Metric), 
      scales = 'free', 
      space = 'fixed') + 
    geom_col(width = 0.85, color = 'black', position = "identity") + 
    theme(
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", size = 1, fill = NA),
      plot.title = element_text(hjust = 0.5, face = 'bold'),
      axis.title.x = element_text(colour = "#000000", size = 14),
      axis.title.y = element_text(colour = "#000000", size = 14),
      axis.text.x  = element_text(colour = "#000000", size = 12, vjust = 0.5, angle = 90, hjust = 1),
      axis.text.y  = element_text(colour = "#000000", size = 12),
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      strip.text.y = element_text(size = 13),
      strip.text.x = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text = element_text(size = 13),
      panel.spacing.x = unit(.6, "lines")
    ) + 
    NoLegend() + 
    xlab("") + 
    ylab(expression('-log'[10]*'P'))
  
  assign(paste0('deg_', i, '_plt'), deg_plot)
  
}

a_plt <- plot_grid(dev_stage_plt,  deg_gen_plt, labels = 'AUTO', label_size = 20, scale = c(0.9,1))

plot_grid(a_plt,  deg_a_plt, ncol = 1, labels = c('', 'C'), label_size = 20, rel_heights = c(1,1.5))
