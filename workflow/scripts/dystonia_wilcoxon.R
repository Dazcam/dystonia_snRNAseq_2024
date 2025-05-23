#--------------------------------------------------------------------------------------
#
#     Dystonia - Run Wilcoxon Tests on Stiletti Superclusters
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# Top 10% defined in TDEP: Note that the values are the top 10% of the cell type not the dataset 

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

source(paste0(script_dir, 'dystonia_gene_lookup.R'))

superclust_spec_scores <- read_tsv(paste0(specificity_dir, 'Siletti_Supercluster_expression_specificity_TDEP_label.tsv.gz'))

# All genes in all clusters have a specificity score > 0
superclust_spec_scores |>
  filter(specificity == 0)

# Top 
superclust_spec_scores |>
  group_by(Supercluster) |>
  dplyr::count()

superclust_spec_scores |>
  group_by(Supercluster, TDEP) |>
  dplyr::count() |>
  print(n = Inf)


# Initialize an empty tibble
superclust_wlcxn_tbl <- tibble()
p_values <- c()

for (cell_pop in unique(superclust_spec_scores$Supercluster)) {
  
  missing_genes_hgnc <- ""
  
  message('\nRunning Wilcoxon Rank sum test for: ', cell_pop) 
  
  # Filter data for the current Supercluster and determine gene presence
  cell_type_scores <- superclust_spec_scores |>
    filter(Supercluster == cell_pop) |>
    mutate(dyst_gene = ifelse(ENSGID %in% dystonia_genes_44_tbl$ensembl_hg38, 1, 0))
  
  # Check for unexpressed genes
  message('Number for genes with specificty values = 0: ', sum(cell_type_scores$specificity == 0))
  
  genes_exp_in_dataset <- cell_type_scores |>
    filter(ENSGID %in% dystonia_genes_44_tbl$ensembl_hg38) |>
    pull(ENSGID) |>
    unique()
  
  message(length(genes_exp_in_dataset), ' dystonia genes are expressed each in cell type')
  #message(sum(cell_type_scores$ENSGID %in% KMT2B_hg19_encoding), ' extra gene in cell type')
  
  # If any genes are missing print gene names
  if (length(genes_exp_in_dataset) != 44) {
    
    message('Missing genes:')
    missing_genes <- setdiff(dystonia_genes_44_tbl$ensembl_hg38, genes_exp_in_dataset)
    missing_genes_hgnc <- lookup_hg38 |> 
      filter(ensembl_gene_id %in% missing_genes) |>
      arrange(hgnc_symbol) |>
      print() |>
      pull(hgnc_symbol) |>
      paste(collapse = ';')}
  
  # Extract specificity scores for dystonia genes and non-dystonia genes
  cell_type_scores <- superclust_spec_scores |>  
    filter(Supercluster == cell_pop) |>  
    mutate(dyst_gene = ifelse(ENSGID %in% dystonia_genes_44_tbl$ensembl_hg38, 1, 0))  
  
  non_dyst_gene_spec_scores  <- cell_type_scores %>%  
    filter(dyst_gene == 0) %>%  
    pull(specificity)  
  
  dyst_gene_spec_scores <- cell_type_scores %>%  
    filter(dyst_gene == 1) %>%  
    pull(specificity)  
  
  message('Number of non-dyst gene scores: ', length(non_dyst_gene_spec_scores))
  message('Number of dyst gene scores: ', length(dyst_gene_spec_scores))
  
  # Ensure both groups have values before running Wilcoxon test  
  if (length(dyst_gene_spec_scores) > 0 & length(non_dyst_gene_spec_scores) > 0) {  
    wilcox_result <- wilcox.test(dyst_gene_spec_scores, non_dyst_gene_spec_scores,  
                                 alternative = 'greater', paired = FALSE)  
    
    p_value <- wilcox_result$p.value  
    p_values <- c(p_values, p_value)  # Store for correction  
    
    # Store results in a tibble  
    superclust_wlcxn_tbl <- superclust_wlcxn_tbl |>  
      bind_rows(tibble(Cluster = cell_pop, 
                       n_genes = length(dyst_gene_spec_scores), 
                       n_genes_bckgrnd = length(non_dyst_gene_spec_scores),  
                       p_value = p_value,
                       missing_genes = missing_genes_hgnc))  
  } else {  
    message("Skipping test due to insufficient data")  
  }  
  
  message('P-value: ', wilcox_result$p.value)
  
}

# Apply multiple testing corrections
superclust_wlcxn_tbl <- superclust_wlcxn_tbl |>  
  mutate(BF = p.adjust(p_value, method = "bonferroni", n = 31),  
         FDR = p.adjust(p_value, method = "fdr")) |> 
  write_tsv(paste0(results_dir, '05tables/dystonia_wilcoxon_results.tsv'))
  

# View results
print(superclust_wlcxn_tbl, n = Inf)

test <- superclust_wlcxn_tbl |>
  mutate(Cluster = str_replace_all(Cluster, '_', ' ')) |>
  mutate(Cluster = factor(Cluster, 
                            levels = c("Vascular", "Upper rhombic lip", "Upper layer intratelencephalic", 
                                       "Thalamic excitatory", "Splatter", "Oligodendrocyte precursor", 
                                       "Oligodendrocyte", "Miscellaneous", "Midbrain derived inhibitory", 
                                       "Microglia", "Medium spiny neuron", "Mammillary body", "MGE interneuron", 
                                       "Lower rhombic lip", "LAMP5 LHX6 and Chandelier", "Hippocampal dentate gyrus", 
                                       "Hippocampal CA4", "Hippocampal CA1 3", "Fibroblast", "Ependymal", 
                                       "Eccentric medium spiny neuron", "Deep layer near projecting", 
                                       "Deep layer intratelencephalic", "Deep layer corticothalamic and 6b", 
                                       "Committed oligodendrocyte precursor", "Choroid plexus", "Cerebellar inhibitory", 
                                       "CGE interneuron", "Bergmann glia", "Astrocyte", "Amygdala excitatory"
                            ))) |>
  mutate(P_fdr_bool = ifelse(FDR < 0.05, TRUE, FALSE)) |>
  mutate(neglog10p = -log10(p_value)) |>
  ggplot(aes(x = Cluster, y = neglog10p)) +
  geom_col(aes(fill = P_fdr_bool), color = 'black') +
  coord_flip() +  # Flip the coordinates for better readability
  labs(x = "", y = '-log10P') +
  
  theme_minimal() +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 13),
        axis.text.x =  element_text(size = 15, colour = 'black'),
        axis.text.y =  element_text(size = 15, colour = 'black'),
        axis.title.x = element_text(size = 15, vjust = -1),
  ) +
  Seurat::NoLegend() +
  xlab("") + 
  ylab(expression('-log'[10]*'(P)')) +
  ylim(0, 3)


# Compute Spearman and Pearson correlations
# cor_spearman <- cor.test(superclust_wlcxn_tbl$n_genes, superclust_wlcxn_tbl$FDR, method = "spearman")
# cor_pearson <- cor.test(superclust_wlcxn_tbl$n_genes, superclust_wlcxn_tbl$FDR, method = "pearson")
# 
# # Print results
# print(cor_spearman)
# print(cor_pearson)
# 
# # Scatter plot with regression line
# ggplot(superclust_wlcxn_tbl, aes(x = n_genes, y = FDR)) +
#   geom_point() +
#   geom_smooth(method = "lm", color = "blue", se = TRUE) +
#   labs(title = "Correlation between Number of Genes and FDR",
#        x = "Number of Genes",
#        y = "FDR") +
#   theme_minimal()
# 
# dystonia_genes_44_tbl |>
#   dplyr::select(hgnc_hg38, ensembl_hg38) |>
#   filter(!hgnc_hg38 %in% c('GCH1', 'TH', 'SPR', 'DDC', 'SLC6A3')) |>
#   pull(ensembl_hg38) |>
#   cat(sep = '\n')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------