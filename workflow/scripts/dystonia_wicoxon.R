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

#source(paste0(script_dir, 'fetal_col_recode.R'))
source(paste0(script_dir, 'dystonia_functions.R'))

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
  
  message('\nRunning Wilcoxon Rank sum test for: ', cell_pop) 
  
  # Filter data for the current Supercluster and determine gene presence
  cell_type_scores <- superclust_spec_scores |>
    filter(Supercluster == cell_pop) |>
    mutate(chick_gene = ifelse(ENSGID %in% dystonia_genes_44_tbl$ensembl_hg38, 1, 0))
  
  # Check for unexpressed genes
  message('Number for genes with specificty values = 0: ', sum(cell_type_scores$specificity == 0))
  
  genes_exp_in_dataset <- cell_type_scores |>
    filter(ENSGID %in% dystonia_genes_44_tbl$ensembl_hg38) |>
    pull(ENSGID) |>
    unique()
  
  message(length(genes_exp_in_dataset), ' dystonia genes are expressed each in cell type')
  #message(sum(cell_type_scores$ENSGID %in% KMT2B_hg19_encoding), ' extra gene in cell type')
  
  # If any genes are missing print gene names
  if (length(genes_exp_in_dataset) != 8) {
    
    message('Missing genes:')
    missing_genes <- setdiff(dystonia_genes_44_tbl$ensembl_hg38, genes_exp_in_dataset)
    lookup_hg38 |> 
      filter(ensembl_gene_id %in% missing_genes) |> 
      print()}
  
  # Extract specificity scores for chick genes and non-chick genes
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
                       p_value = p_value))  
  } else {  
    message("Skipping test due to insufficient data")  
  }  
  
  message('P-value: ', wilcox_result$p.value)
  
}

# Apply multiple testing corrections
superclust_wlcxn_tbl <- superclust_wlcxn_tbl |>  
  mutate(BF = p.adjust(p_value, method = "bonferroni", n = 31),  
         FDR = p.adjust(p_value, method = "fdr"))  

# View results
print(superclust_wlcxn_tbl, n = Inf)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------