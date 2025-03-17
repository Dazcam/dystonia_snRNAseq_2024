#--------------------------------------------------------------------------------------
#
#     Dystonia - Run EWCE GSEA Tests on Stiletti Superclusters
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# Note that we need to reverse engineer the EWCE ctd object from the Stiletti 
# specificity scores to run the bootstrap_enrichment_test()

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
  source(paste0(root_dir, 'workflow/scripts/dystonia_gene_lookup.R'))
  
} else {
  
  source('scripts/dystonia_functions.R')
  source('scripts/dystonia_Renvs.R')
  source('scripts/dystonia_gene_lists.R')

  
}

library(preprocessCore)  # For quantile normalization
# Load data - 406,046 rows initially
set.seed(1234)
superclust_spec_scores <- read_tsv(paste0(specificity_dir, 'Siletti_Supercluster_expression_specificity_TDEP_label.tsv.gz'))
lookup_hg38 <- get_biomart_gene_lookup()


# Reshape to wide format and check row sums
specificity_check <- superclust_spec_scores %>%
  dplyr::select(ENSGID, Supercluster, specificity) %>%
  pivot_wider(names_from = Supercluster, 
              values_from = specificity, 
              values_fill = 0) %>%
  # Calculate row sums for each gene
  rowwise() %>%
  mutate(row_sum = sum(c_across(-ENSGID), na.rm = TRUE)) %>%
  ungroup() %>%
  # Summarize how consistent the sums are
  summarise(mean_row_sum = mean(row_sum, na.rm = TRUE),
            sd_row_sum = sd(row_sum, na.rm = TRUE),
            min_row_sum = min(row_sum, na.rm = TRUE),
            max_row_sum = max(row_sum, na.rm = TRUE),
            n_rows_not_1 = sum(abs(row_sum - 1) > 0.01))  # Tolerance of 0.01 for floating-point precision

# Print the summary
print(specificity_check)

quantile_check <- superclust_spec_scores %>%
  group_by(Supercluster) %>%
  summarise(q_values = list(quantile(specificity, probs = seq(0, 1, 0.25), na.rm = TRUE))) %>%
  unnest_wider(q_values, names_sep = "_") %>%
  print(n = Inf)

# Reshape specificity scores into a matrix - 404,256 rows after join and filters
# bootstrap_enrichment_test() requires gene symbols
specificity_matrix <- superclust_spec_scores %>%
  left_join(lookup_hg38, by = join_by(ENSGID == ensembl_gene_id)) %>%
  dplyr::select(external_gene_name, Supercluster, specificity) %>%
  drop_na() |>
  filter(!external_gene_name == '') |>
  distinct(external_gene_name, Supercluster, .keep_all = TRUE) |>
  pivot_wider(names_from = Supercluster, 
              values_from = specificity, 
              values_fill = 0) %>%
  column_to_rownames(var = "external_gene_name") %>%
  as.matrix()


# Create pseudo-expression matrix
mean_exp_matrix <- specificity_matrix

# Compute specificity quantiles for subset
compute_quantiles_with_report <- function(matrix, numberOfBins) {
  specificity_quantiles <- matrix(NA, nrow = nrow(matrix), ncol = ncol(matrix))
  colnames(specificity_quantiles) <- colnames(matrix)
  rownames(specificity_quantiles) <- rownames(matrix)
  
  for (col in 1:ncol(matrix)) {
    supercluster <- colnames(matrix)[col]
    x <- matrix[, col]
    q <- quantile(x, probs = seq(0, 1, length.out = numberOfBins + 1), na.rm = TRUE)
    unique_breaks <- unique(q)
    cat(sprintf("Supercluster: %s\n", supercluster))
    cat(sprintf("  Unique values in specificity: %d\n", length(unique(x))))
    cat(sprintf("  Unique quantile breaks: %d (out of %d expected)\n", length(unique_breaks), numberOfBins + 1))
    if (length(unique_breaks) <= 1) {
      cat("  Warning: Too few unique breaks, assigning all to bin 1\n")
      specificity_quantiles[, col] <- rep(1, length(x))
    } else {
      quantiles <- cut(x, 
                       breaks = unique_breaks, 
                       include.lowest = TRUE, 
                       labels = FALSE)
      specificity_quantiles[, col] <- ifelse(is.na(quantiles), 1, quantiles)
    }
    cat(sprintf("  NAs in quantiles: %d\n", sum(is.na(specificity_quantiles[, col]))))
    cat("\n")
  }
  return(specificity_quantiles)
}

numberOfBins <- 40
specificity_quantiles <- compute_quantiles_with_report(specificity_matrix, numberOfBins)

# Load the template CTD
ctd <- ewceData::ctd()

# Update Level 1 with subset data
ctd[[1]]$mean_exp <- mean_exp_matrix
ctd[[1]]$specificity <- specificity_matrix
ctd[[1]]$specificity_quantiles <- specificity_quantiles
ctd[[1]]$median_exp <- mean_exp_matrix # Placeholder
ctd[[1]]$median_specificity <- specificity_matrix  # Placeholder

# Assume annot is a vector of cell type labels for consistency with 31 superclusters
n_superclusters <- ncol(specificity_matrix)  # 31
# For simplicity, repeat each supercluster name an arbitrary number of times (e.g., 100) to exceed matrix columns
ctd[[1]]$annot <- rep(colnames(specificity_matrix), each = 100)[1:3100]  # 31 * 100 = 3100, close to 3005

# Update cell_ordering
ctd[[1]]$cell_ordering <- colnames(specificity_matrix)

# Preserve plotting
if (!"plotting" %in% names(ctd[[1]])) {
  ctd[[1]]$plotting <- list()
}

# Define gene set and background (use Ensembl IDs for now)
dystonia_genes_44 <- dystonia_genes_44_hg38$hgnc_hg38  # Replace with your gene set
bg <- rownames(specificity_matrix)

# Run the bootstrap enrichment test
results <- bootstrap_enrichment_test(
  sct_data = ctd,
  hits = dystonia_genes_44,
  bg = bg,
  genelistSpecies = "human",
  sctSpecies = "human",
  reps = 10000,
  no_cores = 7,
  annotLevel = 1,
  store_gene_data = FALSE
)

knitr::kable(results$results)

try({
  plot_list <- EWCE::ewce_plot(total_res = results$results,
                               mtc_method ="BH",
                               ctd = ctd) 
  # print(plot_list$plain)
})


#### Does quantile normalisation have an effect?  ####
specificity_matrix_qn <- normalize.quantiles(specificity_matrix)
rownames(specificity_matrix_qn) <- rownames(specificity_matrix)
colnames(specificity_matrix_qn) <- colnames(specificity_matrix)
mean_exp_matrix_qn <- specificity_matrix_qn
specificity_quantiles_qn <- compute_quantiles_with_report(specificity_matrix_qn, 40)

ctd_qn <- ewceData::ctd()
ctd_qn[[1]]$mean_exp <- mean_exp_matrix_qn
ctd_qn[[1]]$specificity <- specificity_matrix_qn
ctd_qn[[1]]$specificity_quantiles <- specificity_quantiles_qn
ctd_qn[[1]]$annot <- rep(colnames(specificity_matrix_qn), each = 100)[1:3100]
ctd_qn[[1]]$cell_ordering <- colnames(specificity_matrix_qn)
if (!"plotting" %in% names(ctd_qn[[1]])) ctd_qn[[1]]$plotting <- list()

results_qn <- bootstrap_enrichment_test(
  sct_data = ctd_qn,
  hits = dystonia_genes_44,
  bg = bg,
  genelistSpecies = "human",
  sctSpecies = "human",
  reps = 10000,
  no_cores = 7,
  annotLevel = 1,
  store_gene_data = FALSE
)


#### Does gene size have an effect?  #####
results_qn_size <- bootstrap_enrichment_test(
  sct_data = ctd_qn,
  hits = dystonia_genes_44,
  bg = bg,
  genelistSpecies = "human",
  sctSpecies = "human",
  reps = 10000,
  no_cores = 7,
  annotLevel = 1,
  geneSizeControl = TRUE,
  store_gene_data = FALSE
)

#### Conditional normal ####
results_cond_MSN_orig <- bootstrap_enrichment_test(
  sct_data = ctd,
  hits = dystonia_genes_44,
  bg = bg,
  genelistSpecies = "human",
  sctSpecies = "human",
  reps = 10000,
  no_cores = 7,
  annotLevel = 1,
  #  geneSizeControl = TRUE,
  store_gene_data = FALSE,
  controlledCT = 'Medium_spiny_neuron'
)

results_cond_Splat_orig <- bootstrap_enrichment_test(
  sct_data = ctd,
  hits = dystonia_genes_44,
  bg = bg,
  genelistSpecies = "human",
  sctSpecies = "human",
  reps = 10000,
  no_cores = 7,
  annotLevel = 1,
  #  geneSizeControl = TRUE,
  store_gene_data = FALSE,
  controlledCT = 'Splatter'
)

results_cond_EMSN_orig <- bootstrap_enrichment_test(
  sct_data = ctd_qn,
  hits = dystonia_genes_44,
  bg = bg,
  genelistSpecies = "human",
  sctSpecies = "human",
  reps = 10000,
  no_cores = 7,
  annotLevel = 1,
  #  geneSizeControl = TRUE,
  store_gene_data = FALSE,
  controlledCT = 'Eccentric_medium_spiny_neuron'
)

#### Conditional qn ####
results_cond_MSN_qn <- bootstrap_enrichment_test(
  sct_data = ctd_qn,
  hits = dystonia_genes_44,
  bg = bg,
  genelistSpecies = "human",
  sctSpecies = "human",
  reps = 10000,
  no_cores = 7,
  annotLevel = 1,
#  geneSizeControl = TRUE,
  store_gene_data = FALSE,
  controlledCT = 'Medium_spiny_neuron'
)

results_cond_Splat_qn <- bootstrap_enrichment_test(
  sct_data = ctd_qn,
  hits = dystonia_genes_44,
  bg = bg,
  genelistSpecies = "human",
  sctSpecies = "human",
  reps = 10000,
  no_cores = 7,
  annotLevel = 1,
  #  geneSizeControl = TRUE,
  store_gene_data = FALSE,
  controlledCT = 'Splatter'
)

results_cond_EMSN_qn <- bootstrap_enrichment_test(
  sct_data = ctd_qn,
  hits = dystonia_genes_44,
  bg = bg,
  genelistSpecies = "human",
  sctSpecies = "human",
  reps = 10000,
  no_cores = 7,
  annotLevel = 1,
  #  geneSizeControl = TRUE,
  store_gene_data = FALSE,
  controlledCT = 'Eccentric_medium_spiny_neuron'
)

#### Reporting  ####
head(results$results)
head(results_cond_MSN_orig$results)
head(results_cond_EMSN_orig$results)
head(results_cond_Splat_orig$results)

head(results_qn$results)
head(results_cond_MSN_qn$results)
head(results_cond_EMSN_qn$results)
head(results_cond_Splat_qn$results)

p1 <- EWCE::ewce_plot(total_res = results$results, mtc_method ="BH", ctd = ctd)$plain + ggtitle('specificity')
p2 <- EWCE::ewce_plot(total_res = results_cond_MSN_orig$results, mtc_method ="BH", ctd = ctd)$plain + ggtitle('specificity condition on MSN')
p3 <- EWCE::ewce_plot(total_res = results_cond_EMSN_orig$results, mtc_method ="BH", ctd = ctd)$plain + ggtitle('specificity condition on EMSN')
p4 <- EWCE::ewce_plot(total_res = results_cond_Splat_orig$results, mtc_method ="BH", ctd = ctd)$plain + ggtitle('specificity condition on Splatter')
cowplot::plot_grid(p1, p2, p3, p4, ncol = 2)

p5 <- EWCE::ewce_plot(total_res = results_qn$results, mtc_method ="BH", ctd = ctd)$plain + ggtitle('quantile normalised')
p6 <- EWCE::ewce_plot(total_res = results_cond_MSN_qn$results, mtc_method ="BH", ctd = ctd)$plain + ggtitle('quantile normalised condition on MSN')
p7 <- EWCE::ewce_plot(total_res = results_cond_EMSN_qn$results, mtc_method ="BH", ctd = ctd)$plain + ggtitle('quantile normalised condition on EMSN')
p8 <- EWCE::ewce_plot(total_res = results_cond_Splat_qn$results, mtc_method ="BH", ctd = ctd)$plain + ggtitle('quantile normalised condition on Splat')
cowplot::plot_grid(p5, p6, p7, p8, ncol = 2)

plt <- EWCE::plot_ctd(ctd = ctd,
                      level = 1,
                      genes = dystonia_genes_44,
                      metric = "specificity") +
  theme(strip.text = element_text(size = 4))

####

library(EWCE)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

plot_specificity_heatmap <- function(ctd, genes_of_interest, cell_types = NULL) {
  
  # Extract specificity matrix from ctd object (assuming it's in level 1)
  spec_matrix <- ctd[[1]]$specificity
  
  # If specific cell types aren't specified, use all
  if(is.null(cell_types)) {
    cell_types <- colnames(spec_matrix)
  }
  
  # Subset the specificity matrix
  spec_subset <- spec_matrix[rownames(spec_matrix) %in% genes_of_interest, 
                             colnames(spec_matrix) %in% cell_types, 
                             drop = FALSE]
  
  # Remove rows with all zeros or NA
  spec_subset <- spec_subset[rowSums(spec_subset, na.rm = TRUE) > 0, , drop = FALSE]
  
  # Get BuPu colors (excluding the lightest ones for better contrast)
  #bupu_colors <- brewer.pal(9, "BuPu")[3:9]  # Start from 3rd color for darker tones
  # col_fun = colorRamp2(seq(min(spec_subset), max(spec_subset), length.out = 8), 
  #                      c("white", bupu_colors))
  col_fun = colorRamp2(c(0, 0.1, 0.3, max(spec_subset)), 
                       c("white", "lightblue", "purple", "black"))
  

  # # Create the heatmap
  ht = Heatmap(spec_subset,
          name = "Specificity",
          col = col_fun,
          cluster_rows = TRUE,
          cluster_columns = TRUE,
          clustering_distance_rows = "euclidean",
          clustering_distance_columns = "euclidean",
          clustering_method_rows = "complete",
          clustering_method_columns = "complete",
          row_names_side = "left",
          row_dend_side = "right",
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10),
          row_dend_width = unit(10, "mm"),
          column_dend_height = unit(10, "mm"),
          rect_gp = gpar(col = "black", lwd = 1),
          heatmap_legend_param = list(
            title = "Specificity",
            direction = "vertical"
          ),
          width = ncol(spec_subset) * unit(4, "mm"),
          height = nrow(spec_subset) * unit(4, "mm"))

  
  # Draw the heatmap
  draw(ht, heatmap_legend_side = "right")
}

plot_specificity_heatmap(ctd, dystonia_genes_44)
spec_matrix[rownames(spec_matrix) %in% genes_of_interest, , drop = FALSE]
test <- ctd[[1]]$specificity[rownames(ctd[[1]]$specificity) %in% dystonia_genes_44, ,drop = FALSE]
bupu_colors <- brewer.pal(9, "BuPu")[3:9]

colorRamp2(seq(min(ctd[[1]]$specificity), max(ctd[[1]]$specificity), length.out = 11), 
           c("white", bupu_colors))