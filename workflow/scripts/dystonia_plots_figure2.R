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
data <- read_delim(paste0(bulk_dir, 'FUMA_A/gtex_v8_ts_avg_log2TPM_exp.txt')) |>
  rename_with(~ str_replace_all(., "_", "-")) |>
  mutate(Gene = case_when(
    symbol %in% c("DNAJC12", "GNAO1", "ATP1A3", "HPCA", "TUBB4A") ~ "Set A",   # Assign Set A to gene X
    symbol %in% c("ADCY5", "PRRT2", "KCTD17") ~ "Set B",   # Assign Set B to gene Y
    symbol %in% c("SLC2A1", "SPR", "SGCE", "PNKD", "PRKRA", "VPS16") ~ "Set C",   # Assign Set C to gene Z
    TRUE ~ "Other gene"                  # Default label for other genes
  )) 
  
# Example matrix; replace this with your actual data
data_matrix <- as.matrix(data[, 3:56]) 
data_matrix <- apply(data_matrix, 2, as.numeric) 
rownames(data_matrix) <- data$symbol 

# Create tissue annotations using dplyr
tissue_annotations <- colnames(data_matrix) %>%
  as_tibble() %>%
  dplyr::rename(col_name = value) %>%  # Rename to something meaningful
  mutate(Tissue = case_when(
    col_name %in% c("Brain-Frontal-Cortex-BA9", "Brain-Cortex", "Brain-Anterior-cingulate-cortex-BA24") ~ "Tissue 1",   # Assign Set A to gene X
    col_name %in% c("Brain-Cerebellum", "Brain-Cerebellar-Hemisphere") ~ "Tissue 2",   # Assign Set B to gene Y
    col_name %in% c("Brain-Caudate-basal-ganglia", "Brain-Putamen-basal-ganglia", "Brain-Nucleus-accumbens-basal-ganglia") ~ "Tissue 3",
    TRUE ~ "Other Tissue"  # Label for any other tissue
  )) %>%
  select(-col_name) 

# Create a data frame for column annotations
annotation_col <- as.data.frame(tissue_annotations)
rownames(annotation_col) <- colnames(data_matrix)

col_range <- 1:ncol(data_matrix)  # Use all columns (or a specific subset of columns)
row_range <- 1:nrow(data_matrix) 

# Create heatmap with both row and column annotations
ht <- ComplexHeatmap::pheatmap(data_matrix, 
               row_dend_side = "right",       # Place dendrogram on the right side
               row_names_side = "left",
               annotation_row = data.frame(Gene = data$Gene),  # Add annotation for gene sets
               annotation_col = annotation_col,  # Add column annotations
               annotation_colors = list(Gene = c("Set A" = "blue", "Set B" = "green", "Set C" = "red", "Other gene" = "lightgrey"),
                                        Tissue = c("Tissue 1" = "purple", "Tissue 2" = "orange", "Tissue 3" = "pink", "Other Tissue" = "lightgrey")))  # Set colors for each tissue type

ComplexHeatmap::draw(ht) 
heatmap_name <- ht@name  # Retrieve the internal name used by ComplexHeatmap 

# Use the internal name in decorate_heatmap_body and adjust coordinates
decorate_heatmap_body(heatmap_name, { 
  
  grid.rect(x = unit((col_range[1] - 2.5) / ncol(data_matrix), "npc") + unit(0.5 / ncol(data_matrix), "npc"),  # X start 
            y = unit(1 - (row_range[1] - 1.5) / nrow(data_matrix), "npc") - unit(0.5 / nrow(data_matrix), "npc"), # Y start 
            width = unit((diff(col_range) + 12) / ncol(data_matrix), "npc"),   # Adjust width 
            height = unit(diff(row_range) / nrow(data_matrix), "npc"),  # Adjust height 
            just = c("left", "top"), 
            gp = gpar(col = "darkgreen", fill = NA, lwd = 2))  # Transparent center with red border 
  
  grid.rect(x = unit((col_range[1] - 2.5) / ncol(data_matrix), "npc") + unit(0.5 / ncol(data_matrix), "npc"),  # X start 
            y = unit(1 - (row_range[1] + 3.5) / nrow(data_matrix), "npc") - unit(0.5 / nrow(data_matrix), "npc"), # Y start 
            width = unit((diff(col_range) + 12) / ncol(data_matrix), "npc"),   # Adjust width 
            height = unit((diff(row_range) - 2) / nrow(data_matrix), "npc"),  # Adjust height 
            just = c("left", "top"), 
            gp = gpar(col = "darkgreen", fill = NA, lwd = 2)) 
  
  grid.rect(x = unit((col_range[1] - 2.5) / ncol(data_matrix), "npc") + unit(0.5 / ncol(data_matrix), "npc"),  # X start 
            y = unit(1 - (row_range[1] + 6.5) / nrow(data_matrix), "npc") - unit(0.5 / nrow(data_matrix), "npc"), # Y start 
            width = unit((diff(col_range) + 12) / ncol(data_matrix), "npc"),   # Adjust width 
            height = unit(diff(row_range) / nrow(data_matrix), "npc"),  # Adjust height 
            just = c("left", "top"), 
            gp = gpar(col = "darkgreen", fill = NA, lwd = 2)) 
  
})

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------