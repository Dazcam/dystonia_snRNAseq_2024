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
#fuma_dir_pre_review1 <- 'FUMA_gene2func555089_111124'
fuma_dir <- 'FUMA_44_genes_050225/'

#paste0(bulk_dir, 'FUMA_A/gtex_v8_ts_avg_log2TPM_exp.txt') # Old: no KMT2B
data <- read_delim(paste0(bulk_dir, fuma_dir, 'gtex_v8_ts_avg_log2TPM_exp.txt')) |>
  rename_with(~ str_replace_all(., "_", "-")) |>
  mutate(symbol = recode(symbol, "C9orf3" = "AOPEP")) |>
  mutate(symbol = recode(symbol, "BZRAP1" = "TSPOAP1")) |>
  arrange(symbol) |>
  mutate(Gene = case_when(
    symbol %in% c("ADCY5", "PRRT2", "KCTD17", "SLC2A1") ~ "Set A",
    symbol %in% c("SGCE", "PNKD", "PRKRA", "VPS16", "KMT2B") ~ "Set B", 
    symbol %in% c("DNAJC12", "GNAO1", "ATP1A3", "HPCA", "TUBB4A") ~ "Set C", 
    TRUE ~ ""                  # Default label for other genes
  )) 
  
# Example matrix; replace this with your actual data
data_matrix <- as.matrix(data[, 3:56]) 
data_matrix <- apply(data_matrix, 2, as.numeric) 
rownames(data_matrix) <- data$symbol 

# Create tissue annotations using dplyr
tissue_annotations <- colnames(data_matrix) %>%
  as_tibble() %>%
  dplyr::rename(col_name = value) %>%  # Rename to something meaningful
  mutate(Region = case_when(
    col_name %in% c("Brain-Frontal-Cortex-BA9", "Brain-Cortex", "Brain-Anterior-cingulate-cortex-BA24") ~ "Cortex", 
    col_name %in% c("Brain-Cerebellum", "Brain-Cerebellar-Hemisphere") ~ "Cerebellum",  
    col_name %in% c("Brain-Caudate-basal-ganglia", "Brain-Putamen-basal-ganglia", "Brain-Nucleus-accumbens-basal-ganglia") ~ "Striatum",
    col_name %in% c("Brain-Amygdala", "Brain-Cerebellar-Hemisphere") ~ "Amygdala",  
    col_name %in% c("Brain-Cerebellum", "Brain-Substantia-nigra") ~ "Substantia Nigra",  
    col_name %in% c("Brain-Cerebellum", "Brain-Hippocampus") ~ "Hippocampus",  
    TRUE ~ ""  # Label for any other tissue
  )) %>%
  dplyr::select(-col_name) 

# Create a data frame for column annotations
annotation_col <- as.data.frame(cbind(tissue_annotations))
rownames(annotation_col) <- colnames(data_matrix)

col_range <- 1:ncol(data_matrix)  # Use all columns (or a specific subset of columns)
row_range <- 1:nrow(data_matrix) 

colnames(data_matrix) <- str_replace_all(colnames(data_matrix), '-', ' ')

# Create heatmap with both row and column annotations
ht <- ComplexHeatmap::pheatmap(data_matrix, 
               row_dend_side = "right",       # Place dendrogram on the right side
               row_names_side = "left",
               name = 'Expression',
               clustering_method = "average",
               annotation_col = annotation_col,  # Add column annotations
               annotation_colors = list(Region = c("Cortex" = "purple", "Cerebellum" = "orange", "Striatum" = "pink", 
                                                   "Amygdala" = "yellow", "Substantia Nigra"  = "cyan" , "Hippocampus"  = "green")),
               border = FALSE,                
               fontsize_row = 12,             
               fontsize_col = 12,             
               border_color = "black",        
               fontsize = 12,
               cellwidth = 12,               
               cellheight = 12
)
               
ComplexHeatmap::draw(ht)
heatmap_name <- ht@name  # Retrieve the internal name used by ComplexHeatmap 
ht@column_title <- ""

# Use the internal name in decorate_heatmap_body and adjust coordinates
# decorate_heatmap_body(heatmap_name, { 
#   
#   grid.rect(x = unit((col_range[1] - 1.5) / ncol(data_matrix), "npc") + unit(0.5 / ncol(data_matrix), "npc"),  # X start 
#             y = unit(1 - (row_range[1] + 2.5) / nrow(data_matrix), "npc") - unit(0.5 / nrow(data_matrix), "npc"), # Y start 
#             width = unit((diff(col_range) + 14) / ncol(data_matrix), "npc"),   # Adjust width 
#             height = unit((diff(row_range) + 3) / nrow(data_matrix), "npc"),  # Adjust height 
#             just = c("left", "top"), 
#             gp = gpar(col = "darkgreen", fill = NA, lwd = 2))  # Transparent center with red border 
#   
#   grid.rect(x = unit((col_range[1] - 1.5) / ncol(data_matrix), "npc") + unit(0.5 / ncol(data_matrix), "npc"),  # X start 
#             y = unit(1 - (row_range[1] + 9.5) / nrow(data_matrix), "npc") - unit(0.5 / nrow(data_matrix), "npc"), # Y start 
#             width = unit((diff(col_range) + 14) / ncol(data_matrix), "npc"),   # Adjust width 
#             height = unit((diff(row_range) + 4) / nrow(data_matrix), "npc"),  # Adjust height 
#             just = c("left", "top"), 
#             gp = gpar(col = "darkgreen", fill = NA, lwd = 2)) 
#   
#   grid.rect(x = unit((col_range[1] - 1.5) / ncol(data_matrix), "npc") + unit(0.5 / ncol(data_matrix), "npc"),  # X start 
#             y = unit(1 - (row_range[1] + 18.5) / nrow(data_matrix), "npc") - unit(0.5 / nrow(data_matrix), "npc"), # Y start 
#             width = unit((diff(col_range) + 14) / ncol(data_matrix), "npc"),   # Adjust width 
#             height = unit((diff(row_range) + 4) / nrow(data_matrix), "npc"),  # Adjust height 
#             just = c("left", "top"), 
#             gp = gpar(col = "darkgreen", fill = NA, lwd = 2)) 
#   
# })

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------