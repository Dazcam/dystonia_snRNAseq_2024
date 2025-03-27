#--------------------------------------------------------------------------------------
#
#    Dystonia - EWCE
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# Create EWCE CTD object and check for enrichment of 44 dystonia genes

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

##  Load data  ------------------------------------------------------------------------
set.seed(1234)
message('Loading data ...')
adult_object <- readRDS(paste0(R_dir, '03seurat_', region, '.rds'))

# Recode cluster IDs - Make sure were using whole, joined GeX gene matrix -------------
message('\nChanging to RNA object ...\n')
DefaultAssay(adult_object) <- 'RNA'
adult_object <- JoinLayers(adult_object)
message('\nSetting idents ...\n')
Idents(adult_object) <- adult_object$cellIDs
message('Any NAs in Idents: ', anyNA(Idents(adult_object)))


## Prepare annotations  ---------------------------------------------------------------
message('Prepping annotations ...')
annotations <- adult_object@meta.data |>
  dplyr::select(level2 = cellIDs) |>
  mutate(level1 = level2) |>
  mutate(level1 = case_when(
    str_detect(level2, "ExN|UBC") ~ str_replace(level1, "^([A-Za-z]+(?:_[a-z]+)?-[A-Za-z]+).*", "\\1"),
    str_detect(level2, "InN") ~ str_replace(level1, "^([A-Za-z]+(?:_[a-z]+)?-[A-Za-z]+).*", "\\1"),
    TRUE ~ level2
  )) 
head(annotations)
      
## Prep Exp matrix and create CTD object  ---------------------------------------------
message('Prepping GeX matrix ...')
gex_mat <- as(object = adult_object[["RNA"]]$counts, Class = "dgCMatrix")
gex_mat_norm <- EWCE::sct_normalize(gex_mat)
gex_mat_filt <- EWCE::drop_uninformative_genes(
  exp = gex_mat_norm,
  input_species = "human",
  output_species = "human",
  level2annot = annotations$level2,
  no_cores = 7)

message('Creating ctd ...')
ctd <- EWCE::generate_celltype_data(exp = gex_mat_filt, 
                                    annotLevels = annotations, 
                                    groupName = region,
                                    savePath = R_dir,
                                    no_cores = 7)


## Run GSE Tests  ---------------------------------------------------------------------
message('Running bootstrap GSE lvl 1 ...')
results_lvl1 <- bootstrap_enrichment_test(
  sct_data = ctd,
  hits = dystonia_genes,
  bg = rownames(ctd[[1]]$specificity), # All genes
  genelistSpecies = "human",
  sctSpecies = "human",
  reps = 10000,
  no_cores = 7,
  annotLevel = 1,
  store_gene_data = FALSE
)

message('Running bootstrap GSE lvl 2 ...')
results_lvl2 <- bootstrap_enrichment_test(
  sct_data = ctd,
  hits = dystonia_genes,
  bg = rownames(ctd[[1]]$specificity), # All genes
  genelistSpecies = "human",
  sctSpecies = "human",
  reps = 10000,
  no_cores = 7,
  annotLevel = 2,
  store_gene_data = FALSE
)

## Extract tables and plots  ----------------------------------------------------------
message('Running bootstrap GSE lvl 2 ...')
write_tsv(results_lvl1$results |> as_tibble(), paste0(table_dir, region, '_adult_ewce_lvl1.tsv'))
write_tsv(results_lvl2$results |> as_tibble(), paste0(table_dir, region, '_adult_ewce_lvl2.tsv'))


lvl1_plt <- EWCE::ewce_plot(total_res = results_lvl1$results, mtc_method = "BH")$plain +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"))
lvl2_plt <- EWCE::ewce_plot(total_res = results_lvl2$results, mtc_method = "BH")$plain +
  theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"))
comb_plt <- plot_grid(lvl1_plt + ggtitle('Level 1'), lvl2_plt + ggtitle('Level 2'))
ggsave(paste0(fig_dir, region, '_ewce.png'), comb_plt, width = 7, height = 5, 
       dpi = 300, units = "in", bg = "white")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------