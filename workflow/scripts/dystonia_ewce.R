#--------------------------------------------------------------------------------------
#
#    Dystonia - EWCE
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# Create EWCE CTD object and check for enrichment of 44 dystonia genes in adult 
# and fetal celltypes across regions

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

# Params
cores <- 8
run_norm <- T
make_plots <- F

# Adult cell types
if (exists("snakemake")) { 

  ##  Load data  ----------------------------------------------------------------------
  set.seed(1234)
  message('Loading data ...')
  adult_object <- readRDS(paste0(R_dir, '03seurat_', region, '.rds'))
  message("Dystonia genes: ", length(dystonia_genes))
  message(paste0(dystonia_genes, collapse = ', '))

  # Recode cluster IDs - Make sure were using whole, joined GeX gene matrix -----------
  message('\nChanging to RNA object ...\n')
  DefaultAssay(adult_object) <- 'RNA'
  adult_object <- JoinLayers(adult_object)
  message('\nSetting idents ...\n')
  Idents(adult_object) <- adult_object$cellIDs
  message('Any NAs in Idents: ', anyNA(Idents(adult_object)))


  ## Prepare annotations  -------------------------------------------------------------
  message('Prepping annotations ...')
  annotations <- adult_object@meta.data |>
    select(level2 = cellIDs) |>
    mutate(level1 = level2) |>
    mutate(level1 = case_when(
      str_detect(level2, "ExN|UBC") ~ str_replace(level2, "-?(ExN|UBC).*", "-ExN"),
      str_detect(level2, "InN") ~ str_replace(level2, "-?(InN).*", "-InN"),
      str_detect(level2, "DaN") ~ str_replace(level2, "-?(DaN).*", "-DaN"),
      str_detect(level2, "Olig") ~ str_replace(level2, "-?(DaN).*", "-Olig"),
      TRUE ~ level2
    )) |>
    relocate(level1)
  head(annotations)
        
  ## Prep Exp matrix and create CTD object  -------------------------------------------
  message('Prepping GeX matrix ...')
  gex_mat <- as(object = adult_object[["RNA"]]$counts, Class = "dgCMatrix")
  
  ## To get FCX through I may have to consider following
  # Pre-filter genes
  # keep_genes <- rowSums(gex_mat > 0) >= 10
  # gex_mat <- gex_mat[keep_genes, ]
  # message("Genes after pre-filtering: ", nrow(gex_mat))
  
  # Downsample cells
  if (region == 'fcx') {
    cells_to_keep <- sample(colnames(gex_mat), size = 150000, replace = FALSE)
    gex_mat <- gex_mat[, cells_to_keep]
    annotations <- annotations[cells_to_keep, , drop = FALSE]
    message("Cells after downsampling: ", ncol(gex_mat))}
  
  if (run_norm) {
    # Normalisation: Note fcx fails here: Error: cannot allocate vector of size 73.4 Gb
    # Tried 8 cores, 370GB; 40 cores, 380GB (latter is pretty much max for a single core)
    message('Normalizing with SCTransform...')
    gex_mat <- EWCE::sct_normalize(gex_mat)
    message("Memory after normalization: ", format(object.size(gex_mat), units = "GB"))
    gc()}
  
  # Drop uninformative genes
  message('Dropping uninformative genes...')
  gex_mat_filt <- EWCE::drop_uninformative_genes(
    exp = gex_mat,
    input_species = "human",
    output_species = "human",
    sctSpecies_origin = "human",
    level2annot = annotations$level2,
    no_cores = cores)
  message("Genes after filtering: ", nrow(gex_mat_filt))
  message("Memory after filtering: ", format(object.size(gex_mat_filt), units = "GB"))
  rm(gex_mat)
  gc()
  
  
  message('Creating ctd ...')
  ctd <- EWCE::generate_celltype_data(exp = gex_mat_filt, 
                                      annotLevels = annotations, 
                                      groupName = ftl_region,
                                      savePath = R_dir,
                                      no_cores = 7)
  rm(gex_mat_filt)
  gc()
  load(ctd)
  
  ## Run GSE Tests  -------------------------------------------------------------------
  message('Running bootstrap GSE lvl 1 ...')
  results_lvl1 <- bootstrap_enrichment_test(
    sct_data = ctd,
    hits = dystonia_genes,
    bg = rownames(ctd[[1]]$specificity), # All genes
    genelistSpecies = "human",
    sctSpecies = "human",
    reps = 10000,
    no_cores = cores,
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
    no_cores = cores,
    annotLevel = 2,
    store_gene_data = FALSE
  )
  
  ## Extract tables and plots  --------------------------------------------------------
  message('Running bootstrap GSE lvl 2 ...')
  write_tsv(results_lvl1$results |> as_tibble(), paste0(table_dir, region, '_adult_ewce_lvl1.tsv'))
  write_tsv(results_lvl2$results |> as_tibble(), paste0(table_dir, region, '_adult_ewce_lvl2.tsv'))
  
}

# Fetal cell types
if (Sys.info()[["nodename"]] == "Darrens-iMac-2.local") {
 
  for (ftl_region in c('fcx', 'ge', 'cer')) {
    
    message('Running EWCE for ', ftl_region , ' ...\n')
    fetal_object <- readRDS(paste0(fetal_dir , 'seurat_', ftl_region, '_fetal.rds'))
    
    message('Prepping annotations ...\n')
    annotations <- fetal_object@meta.data |>
      dplyr::select(level2 = cellIDs) |>
      mutate(level1 = level2)
    head(annotations)
    
    message('Normalizing with SCTransform...\n')
    gex_mat <- EWCE::sct_normalize(fetal_object[["RNA"]]$counts)
    
    message('\nDropping uninformative genes...\n')
    gex_mat_filt <- EWCE::drop_uninformative_genes(
      exp = gex_mat,
      input_species = "human",
      output_species = "human",
      sctSpecies_origin = "human",
      level2annot = annotations$level2,
      no_cores = cores)
    message("Genes after filtering: ", nrow(gex_mat_filt))
    message("Memory after filtering: ", format(object.size(gex_mat_filt), units = "GB"))
    rm(gex_mat)
    gc()
    
    message('\nCreating ctd ...\n')
    ctd <- EWCE::generate_celltype_data(exp = gex_mat_filt, 
                                        annotLevels = annotations, 
                                        groupName = ftl_region,
                                        savePath = R_dir,
                                        no_cores = cores)
    rm(gex_mat_filt)
    
    load(ctd)
    
    message('\nRunning bootstrap GSE lvl 2 ...\n')
    results_lvl2 <- bootstrap_enrichment_test(
      sct_data = ctd,
      hits = dystonia_genes,
      bg = rownames(ctd[[1]]$specificity), # All genes
      genelistSpecies = "human",
      sctSpecies = "human",
      reps = 10000,
      no_cores = cores,
      annotLevel = 2,
      store_gene_data = FALSE
    )
    
    message('\nWriting results ...\n')
    write_tsv(results_lvl2$results |> as_tibble(), paste0(table_dir, ftl_region, '_fetal_ewce_lvl2.tsv'))
    
  }
   
}
  
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------