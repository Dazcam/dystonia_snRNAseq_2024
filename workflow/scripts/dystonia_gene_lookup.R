#--------------------------------------------------------------------------------------
#
#    Dystonia genes lookup
#
#--------------------------------------------------------------------------------------

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

lookup_hg19 <- get_biomart_gene_lookup('hg19')
lookup_hg38 <- get_biomart_gene_lookup('hg38')

# Get lookup table for initial 25 gene list  ------------------------------------------
dystonia_genes_25 <- read_excel(paste0(data_dir, 'sheets/Dystonia_Genes_Clinical_5.0.xlsx'), 
                                range = 'D1:D26') %>%
  arrange(GeneName) %>%
  pull(GeneName)

# Note 2 Ensembl IDs for ATP1A3 in hg19 and different encodings for KMT2B in hg19 and hg38
# Used ENSG00000105409 which is in both hg19 and hg38 (other one is ENSG00000272584 ony in hg19)
#KMT2B_hg19_encoding <- "ENSG00000105663" # Encoded as a processed_transcript in hg19
#KMT2B_hg38_encoding <- "ENSG00000272333" # Encoded as a protein_coding gene in hg38
dystonia_genes_25_hg19 <- lookup_hg19 |> 
  filter(hgnc_symbol %in% dystonia_genes_25) |>
  arrange(hgnc_symbol) |>
  dplyr::select(hgnc_hg19 = hgnc_symbol, 
         ensembl_hg19 = ensembl_gene_id, 
         bioytpe_hg19 = gene_biotype) |>
  distinct() |>
  print(n = Inf) 

dystonia_genes_25_hg38 <- lookup_hg38 |>
  filter(hgnc_symbol %in% dystonia_genes_25) |>
  #filter(ensembl_gene_id %in% dystonia_genes_hg19_ensembl) |>
  arrange(hgnc_symbol) |>
  dplyr::select(hgnc_hg38 = hgnc_symbol, 
                ensembl_hg38 = ensembl_gene_id, 
                bioytpe_hg38 = gene_biotype) |>
  distinct() |>
  print(n = Inf) 

dystonia_genes_25_tbl <- dystonia_genes_25_hg19 |>
  inner_join(dystonia_genes_25_hg38, by = join_by(hgnc_hg19 == hgnc_hg38)) |>
  filter(!ensembl_hg19 == 'ENSG00000272584') |>
  print(n = Inf) 

# Get lookup table for larger 44 gene list  -------------------------------------------
dystonia_genes_44 <- read_excel(paste0(data_dir, 'sheets/Updated_Dystonia_Gene_List_040225.xlsx'), 
                                range = 'J3:J46', col_names = 'GeneName') %>%
  arrange(GeneName) %>%
  separate(GeneName, into = c('GeneName', 'Alias')) %>%
  pull(GeneName)

dystonia_genes_44_alias <- read_excel(paste0(data_dir, 'sheets/Updated_Dystonia_Gene_List_040225.xlsx'), 
                                range = 'J3:J46', col_names = 'GeneName') %>%
  arrange(GeneName) %>%
  separate(GeneName, into = c('GeneName', 'Alias')) %>%
  drop_na() 

# ATP1A3, BCAP31, TIMM8A two ensembl ids in hg19 
# AOPEP (alias c9orf3) not in hg19
dystonia_genes_44_hg19 <- lookup_hg19 |> 
  filter(hgnc_symbol %in% dystonia_genes_44) |>
  arrange(hgnc_symbol) |>
  dplyr::select(hgnc_hg19 = hgnc_symbol, 
                ensembl_hg19 = ensembl_gene_id, 
                bioytpe_hg19 = gene_biotype) |>
  distinct() |>
  print(n = Inf) 

# SLC6A3, SQSTM1 two ensembl ids in hg38
# AOPEP in hg38
dystonia_genes_44_hg38 <- lookup_hg38 |>
  filter(hgnc_symbol %in% dystonia_genes_44) |>
  #filter(ensembl_gene_id %in% dystonia_genes_hg19_ensembl) |>
  arrange(hgnc_symbol) |>
  dplyr::select(hgnc_hg38 = hgnc_symbol, 
                ensembl_hg38 = ensembl_gene_id, 
                bioytpe_hg38 = gene_biotype) |>
  distinct() |>
  print(n = Inf) 

dystonia_genes_44_tbl <- dystonia_genes_44_hg38 |> 
  left_join(dystonia_genes_44_hg19, by = join_by(hgnc_hg38 == hgnc_hg19), relationship = "many-to-many") |> 
  dplyr::filter(!(ensembl_hg19 %in% c('ENSG00000272584', # ATP1A3: kept ENSG00000105409
                                      'ENSG00000267977', # BCAP31: kept ENSG00000185825
                                      'ENSG00000268249'  # TIMM8A: kept ENSG00000126953
                                      )) | is.na(ensembl_hg19)) |> 
  dplyr::filter(!(ensembl_hg38 %in% c('ENSG00000276996', # SLC2A1: kept ENSG00000142319
                                      'ENSG00000284099'  # SQSTM1: kept ENSG00000161011
                                      )) | is.na(ensembl_hg19)) |> 
  print(n = Inf)

# Genes for additional GO analyses
omit_genes <- c('SPR', 'GCH1', 'TH')
dystonia_genes_44_tbl |> 
  filter(!hgnc_hg38 %in% omit_genes)|>
  pull(ensembl_hg38) |>
  cat(sep = '\n')

lookup_hg38 |>
  filter(hgnc_symbol %in% dystonia_genes_44_alias$Alias) |>
  #filter(ensembl_gene_id %in% dystonia_genes_hg19_ensembl) |>
  arrange(hgnc_symbol) |>
  dplyr::select(hgnc_hg38 = hgnc_symbol, 
                ensembl_hg38 = ensembl_gene_id, 
                bioytpe_hg38 = gene_biotype) |>
  distinct() |>
  print(n = Inf) 

dystonia_genes <- dystonia_genes_44_tbl |> pull(hgnc_hg38)
dput(dystonia_genes)
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
