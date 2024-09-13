#--------------------------------------------------------------------------------------
#
#    Dystonia - Prepare Brain Span Data
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# exp_tbl: rows = genes; cols = samples; 1st col is the row number
# col_meta_tbl: the genes are listed in the same order as the rows in expression_matrix.csv
# row_meta_tbl: the samples are listed in the same order as the columns in expression_matrix

# Samples: 524
# Genes: 52376

##  Load Packages  --------------------------------------------------------------------
library(dlookr)
library(GWalkR)
library(tidyverse)
library(tabulapdf) # Extract tables from pdf
library(gtsummary)
library(readxl)
library(pheatmap)
library(cowplot)
library(reshape2)
library(performance) # LM diagnostic plots
library(explore)

root_dir <- '~/Desktop/dystonia_snRNAseq_2024/'
results_dir <- paste0(root_dir, 'results/')
bulk_dir <- paste0(results_dir, '02Bulk_data/')
resources_dir <- paste0(root_dir, 'resources/sheets/')
brain_span_dir <- paste0(root_dir, 'resources/public_data/brain_span/genes_matrix_csv/')
dystonia_genes <- c("EIF2AK2", "GNAO1", "KCTD17", "TUBB4A", "ATP1A3", "KMT2B", 
                    "DNAJC12", "KCNA1", "SPR", "SLC2A1", "HPCA", "PNKD", "SGCE",
                    "THAP1", "GCH1", "ANO3", "TOR1A", "GNAL", "KCNMA1", "PRRT2",
                    "ADCY5", "TH", "PRKRA", "SCN8A", "VPS16")

## Set variables  ---------------------------------------------------------------------
root_dir <- '~/Desktop/dystonia_snRNAseq_2024/'
sheets_dir <- paste0(root_dir, 'resources/sheets/')
brain_span_dir <- paste0(root_dir, 'resources/public_data/brain_span/genes_matrix_csv/') 
dystonia_genes <- c("EIF2AK2", "GNAO1", "KCTD17", "TUBB4A", "ATP1A3", "KMT2B", 
                    "DNAJC12", "KCNA1", "SPR", "SLC2A1", "HPCA", "PNKD", "SGCE",
                    "THAP1", "GCH1", "ANO3", "TOR1A", "GNAL", "KCNMA1", "PRRT2",
                    "ADCY5", "TH", "PRKRA", "SCN8A", "VPS16")
region_levels <- c("PFC", "PMSC", "NPFC", "Str", "Cer", "Hip", "Tha", "Amy")
dev_stage_levels <- c("EarlyFetal", "LateFetal", "MidFetal", "Infancy", "Childhood", 
                      "Adolescence", "Adulthood")

## Read in data  ----------------------------------------------------------------------
reads_tbl <- read_csv(paste0(brain_span_dir, 'expression_matrix.csv'), col_names = FALSE)
col_meta_tbl <- read_csv(paste0(brain_span_dir, 'columns_metadata.csv')) |>
  mutate(col_age = str_replace_all(age, c(' pcw' = 'PCW',  # Recode age to match across tbls
                                          ' mos' = 'M', 
                                          ' yrs' = 'Y'))) |>
  rename(col_donor_name = donor_name,
         col_gender = gender) |>
  relocate(column_num, col_donor_name, structure_acronym, structure_name, col_age) |>
  rename(col_meta_id = structure_acronym,
         col_meta_region_name = structure_name) |> 
  mutate(peall_region = case_match(
    col_meta_id ,
    'A1C' ~ 'NPFC',
    'AMY' ~ 'Amy',
    'CB' ~ 'Cer',
    'CBC' ~ 'Cer',
    'CGE' ~ 'Str',
    'DFC' ~ 'PFC',
    'DTH' ~ 'Tha',
    'HIP' ~ 'Hip',
    'IPC' ~ 'NPFC',
    'ITC' ~ 'NPFC',
    'LGE' ~ 'Str',
    'M1C' ~ 'PMSC',
    'M1C-S1C' ~ 'PMSC',
    'MD' ~ 'Tha',
    'MFC' ~ 'PFC',
    'MGE' ~ 'Str',
    'OFC' ~ 'PFC',
    'Ocx' ~ 'NPFC',
    'PCx' ~ 'PMSC', 
    'S1C' ~ 'PMSC',
    'STC' ~ 'NPFC',
    'STR' ~ 'Str',
    'TCx' ~ 'NPFC',  
    'URL' ~ 'Cer',
    'V1C' ~ 'NPFC',
    'VFC' ~ 'PFC')) 
  
row_meta_tbl <- read_csv(paste0(brain_span_dir, 'rows_metadata.csv'))
pdf_tbl <- extract_tables(paste0(brain_span_dir, 'Human_Brain_Seq_Stages_Sample_Metadata.pdf'))
ethnicity_tbl <- bind_rows(pdf_tbl[[1]] |> select(-PMI), pdf_tbl[[2]] |> select(-PMI)) |>
  rename(donor_name = `Internal ID`) |>
  mutate(eth_age = str_replace(Age, ' ', '')) |>  # Recode age to match across tbls
  mutate(eth_age = str_replace(eth_age, 'PM', 'M')) |>
  mutate(age = str_replace(eth_age, '22PCW', '24PCW')) |> # Mismatch in col, rin, eth
  select(-`External ID`, -Age, -pH) |>
  rename(eth_gender = Gender, 
         ethnicity = Ethn.,
         eth_donor_id = donor_name)
  
## Create rin tbl from groups of PDFs  ------------------------------------------------
# Raw rin_tbl: 575 rows. col_meta_tbl: 524 rows.
rin_tbl <- tibble(brain_code = character(),
                  donor_id = character(),
                  age = character(),
                  region = character(),
                  hemisphere = character(),
                  rin = numeric(),
                  dissection_score = numeric())

for (i in seq(3, 17, 1)) {
  
  tbl_raw <- pdf_tbl[[i]]
  tbl_clean <- tbl_raw |>
    tail(-2) |>
    select(...3:...9) |>
    rename(brain_code = ...3,
           donor_id = ...4,
           age = ...5,
           region = ...6,
           hemisphere = ...7,
           rin = ...8,
           dissection_score = ...9) |>
    drop_na()
  
  rin_tbl <- bind_rows(rin_tbl, tbl_clean)
  
}

# Recode to match col_meta_tbl for join: 1st attempt  ---------------------------------
# Entries for 35PCW and 37PCW in col_meta_tbl but not rin_tbl (2 donors / 18 dissections)
# Donor H376.IV.50 24PCW in col tbl, 22PCW in rin an eth tbl: Set to 24PCW (all 2nd tri anyway)
# Recode region to match col_meta_tbl
rin_tbl <- rin_tbl |> 
  mutate(age = str_replace(age, '12M', '1Y')) |> # Recode months to year
  mutate(age = str_replace(age, '22PCW', '24PCW')) |> # Mismatch in col, rin, eth
  rename(rin_age = age,
         rin_donor_id = donor_id,
         rin_region = region) |>
  select(-dissection_score, -hemisphere, -brain_code) |>
  mutate(rin_row_num = row_number()) |>
  mutate(rin_region_recode = case_match(
    rin_region,
    'A1C' ~ 'A1C',
    'AMY' ~ 'AMY',
    'C-PC (M1C/S1C' ~ 'M1C-S1C',
    'CBC' ~ 'CBC',
    'CGE' ~ 'CGE',
    'DFC' ~ 'DFC',
    'DTH' ~ 'DTH',
    'FC (DFC)' ~ 'DFC',
    'FC (MFC)' ~ 'MFC',
    'FC (OFC)' ~ 'OFC',
    'FC (VFC)' ~ 'VFC',
    'HIP' ~ 'HIP',
    'IPC' ~ 'IPC',
    'ITC' ~ 'ITC',
    'LGE' ~ 'LGE',
    'M1C' ~ 'M1C',
    'M1C/S1C' ~ 'M1C-S1C',
    'MD' ~ 'MD',
    'MFC' ~ 'MFC',
    'MGE' ~ 'MGE',
    'OC' ~ 'Ocx',
    'OC (V1C)' ~ 'Ocx',     # Ambiguous
    'OFC' ~ 'OFC',
    'PC (IPC)' ~ 'PCx',  # Ambiguous
    'S1C' ~ 'S1C',
    'STC' ~ 'STC',
    'STR' ~ 'STR',
    'TC (A1C/STC)' ~ 'TC',  # Ambiguous
    'TC (ITC)' ~ 'TC', # Ambiguous
    'URL' ~ 'URL',
    'V1C' ~ 'V1C',
    'VFC' ~ 'VFC')) |>
  relocate(rin_row_num, rin_donor_id, rin_region, rin_region_recode)

# Join col meta and rin tbls  ---------------------------------------------------------
# Only 451 samples after first join: Check NAs and recoding
# 481 samples in KPs tbl
inner_join_tbl <- col_meta_tbl |>
  inner_join(rin_tbl, by = join_by(col_donor_name == rin_donor_id,
                                   col_meta_id  == rin_region_recode)) 

# 124 samples in rin tbl not in col tbl
right_join_tbl <- col_meta_tbl |>
  right_join(rin_tbl, by = join_by(col_donor_name == rin_donor_id,
                                   col_meta_id  == rin_region_recode)) |>
  relocate(column_num, rin_row_num, col_donor_name, col_meta_id, col_meta_region_name) |>
  filter(if_any(everything(), is.na)) |>
  select(rin_row_num, col_donor_name, col_meta_id) |>
  rename(rin_region_recode = col_meta_id)

# 73 samples in col tbl but not in rin tbl
left_join_tbl <- col_meta_tbl |>
  left_join(rin_tbl, by = join_by(col_donor_name == rin_donor_id,
                                   col_meta_id  == rin_region_recode)) |>
  filter(if_any(everything(), is.na)) |>
  relocate(column_num, rin_row_num, col_donor_name, col_meta_id, col_meta_region_name) 


# Pull unique donors from rin table as those donors def in both rin and col tbls
donors_rin <- rin_tbl |> distinct(rin_donor_id) |> pull()
donors_col <- col_meta_tbl |> distinct(col_donor_name) |> pull()

# Walk through rin and col donors to catch mismatches
for (i in donors_rin) {
  
  # Get missing row indexes for col_meta_tbl and rin_tbl
  
  col_missing <-  left_join_tbl |> 
    filter(col_donor_name == i) |>
    pull(column_num)
  
  rin_missing <-  right_join_tbl |> 
    filter(col_donor_name == i) |>
    pull(rin_row_num) 
  
  message('Col Meta missing:')
  col_meta_tbl |> 
    filter(col_donor_name == i) |>
    filter(column_num %in% col_missing) |>
    print(n = Inf)
  
  message('Rin missing:')
  rin_tbl |> 
    filter(rin_donor_id == i) |>
    filter(rin_row_num %in% rin_missing) |>
    print(n = Inf)
  
  message('Col join missing:')
  print(left_join_tbl |> filter(col_donor_name == i), n = Inf)
  
  message('Rin join missing:')
  print(right_join_tbl |> filter(col_donor_name == i), n = Inf)
  
}

## Option to recode variables  --------------------------------------------------------
read_xlsx(paste0(sheets_dir, 'dystonia_col_rin_table_mismatch.xlsx'))

## Remove RIN: 438 rows; 39 donors; Join with ethinicity tbl  -------------------------
all_join_tbl <- inner_join_tbl |>
  filter(rin >= 7) |>
  inner_join(ethnicity_tbl, by = join_by(col_donor_name == eth_donor_id)) 
  
# Compare the brain region annotations and counts across col and rin tables
col_region_cnts <- col_meta_tbl |> 
  group_by(col_meta_id, col_meta_region_name, peall_region) |> 
  count(name = 'col_meta_n') 

rin_region_cnts <- rin_tbl |> 
  group_by(rin_region, rin_region_recode) |> 
  count(name = 'rin_meta_n')

## Annotate gene meta to expr table  --------------------------------------------------
# Need to recode MLL2 to KMT2B (MLL4 also present in data)
reads_tbl2 <- reads_tbl |>
  inner_join(row_meta_tbl, by = join_by('X1' == 'row_num')) |>
  filter(gene_symbol %in% c(dystonia_genes, 'MLL2')) |>
  mutate(gene_symbol = recode(gene_symbol, "MLL2" = "KMT2B")) |>
  relocate(gene_symbol, gene_id, ensembl_gene_id, entrez_id) |>
  dplyr::select(-gene_id, -ensembl_gene_id, -entrez_id, -X1) |>
  as.data.frame()

## Transpose df to join col_meta_data -------------------------------------------------
rownames(reads_tbl2) <- reads_tbl2$gene_symbol
expr_tbl <- reads_tbl2 |>
  dplyr::select(-gene_symbol) |>
  t() |>
  as.data.frame() |>
  as_tibble() |>
  mutate(column_num = row_number()) |>
  inner_join(all_join_tbl, by = 'column_num') |>
  relocate(column_num:eth_age) |>
  mutate(dev_stage = case_match(col_age,
    c('8PCW', '9PCW') ~ 'EarlyFetal',
    c('12PCW', '13PCW', '16PCW', '17PCW', '19PCW', '21PCW') ~ 'MidFetal',
    c('24PCW', '25PCW', '26PCW') ~ 'LateFetal', # Check donor H376.IV.50 recode
    c('4M', '10M') ~ 'Infancy',
    c('1Y', '2Y', '3Y', '4Y', '8Y', '11Y') ~ 'Childhood',
    c('13Y', '15Y', '18Y', '19Y') ~ 'Adolescence',
    c('21Y', '23Y', '30Y', '36Y', '37Y', '40Y') ~ 'Adulthood')) |>
  select(donor = col_donor_name,
         region = peall_region,
         dev_stage,
         ethnicity,
         rin,
         sex = col_gender,
         EIF2AK2:VPS16) |>
  mutate(region = factor(region, levels = region_levels)) |>
  mutate(dev_stage = factor(dev_stage, levels = dev_stage_levels))

expr_tbl |> select(region) |> group_by(region) |> count() 
expr_tbl |> select(dev_stage) |> group_by(dev_stage) |> count()

saveRDS(expr_tbl, paste0(bulk_dir, 'dystonia_brain_span_clean_tbl.rds'))
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------