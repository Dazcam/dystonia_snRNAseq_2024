#--------------------------------------------------------------------------------------
#
#    Dystonia - gene lists and colur schemes
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------
# Striatal scRNAseq (mouse). Anderson et al (2023). PMID: 37270616

## General gene list  -----------------------------------------------------------------
general_genes <- c('SLC17A7', 'SLC17A6', 'SLC17A8', # VGLUT1-3
                   'SLC6A1', 'SLC6A13', 'SLC6A11', 'SLC6A12', # GABA transporters
                   'SST', 'NPY', 'GAD1', 'GAD2', 'PVALB', 'CALB2', 'VIP', # InN markers
                   'C3', 'C1QB',                           # MG markers
                   'AQP4', 'SOX9', 'GFAP',        # Astrocytes 
                   'OLIG1', 'OLIG2', 'MBP',                # Oligodendrocytes
                   'PDGRFA', 'PMP2',
                   'EOMES', 'EBF1', 'ABCB1') 

## Str  -------------------------------------------------------------------------------
str_genes <- c('DRD1', 'DRD2', 'TAC1', 'PENK',
               'FOXP1', 'MYT1L', 'MEIS2', 'CSDE1',
               'SOX11', 'BCL11B', 'YBX1', 'EBF1',
               'RARB', 'NR1D1', 'TEF', 'NOS1',
               'SST', 'NPY', 'GAD1', 'GAD2', 'PVALB', 'CALB2', 'VIP',
               'CCK', 'TH', 'CHAT', 'CALB1', 'SOX4', 'DLX2',
               'LHX6', 'LHX8', 'NKX2.1', 'NR2F2')

str_final_genes <- factor(c('GAD1', 'GAD2', 'SLC6A1', 'DRD1', 'DRD2', 
                            'SST', 'NPY', 'CALB2', 'VIP', 'NPY', 'CHAT', 
                            'SLC17A7','GFAP', 'AQP4', 'SOX9', 'OLIG1', 
                            'OLIG2', 'MBP', 'OLIG1', 'OLIG2', 'C3', 'C1QB'))

str_umap_cols_recode <- c("Str-adult-InN-1" = '#708238', "Str-adult-InN-2" = '#31C53F', "Str-adult-InN-3" = '#3CBB75FF', 
                          "Str-adult-InN-4" = '#9DC183', "Str-adult-InN-5" = '#006400', "Str-adult-InN-6" = '#95D840FF', 
                          "Str-adult-InN-7" = '#2FF18B', "Str-adult-InN-8" = '#B7FFB7', "Str-adult-InN-9" = '#3CBB75FF',
                          "Str-adult-Ast-1" = '#FF5959', "Str-adult-Ast-2" = '#FF5959', "Str-adult-Olig-1" = '#FDE725FF', 
                          "Str-adult-OPC" = '#FDE725FF', "Str-adult-MG" = '#F58231', "Str-adult-ExN" = '#00B6EB', 
                          "Str-adult-Misc-1" = '#CCCCCC')

str_vln_cols_recode <- c("Str-adult-InN-1" = '#3CBB75FF', "Str-adult-InN-2" = '#3CBB75FF', "Str-adult-InN-3" = '#3CBB75FF', 
                         "Str-adult-InN-4" = '#3CBB75FF', "Str-adult-InN-5" = '#3CBB75FF', "Str-adult-InN-6" = '#3CBB75FF', 
                         "Str-adult-InN-7" = '#3CBB75FF', "Str-adult-InN-8" = '#3CBB75FF', "Str-adult-InN-9" = '#3CBB75FF',
                         "Str-adult-Ast-1" = '#FF5959', "Str-adult-Ast-2" = '#FF5959', "Str-adult-Olig-1" = '#FDE725FF', 
                         "Str-adult-OPC" = '#FDE725FF', "Str-adult-MG" = '#F58231', "Str-adult-ExN" = '#00B6EB', 
                         "Str-adult-Misc-1" = '#CCCCCC')

## Cer -----
cer_genes <- c('MFAP4', 'MGP', 'RBFOX3', 'RELN',         # Granule cells
               'NEFM', 'NEUROD1',                        # Ex cer nuclei cell (eCN)
               'PCP4', 'RORA', 'CA8', 'ITPR1',           # Purkinje cells  
               'GAD1', 'GAD2', 'PAX2',                   # GABA Ns
               'MBP', 'APOD', 'OLIG1', 'OLIG2',          # Oligodendrocytes
               'GDF10', 'PAX3',                          # Bergmann cells
               'AQP4', 'SOX9', 'GFAP',                   # Astrocytes 
               'AIF1', 'C1QB',                           # Microglia 
               'EOMES')                                  # Unipolar brush cell (UBC)

fcx_genes <- c('CUX3', #L2/3
               'RORB', #L4
               'TLE4', #L5/6
               'GFAP', 'SLC1A2', #Ast
               'MBP', #Oligs
               "GLI3", "OLIG1",
               "MKI67", "C3", "ITM2A", "SST", "CALB2", 
               "SCGN", "TLE3", "FEZF2", "CRYM", "LHX2")

fcx_genes <- c('SLC17A7', 'SLC17A6', 'SLC17A8', # VGLUT1-3
               'CUX3', #L2/3
               'RORB', #L4
               'TLE4', #L5/6
               'SLC6A1', 'SLC6A13', 'SLC6A11', 'SLC6A12', # GABA transporters
               'SST', 'NPY', 'GAD1', 'GAD2', 'PVALB', 'CALB2', 'VIP', # InN markers
               'C3', 'C1QB',            # MG markers
               'AQP4', 'SOX9', 'GFAP',         # Astrocytes 
               'OLIG1', 'OLIG2', 'MBP',        # Oligodendrocytes
               'PDGRFA', 'PMP2',
               'EOMES', 'EBF1', 'ABCB1') 

# Colours
greens <- c('#3CBB75FF', '#00FF00A5','#006400', '#B7FFB7', '#10A53DFF',
            '#95D840FF', '#9DC183',  '#708238', '#55C667FF', '#73D055FF',
            '#567D46')

ExN_blues <- c('#76B5C5', '#00BDD2', '#CEE5FD', '#00B6EB', '#ABDBE3',
               '#1E81B0', '#B8D2EB', '#779CBA')

reds <- c('#FAA0A0', '#FF5959', '#F75151', '#EF0029', '#D2042D')

reds <- c('#B200ED',  '#DCBEFF', '#6F2DA8')

# Striatum. ---------------------------------------------------------------------------
# https://cssgradient.io/shades-of-green/
str_umap_cols <- c('0' = '#3CBB75FF','1' = '#31C53F','2' = '#708238', '3' = '#B7FFB7','4' = '#006400',
                   '5' = '#FDE725FF', '6' = '#95D840FF', '7'='#2FF18B', '8' = '#9DC183','9' = '#FF5959', 
                   '10'='#3CBB75FF', '11'='#31C53F', '12'='#708238', '13'='#FDE725FF', '14'='#006400',
                   '15'='#95D840FF', '16'='#00B6EB', '17' = '#F58231', '18'='#2FF18B', '19'='#9DC183', 
                   '20'='#00B6EB', '21'='#0CB702', '22'='#00BE67', '23'='#FF5959', '24'='#7CAE00', 
                   '25' = '#FDE725FF')

str_vln_cols <- c('0' = '#3CBB75FF','1' = '#3CBB75FF','2' = '#3CBB75FF','3' = '#3CBB75FF','4' = '#3CBB75FF',
                  '5' = '#FDE725FF', '6' = '#3CBB75FF', '7'='#3CBB75FF', '8' = '#3CBB75FF','9' = '#FF5959', 
                  '10'='#3CBB75FF', '11'='#3CBB75FF', '12'='#3CBB75FF', '13'='#FDE725FF', '14'='#3CBB75FF',
                  '15'='#3CBB75FF', '16'='#CEE5FD', '17' = '#F58231', '18'='#3CBB75FF', '19'='#3CBB75FF', 
                  '20'='#CEE5FD', '21'='#3CBB75FF', '22'='#3CBB75FF', '23'='#FF5959', '24'='#3CBB75FF', 
                  '25' = '#FDE725FF')

if (exists('seurat_fcx')) {
# Frontal cortex
fcx_clusters_recode <- seurat_sk_fcx@meta.data %>% 
  as_tibble() %>%
  mutate(harmony_clusters_recode = recode(harmony_clusters, 
                                          `0` = "FCX-adult-ExN-1", 
                                          `1` = "FCX-adult-ExN-2",
                                          `2` = "FCX-adult-InN-1",
                                          `3` = "FCX-adult-ExN-3",
                                          `4` = "FCX-adult-Olig-1",
                                          `5` = "FCX-adult-InN-2",
                                          `6` = "FCX-adult-InN-3",
                                          `7` = "FCX-adult-Ast-1",
                                          `8` = "FCX-adult-InN-4",
                                          `9` = "FCX-adult-ExN-4",
                                          `10` = "FCX-adult-ExN-5",
                                          `11` = "FCX-adult-InN-5",
                                          `12` = "FCX-adult-InN-6",
                                          `13` = "FCX-adult-Olig-2",
                                          `14` = "FCX-adult-ExN-6",
                                          `15` = "FCX-adult-ExN-7",
                                          `16` = "FCX-adult-InN-7",
                                          `17` = "FCX-adult-ExN-8",
                                          `18` = "FCX-adult-InN-8",
                                          `19` = "FCX-adult-ExN-9",
                                          `20` = "FCX-adult-ExN-10",
                                          `21` = "FCX-adult-ExN-11",
                                          `22` = "FCX-adult-MG-1",
                                          `23` = "FCX-adult-ExN-12",
                                          `24` = "FCX-adult-InN-9",
                                          `25` = "FCX-adult-ExN-13",
                                          `26` = "FCX-adult-InN-10",
                                          `27` = "FCX-adult-Ast-2",
                                          `28` = "FCX-adult-InN-11",
                                          `29` = "FCX-adult-InN-12",
                                          `30` = "FCX-adult-InN-13",
                                          `31` = "FCX-adult-InN-14",
                                          `32` = "FCX-adult-ExN-14",
                                          `33` = "FCX-adult-Endo-??",
                                          `34` = "FCX-adult-ExN-15",
                                          `35` = "FCX-adult-Olig-3")) %>%
  pull(harmony_clusters_recode)

fcx_umap_cols_recode <- c("FCX-adult-ExN-1" = '#00B6EB' , "FCX-adult-ExN-2" = '#CEE5FD', "FCX-adult-InN-1" = '#3CBB75FF',
                          "FCX-adult-ExN-3" = "#00BDD2", "FCX-adult-Olig-1" = '#FDE725FF', "FCX-adult-InN-2" = '#006400', 
                          "FCX-adult-InN-3" = '#2FF18B', "FCX-adult-Ast-1" = '#FF5959', "FCX-adult-InN-4" = '#95D840FF', 
                          "FCX-adult-ExN-4" = '#779CBA', "FCX-adult-ExN-5" = '#76B5C5', "FCX-adult-InN-5" = '#708238', 
                          "FCX-adult-InN-6" = '#00BE67', "FCX-adult-Olig-2" = '#FDE725FF', "FCX-adult-ExN-6" = '#ABDBE3', 
                          "FCX-adult-ExN-7" = '#1E81B0', "FCX-adult-InN-7" = '#31C53F', "FCX-adult-ExN-8" = '#76B5C5', 
                          "FCX-adult-InN-8" = '#10A53DFF', "FCX-adult-ExN-9" = '#76B5C5', "FCX-adult-ExN-10" = "#00BDD2", 
                          "FCX-adult-ExN-11" = '#CEE5FD', "FCX-adult-MG-1" = '#F58231', "FCX-adult-ExN-12" = '#00B6EB', 
                          "FCX-adult-InN-9" = '#00FF00A5', "FCX-adult-ExN-13" = '#779CBA', "FCX-adult-InN-10" = '#B7FFB7', 
                          "FCX-adult-Ast-2" = '#EF0029', "FCX-adult-InN-11" = '#9DC183', "FCX-adult-InN-12" = '#31C53F', 
                          "FCX-adult-InN-13" = '#3CBB75FF', "FCX-adult-InN-14" = '#006400', "FCX-adult-ExN-14" = '#ABDBE3', 
                          "FCX-adult-Endo" = '#6F2DA8', "FCX-adult-ExN-15" = '#00B6EB', "FCX-adult-Olig-3" = '#FDE725FF')

fcx_vln_cols_recode <- c("FCX-adult-ExN-1" = '#00B6EB' , "FCX-adult-ExN-2" = '#00B6EB', "FCX-adult-InN-1" = '#3CBB75FF',
                         "FCX-adult-ExN-3" = '#00B6EB', "FCX-adult-Olig-1" = '#FDE725FF', "FCX-adult-InN-2" = '#3CBB75FF', 
                         "FCX-adult-InN-3" = '#3CBB75FF', "FCX-adult-Ast-1" = '#FF5959', "FCX-adult-InN-4" = '#3CBB75FF', 
                         "FCX-adult-ExN-4" = '#00B6EB', "FCX-adult-ExN-5" = '#00B6EB', "FCX-adult-InN-5" = '#3CBB75FF', 
                         "FCX-adult-InN-6" = '#3CBB75FF', "FCX-adult-Olig-2" = '#FDE725FF', "FCX-adult-ExN-6" = '#00B6EB', 
                         "FCX-adult-ExN-7" = '#00B6EB', "FCX-adult-InN-7" = '#3CBB75FF', "FCX-adult-ExN-8" = '#00B6EB', 
                         "FCX-adult-InN-8" = '#3CBB75FF', "FCX-adult-ExN-9" = '#00B6EB', "FCX-adult-ExN-10" = '#00B6EB', 
                         "FCX-adult-ExN-11" = '#00B6EB', "FCX-adult-MG-1" = '#F58231', "FCX-adult-ExN-12" = '#00B6EB', 
                         "FCX-adult-InN-9" = '#3CBB75FF', "FCX-adult-ExN-13" = '#00B6EB', "FCX-adult-InN-10" = '#3CBB75FF', 
                         "FCX-adult-Ast-2" = '#FF5959', "FCX-adult-InN-11" = '#3CBB75FF', "FCX-adult-InN-12" = '#3CBB75FF', 
                         "FCX-adult-InN-13" = '#3CBB75FF', "FCX-adult-InN-14" = '#3CBB75FF', "FCX-adult-ExN-14" = '#00B6EB', 
                         "FCX-adult-Endo" = '#6F2DA8', "FCX-adult-ExN-15" = '#00B6EB', "FCX-adult-Olig-3" = '#FDE725FF')

}

# Old colours
cer_colours <- c('#3CBB75FF', '#FAA0A0', '#EF0029', '#76B5C5', '#00FF00A5',
                 '#006400', '#B7FFB7', '#D078FF', '#CEE5FD', '#D2042D', 
                 '#00BDD2', '#00B6EB', '#DCBEFF', '#10A53DFF', '#95D840FF', 
                 '#9DC183', '#FDE725FF', '#1E81B0', '#708238', '#9A6324', 
                 '#F58231')




pfc_colours <- c('#DCBEFF', '#9A6324', '#CEE5FD', '#CEE5FD', '#CEE5FD', 
                 '#CEE5FD', '#CEE5FD', '#3CBB75FF', '#3CBB75FF', '#3CBB75FF', 
                 '#3CBB75FF', '#D078FF', '#F58231', '#CCCCCC', '#FDE725FF', 
                 '#FF5959', '#FF5959')