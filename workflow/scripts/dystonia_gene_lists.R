cer_genes <- c('MFAP4', 'MGP', 'RBFOX3', 'RELN',         # Granule cells
               'NEFM', 'NEUROD1',                        # Ex cer nuclei cell (eCN)
               'PCP4', 'RORA', 'CA8', 'ITPR1',           # Purkinje cells  
               'GAD1', 'GAD2', 'PAX2',                   # GABA Ns
               'MBP', 'APOD', 'OLIG1', 'OLIG2',          # Oligodendrocytes
               'GDF10', 'PAX3',                          # Bergmann cells
               'AQP4', 'SOX9', 'GFAP',                    # Astrocytes 
               'AIF1', 'C1QB',                           # Microglia 
               'EOMES')                                  # Unipolar brush cell (UBC)

fcx_genes <- c('CUX3', #L2/3
               'RORB', #L4
               'TLE4', #L5/6
               'GFAP', 'SLC1A2', #Ast
               'MBP', #Oligs
               "GLI3", "OLIG1",
               "MKI67", "C3", "ITM2A", "SST", "CALB2", 
               "SCGN", "TLE3", "FEZF2", "CRYM", "LHX2"
               )

fcx_genes_cameron <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1",
                       "MKI67", "C3", "ITM2A", "SST", "CALB2", 
                       "SCGN", "TLE3", "FEZF2", "CRYM", "LHX2")
ge_genes_cameron <- c("GAD1", "SLC17A7", "EOMES", "GLI3", "OLIG1", 
                      "MKI67", "C3", "ITM2A", "LHX6", "SIX3", 
                      "PROX1", "TSHZ1", "DLX1", "SCGN")
hip_genes_cameron <- c("NEUROD1", "GRIK4", "EOMES", "GLI3", "OLIG1", 
                       "MKI67", "C3", "ITM2A", "SLC17A6", "ADARB2",
                       "GAD2", "TNC", "PROX1", "RELN", "LHX6")
tha_genes_cameron <- c("EOMES", "GLI3", "OLIG1", "MKI67", "C3", 
                       "ITM2A", "SLC1A6", "LHX9", "TNC", "GAD2", 
                       "PAX6", "SLC17A6")
cer_genes_cameron <- c("GAD1", "EOMES", "GLI3", "OLIG1", "MKI67", 
                       "C3", "ITM2A", "CA8", "ITPR1", "RBFOX3", "RELN")

general_genes <- c('SLC17A7', 'SLC17A6', 'SLC17A8', # VGLUT1-3
                   'SLC6A1', 'SLC6A13', 'SLC6A11', 'SLC6A12', # GABA transporters
                   'SST', 'NPY', 'GAD1', 'GAD2', 'PVALB', 'CALB2', 'VIP', # InN markers
                   'C3', 'C1QB',                           # MG markers
                   'AQP4', 'SOX9', 'GFAP',        # Astrocytes 
                   'OLIG1', 'OLIG2', 'MBP',                # Oligodendrocytes
                   'PDGRFA', 'PMP2',
                   'EOMES', 'EBF1', 'ABCB1') 

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


# Can only be loaded when appropriate object exist
if (exists('seurat_sk_str')) {
  str_clusters_recode <- seurat_sk_str@meta.data %>% 
    as_tibble() %>%
    mutate(harmony_clusters_recode = recode(harmony_clusters, 
                                            `0` = "Str-adult-InN-1", 
                                            `1` = "Str-adult-InN-2",
                                            `2` = "Str-adult-InN-3",
                                            `3` = "Str-adult-InN-4",
                                            `4` = "Str-adult-InN-5",
                                            `5` = "Str-OPC-1",
                                            `6` = "Str-adult-InN-6",
                                            `7` = "Str-adult-InN-7",
                                            `8` = "Str-adult-InN-8",
                                            `9` = "Str-adult-Ast-1",
                                            `10` = "Str-adult-InN-9",
                                            `11` = "Str-adult-InN-10",
                                            `12` = "Str-adult-InN-11",
                                            `13` = "Str-OPC-2",
                                            `14` = "Str-adult-InN-12",
                                            `15` = "Str-adult-InN-13",
                                            `16` = "Str-adult-ExN-1",
                                            `17` = "Str-adult-MG",
                                            `18` = "Str-adult-InN-14",
                                            `19` = "Str-adult-InN-15",
                                            `20` = "Str-adult-ExN-2",
                                            `21` = "Str-adult-InN-16",
                                            `22` = "Str-adult-InN-17",
                                            `23` = "Str-adult-Ast-1",
                                            `24` = "Str-adult-InN-18",
                                            `25` = "Str-OPC-3")) %>%
    pull(harmony_clusters_recode)
  
  str_umap_cols_recode <- c("Str-adult-InN-1" = '#3CBB75FF', "Str-adult-InN-2" = '#31C53F', "Str-adult-InN-3" = '#708238', 
                            "Str-adult-InN-4" = '#B7FFB7', "Str-adult-InN-5" = '#006400', "Str-OPC-1" = '#FDE725FF', 
                            "Str-adult-InN-6" = '#95D840FF', "Str-adult-InN-7" = '#2FF18B', "Str-adult-InN-8" = '#9DC183', 
                            "Str-adult-Ast-1" = '#FF5959', "Str-adult-InN-9" = '#3CBB75FF', "Str-adult-InN-10" = '#31C53F', 
                            "Str-adult-InN-11" = '#708238', "Str-OPC-2" = '#FDE725FF', "Str-adult-InN-12" = '#006400', 
                            "Str-adult-InN-13" = '#95D840FF', "Str-adult-ExN-1" = '#00B6EB', "Str-adult-MG" = '#F58231', 
                            "Str-adult-InN-14" = '#2FF18B', "Str-adult-InN-15" = '#9DC183', "Str-adult-ExN-2" = '#00B6EB', 
                            "Str-adult-InN-16" = '#0CB702', "Str-adult-InN-17" = '#00BE67', "Str-adult-Ast-2" = '#FF5959', 
                            "Str-adult-InN-18" = '#7CAE00', "Str-OPC-3" = '#FDE725FF')
  
  str_vln_cols_recode <- c("Str-adult-InN-1" = '#3CBB75FF', "Str-adult-InN-2" = '#3CBB75FF', "Str-adult-InN-3" = '#3CBB75FF', 
                           "Str-adult-InN-4" = '#3CBB75FF', "Str-adult-InN-5" = '#3CBB75FF', "Str-OPC-1" = '#FDE725FF', 
                           "Str-adult-InN-6" = '#3CBB75FF', "Str-adult-InN-7" = '#3CBB75FF', "Str-adult-InN-8" = '#3CBB75FF', 
                           "Str-adult-Ast-1" = '#FF5959', "Str-adult-InN-9" = '#3CBB75FF', "Str-adult-InN-10" = '#3CBB75FF', 
                           "Str-adult-InN-11" = '#3CBB75FF', "Str-OPC-2" = '#FDE725FF', "Str-adult-InN-12" = '#3CBB75FF', 
                           "Str-adult-InN-13" = '#3CBB75FF', "Str-adult-ExN-1" = '#00B6EB', "Str-adult-MG" = '#F58231', 
                           "Str-adult-InN-14" = '#3CBB75FF', "Str-adult-InN-15" = '#3CBB75FF', "Str-adult-ExN-2" = '#00B6EB', 
                           "Str-adult-InN-16" = '#3CBB75FF', "Str-adult-InN-17" = '#3CBB75FF', "Str-adult-Ast-2" = '#FF5959', 
                           "Str-adult-InN-18" = '#3CBB75FF', "Str-OPC-3" = '#FDE725FF')
}

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