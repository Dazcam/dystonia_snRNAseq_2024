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

cholinergic_genes <- c('CHRNA2', 'CHRNA3', 'CHRNA4', 'CHRNA5', 'CHRNA6',  # Nicotinic (brain)
                       'CHRNA7', 'CHRNA9', 'CHRNA10', 'CHRNB2', 'CHRNB4', # Nicotinic (brain)
                       'CHRM1', 'CHRM2', 'CHRM3', 'CHRM4', 'CHRM5', # Muscarinic (brain)
                       'D1BR', # D5 Dopamine receptor (input for dopamine in Str)
                       'ACHE', # Acetylcholinesterase
                       'SLC18A3', # Vesicular acetylcholine transporter
                       'ZIC4', 'LHX6', 'GBX2', 'FGF8', 'FGF17', 'DBX1') # https://doi.org/10.3389/fnmol.2019.00204 

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

cer_final_genes <- factor(c('GAD1', 'GAD2', 'PVALB', 'SLC17A7',
                            'RBFOX3', 'RELN', 'EOMES', 'PAX3', 'GFAP',
                            'SOX9', 'AQP4', 'OLIG1', 'OLIG2', 'MBP',
                            'C3', 'C1QB', 'EBF1', 'ABCB1', 'APOD'))


cer_umap_cols_recode <- c("Cer-adult-InN-1" = '#95D840FF', "Cer-adult-InN-2" = '#31C53F', "Cer-adult-InN-3" = '#3CBB75FF', 
                          "Cer-adult-ExN" = '#00B6EB', "Cer-adult-UBC" = '#CEE5FD', "Cer-adult-BGli?" = '#6F2DA8', 
                          "Cer-adult-Olig" = '#FDE725FF', "Cer-adult-OPC" = '#FFBF00', "Cer-adult-Ast" = '#FF5959', 
                          "Cer-adult-MG" = '#F58231', "Cer-adult-Endo?" = '#9A6324', "Cer-adult-Pericyte?" = '#CCCCCC',
                          "Cer-adult-Leuko?" = "#F032E6")

cer_vln_cols_recode <- c("Cer-adult-InN-1" = '#3CBB75FF', "Cer-adult-InN-2" = '#3CBB75FF', "Cer-adult-InN-3" = '#3CBB75FF', 
                         "Cer-adult-ExN" = '#00B6EB', "Cer-adult-UBC" = '#00B6EB', "Cer-adult-BGli?" = '#6F2DA8', 
                         "Cer-adult-Olig" = '#FDE725FF', "Cer-adult-OPC" = '#FFBF00', "Cer-adult-Ast" = '#FF5959', 
                         "Cer-adult-MG" = '#F58231', "Cer-adult-Endo?" = '#9A6324', "Cer-adult-Pericyte?" = '#CCCCCC',
                         "Cer-adult-Leuko?" = "#F032E6")

## FCX  ----
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

fcx_final_genes <- factor(c('SLC17A7', 'RORB', 'TLE4', 'GAD1', 'GAD2', 
                            'PVALB', 'SST', 'CALB2','VIP',  'OLIG1', 
                           'OLIG2', 'MBP', 'GFAP', 'SOX9', 'AQP4', 
                            'C3', 'C1QB', 'EBF1', 'ABCB1', 'APOD'))

fcx_umap_cols_recode <- c("Ctx-adult-ExN-1" = '#00B6EB', "Ctx-adult-ExN-2" = '#CEE5FD', "Ctx-adult-ExN-3" = '#598BAF',
                          "Ctx-adult-ExN-4" = '#6693F5', "Ctx-adult-ExN-5" = '#89CFEF', "Ctx-adult-ExN-6" = '#95C8D8', 
                          "Ctx-adult-ExN-7" = '#4682B4', "Ctx-adult-ExN-8" = '#73C2FB', "Ctx-adult-InN-1" = '#708238', 
                          "Ctx-adult-InN-2" = '#31C53F', "Ctx-adult-InN-3" = '#3CBB75FF', "Ctx-adult-InN-4" = '#9DC183', 
                          "Ctx-adult-InN-5" = '#006400', "Ctx-adult-InN-6"= '#95D840FF', "Ctx-adult-InN-7" = '#2FF18B', 
                          "Ctx-adult-InN-8" = '#B7FFB7', "Ctx-adult-Olig" = '#FDE725FF', "Ctx-adult-OPC" = '#FFBF00', 
                          "Ctx-adult-Ast" = '#FF5959', "Ctx-adult-MG" = '#F58231', "Ctx-adult-Undef" = '#CCCCCC')


fcx_vln_cols_recode <- c("Ctx-adult-ExN-1" = '#00B6EB', "Ctx-adult-ExN-2" = '#00B6EB', "Ctx-adult-ExN-3" = '#00B6EB',
                         "Ctx-adult-ExN-4" = '#00B6EB', "Ctx-adult-ExN-5" = '#00B6EB', "Ctx-adult-ExN-6" = '#00B6EB', 
                         "Ctx-adult-ExN-7" = '#00B6EB', "Ctx-adult-ExN-8" = '#00B6EB', "Ctx-adult-InN-1" = '#3CBB75FF', 
                         "Ctx-adult-InN-2" = '#3CBB75FF', "Ctx-adult-InN-3" = '#3CBB75FF', "Ctx-adult-InN-4" = '#3CBB75FF', 
                         "Ctx-adult-InN-5" = '#3CBB75FF', "Ctx-adult-InN-6"= '#3CBB75FF', "Ctx-adult-InN-7" = '#3CBB75FF', 
                         "Ctx-adult-InN-8" = '#3CBB75FF', "Ctx-adult-Olig" = '#FDE725FF', "Ctx-adult-OPC" = '#FFBF00', 
                         "Ctx-adult-Ast" = '#FF5959', "Ctx-adult-MG" = '#F58231', "Ctx-adult-Undef" = '#CCCCCC')
  
# Colours
greens <- c('#3CBB75FF', '#00FF00A5','#006400', '#B7FFB7', '#10A53DFF',
            '#95D840FF', '#9DC183',  '#708238', '#55C667FF', '#73D055FF',
            '#567D46')

ExN_blues <- c('#76B5C5', '#00BDD2', '#CEE5FD', '#00B6EB', '#ABDBE3',
               '#1E81B0', '#B8D2EB', '#779CBA')

reds <- c('#FAA0A0', '#FF5959', '#F75151', '#EF0029', '#D2042D')

reds <- c('#B200ED',  '#DCBEFF', '#6F2DA8')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
