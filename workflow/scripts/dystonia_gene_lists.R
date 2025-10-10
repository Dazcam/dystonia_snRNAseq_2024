#--------------------------------------------------------------------------------------
#
#    Dystonia - gene lists and colour schemes
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

sng_genes <- c('GFAP', 'ORL1', 'GINS3', # Astrocytes (Agarwal)
              'MOG', 'MOBP', 'PALM2', 'LGALS1', 'PPM1G', # ODCs (Agarwal)
              'RGS5', # Endo (Agarwal)
              'CSF1R', # Microglia (Agarwal)
              'VCAN', # OPCs (Agarwal)
              'TH', 'SLC6A3', 'LMX1B', 'KCNJ6', 'NR4A2', # DaNs (Agarwal)
              'GAD1', 'GAD2', 'GABRA1', 'GABRA2', 'RYR2', # InN (Agarwal)
              'AQP4', # Astrocytes (Kamath)
              'SLC17A6', 'RBFOX3', 'SLC18A2', 'GALNTL6', 'RIT2', # Neurons (Kamath)
              'C3', # Microglia (Kamath)
              'DCC', # OPCs (Kamath)
              'FLT1', # Endo (Kamath)
              'PDGFRB', # Pericytes (Kamath)
              'COL1A2') # Fibroblast-like cells (Kamath)

## Str  -------------------------------------------------------------------------------
str_genes <- c('DRD1', 'DRD2', 'TAC1', 'PENK',
               'FOXP1', 'MYT1L', 'MEIS2', 'CSDE1',
               'SOX11', 'BCL11B', 'YBX1', 'EBF1',
               'RARB', 'NR1D1', 'TEF', 'NOS1',
               'SST', 'NPY', 'GAD1', 'GAD2', 'PVALB', 'CALB2', 'VIP',
               'CCK', 'TH', 'CHAT', 'CALB1', 'SOX4', 'DLX2',
               'LHX6', 'LHX8', 'NKX2.1', 'NR2F2')

# See: [Saunders et al. 2018](https://doi.org/10.1016/j.cell.2018.07.028)
small_populations <- c('VTN', 'KCNJ8', 'ABCC9', 'ART3',  # Pericytes
                       'ACTA2', 'RGS5', 'ALDH1A1',      # Smooth muscle cells
                       'SLC38A2', 'SLC4A10', 'SLC26A2', 'SLC47A1', 'FXYD5',  # Fibroblast
                       'ATP1B1', 'COL4A1', 'COL4A2', 'COL15A1', 'COLA1',
                       'COL3A1',
                       'TM4SF1', 'SLC38A5', 'CYTL1', 'BMX', 'MGP',    # Endothelial
                       'FBLN5', 'ELN', 'IGFBP4', 'CLU')

str_final_genes <- factor(c('GAD1', 'GAD2', 'DRD1', 'DRD2', 'CALB2', 
                            'CCK', 'VIP', 'NPY','SST', 'NPY', 
                            'GFAP', 'AQP4', 'FOXJ1', 'OLIG1', 'OLIG2', 
                            'C3'))

str_umap_cols_recode <- c("Str-adult-InN-1" = '#708238', "Str-adult-InN-2" = '#31C53F', "Str-adult-InN-3" = '#3CBB75FF', 
                          "Str-adult-InN-4" = '#9DC183', "Str-adult-InN-5" = '#006400', "Str-adult-InN-6" = '#95D840FF', 
                          "Str-adult-InN-7" = '#2FF18B', "Str-adult-InN-8" = '#B7FFB7', "Str-adult-InN-9" = '#3CBB75FF',
                          "Str-adult-InN-10" = '#3CBB75FF', "Str-adult-Ast" = '#FF5959', "Str-adult-Ependy" = '#c19adf', 
                          "Str-adult-Olig-1" = '#FDE725FF', "Str-adult-OPC" = '#FDE725FF', "Str-adult-MG" = '#F58231',  
                          "Str-adult-Misc-1" = '#CCCCCC')

str_vln_cols_recode <- c("Str-adult-InN-1" = '#3CBB75FF', "Str-adult-InN-2" = '#3CBB75FF', "Str-adult-InN-3" = '#3CBB75FF', 
                         "Str-adult-InN-4" = '#3CBB75FF', "Str-adult-InN-5" = '#3CBB75FF', "Str-adult-InN-6" = '#3CBB75FF', 
                         "Str-adult-InN-7" = '#3CBB75FF', "Str-adult-InN-8" = '#3CBB75FF', "Str-adult-InN-9" = '#3CBB75FF',
                         "Str-adult-InN-10" = '#3CBB75FF', "Str-adult-Ast" = '#FF5959', "Str-adult-Ependy" = '#c19adf',
                         "Str-adult-Olig-1" = '#FDE725FF', "Str-adult-OPC" = '#FDE725FF', "Str-adult-MG" = '#F58231', 
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

cer_final_genes <- factor(c('SLC17A7', 'EOMES', 'RELN', 'GAD1', 'GAD2', 'PVALB', 
                            'GFAP', 'AQP4', 'LINC01727', 'OLIG1', 'OLIG2', 'C3', 
                            'C1QB', 'COL15A1', 'ABCC9', 'RGS5', 'CCL5'))


cer_umap_cols_recode <- c("Cer-adult-ExN" = '#00B6EB', "Cer-adult-UBC" = '#CEE5FD', "Cer-adult-InN-1" = '#95D840FF', 
                          "Cer-adult-InN-2" = '#31C53F', "Cer-adult-InN-3" = '#3CBB75FF', "Cer-adult-BGli" = '#6F2DA8', 
                          "Cer-adult-Olig" = '#FDE725FF', "Cer-adult-OPC" = '#FFBF00', "Cer-adult-Ast" = '#FF5959', 
                          "Cer-adult-MG" = '#F58231', "Cer-adult-Fibro" = '#b57e1d', "Cer-adult-Mural" = '#993300',
                          "Cer-adult-Leuko" = "#F032E6")

cer_vln_cols_recode <- c("Cer-adult-ExN" = '#00B6EB', "Cer-adult-UBC" = '#00B6EB', "Cer-adult-InN-1" = '#3CBB75FF', 
                         "Cer-adult-InN-2" = '#3CBB75FF', "Cer-adult-InN-3" = '#3CBB75FF', "Cer-adult-BGli" = '#6F2DA8', 
                         "Cer-adult-Olig" = '#FDE725FF', "Cer-adult-OPC" = '#FFBF00', "Cer-adult-Ast" = '#FF5959', 
                         "Cer-adult-MG" = '#F58231', "Cer-adult-Fibro" = '#b57e1d', "Cer-adult-Mural" = '#993300',
                         "Cer-adult-Leuko" = "#F032E6")

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
                           'OLIG2', 'GFAP', 'AQP4','C3', 'C1QB', 'CLDN5'))

fcx_umap_cols_recode <- c("FC-adult-ExN-1" = '#00B6EB', "FC-adult-ExN-2" = '#CEE5FD', "FC-adult-ExN-3" = '#598BAF',
                          "FC-adult-ExN-4" = '#6693F5', "FC-adult-ExN-5" = '#89CFEF', "FC-adult-ExN-6" = '#95C8D8', 
                          "FC-adult-ExN-7" = '#4682B4', "FC-adult-ExN-8" = '#73C2FB', "FC-adult-InN-1" = '#708238', 
                          "FC-adult-InN-2" = '#31C53F', "FC-adult-InN-3" = '#3CBB75FF', "FC-adult-InN-4" = '#9DC183', 
                          "FC-adult-InN-5" = '#006400', "FC-adult-InN-6"= '#95D840FF', "FC-adult-InN-7" = '#2FF18B', 
                          "FC-adult-InN-8" = '#B7FFB7', "FC-adult-Olig" = '#FDE725FF', "FC-adult-OPC" = '#FFBF00', 
                          "FC-adult-Ast" = '#FF5959', "FC-adult-MG" = '#F58231', "FC-adult-Endo" = '#9A6324')


fcx_vln_cols_recode <- c("FC-adult-ExN-1" = '#00B6EB', "FC-adult-ExN-2" = '#00B6EB', "FC-adult-ExN-3" = '#00B6EB',
                         "FC-adult-ExN-4" = '#00B6EB', "FC-adult-ExN-5" = '#00B6EB', "FC-adult-ExN-6" = '#00B6EB', 
                         "FC-adult-ExN-7" = '#00B6EB', "FC-adult-ExN-8" = '#00B6EB', "FC-adult-InN-1" = '#3CBB75FF', 
                         "FC-adult-InN-2" = '#3CBB75FF', "FC-adult-InN-3" = '#3CBB75FF', "FC-adult-InN-4" = '#3CBB75FF', 
                         "FC-adult-InN-5" = '#3CBB75FF', "FC-adult-InN-6"= '#3CBB75FF', "FC-adult-InN-7" = '#3CBB75FF', 
                         "FC-adult-InN-8" = '#3CBB75FF', "FC-adult-Olig" = '#FDE725FF', "FC-adult-OPC" = '#FFBF00', 
                         "FC-adult-Ast" = '#FF5959', "FC-adult-MG" = '#F58231', "FC-adult-Endo" = '#9A6324')

### Substantia Nigra 
sng_umap_cols_recode <- c("SNg-adult-ExN" = '#00B6EB', "SNg-adult-InN-1"  = '#3CBB75FF', "SNg-adult-InN-2"  = '#31C53F', 
                          "SNg-adult-InN-3"  = '#006400', "SNg-adult-InN-4"  = '#2FF18B', "SNg-adult-InN-5"  = '#B7FFB7', 
                          "SNg-adult-InN-6"  = '#95D840FF', "SNg-adult-InN-7"  = '#708238', "SNg-adult-Ast" = '#FF5959', 
                          "SNg-adult-Olig-1" = '#FDE725FF', "SNg-adult-Olig-2" = '#FDE725FF', "SNg-adult-OPC" = '#FFBF00' ,
                          "SNg-adult-MG" = '#F58231', "SNg-adult-Leuko" = "#F032E6", "SNg-adult-Endo" = '#9A6324', 
                          "SNg-adult-Fibro" = '#b57e1d', "SNg-adult-Tcell" = '#6F2DA8')

sng_vln_cols_recode <- c("SNg-adult-ExN" = '#00B6EB', "SNg-adult-InN-1"  = '#3CBB75FF', "SNg-adult-InN-2"  = '#3CBB75FF', 
                         "SNg-adult-InN-3"  = '#3CBB75FF', "SNg-adult-InN-4"  = '#3CBB75FF', "SNg-adult-InN-5"  = '#3CBB75FF', 
                         "SNg-adult-InN-6"  = '#3CBB75FF', "SNg-adult-InN-7"  = '#3CBB75FF', "SNg-adult-Ast" = '#FF5959', 
                         "SNg-adult-Olig-1" = '#FDE725FF', "SNg-adult-Olig-2" = '#FDE725FF', "SNg-adult-OPC" = '#FFBF00' ,
                         "SNg-adult-MG" = '#F58231', "SNg-adult-Leuko" = "#F032E6", "SNg-adult-Endo" = '#9A6324', 
                         "SNg-adult-Fibro" = '#b57e1d', "SNg-adult-Tcell" = '#6F2DA8')

sng_final_genes <- c('SLC17A7', 'SLC17A6', 'SLC17A8', 'RBFOX3',
                     'GAD1', 'GAD2', 'SLC18A2', 'SLC6A12', 'SST', 'PVALB', 'CALB2',
                     'GFAP', 'AQP4', 
                     'OLIG1', 'MOG',
                     'VCAN',
                     'C3', 'CSF1R',
                     'FLT1',
                     'COL1A2',
                     "SKAP1")

# Fetal gene sets
fetal_fcx_vln_recode <- c("FC-fetal-ExN-1" = '#00B6EB', "FC-fetal-ExN-2" = '#00B6EB', "FC-fetal-ExN-3" = '#00B6EB', 
                          "FC-fetal-ExN-4" = '#00B6EB', "FC-fetal-ExN-5" = '#00B6EB', "FC-fetal-InN-1" = '#3CBB75FF', 
                          "FC-fetal-InN-2" = '#3CBB75FF', "FC-fetal-InN-3" = '#3CBB75FF', "FC-fetal-InN-4" = '#3CBB75FF', 
                          "FC-fetal-OPC" = '#FDE725FF', "FC-fetal-RG-1" = '#FF5959', "FC-fetal-RG-2" = '#FF5959', 
                          "FC-fetal-MG" = '#F58231', "FC-fetal-CycPro" = '#DCBEFF', "FC-fetal-IP" = '#B200ED', 
                          "FC-fetal-Endo" = '#9A6324', "FC-fetal-N-undef" = '#CCCCCC')

fetal_cer_vln_recode <- c("Cer-fetal-ExN-1" = '#00B6EB', "Cer-fetal-ExN-2" = '#00B6EB', "Cer-fetal-ExN-3" = '#00B6EB', 
                          "Cer-fetal-ExN-4" = '#00B6EB', "Cer-fetal-ExN-5" = '#00B6EB', "Cer-fetal-InN-1" = '#3CBB75FF', 
                          "Cer-fetal-InN-2" = '#3CBB75FF', "Cer-fetal-InN-3" = '#3CBB75FF', "Cer-fetal-InN-4" = '#3CBB75FF', 
                          "Cer-fetal-InN-5" = '#3CBB75FF', "Cer-fetal-InN-6" = '#3CBB75FF', "Cer-fetal-InN-7" = '#3CBB75FF',
                          "Cer-fetal-InN-8" = '#3CBB75FF', "Cer-fetal-OPC" = '#FDE725FF', "Cer-fetal-RG-1" = '#FF5959', 
                          "Cer-fetal-RG-2" = '#FF5959', "Cer-fetal-RG-3" = '#FF5959', "Cer-fetal-MG" = '#F58231',
                          "Cer-fetal-CycPro" = '#DCBEFF',  "Cer-fetal-IP" = '#B200ED', "Cer-fetal-Endo" = '#9A6324')

fetal_ge_vln_recode <- c("GE-fetal-InN-1" = '#3CBB75FF', "GE-fetal-InN-2" = '#3CBB75FF', "GE-fetal-InN-3" = '#3CBB75FF',
                         "GE-fetal-InN-4" = '#3CBB75FF', "GE-fetal-InN-5" = '#3CBB75FF', "GE-fetal-InN-6" = '#3CBB75FF', 
                         "GE-fetal-InN-7" = '#3CBB75FF', "GE-fetal-RG-1" = '#FF5959', "GE-fetal-RG-2" = '#FF5959', 
                         "GE-fetal-RG-3" = '#FF5959', "GE-fetal-CycPro" = '#DCBEFF')

# Colours
greens <- c('#3CBB75FF', '#00FF00A5','#006400', '#B7FFB7', '#10A53DFF',
            '#95D840FF', '#9DC183',  '#708238', '#55C667FF', '#73D055FF',
            '#567D46')

ExN_blues <- c('#76B5C5', '#00BDD2', '#CEE5FD', '#00B6EB', '#ABDBE3',
               '#1E81B0', '#B8D2EB', '#779CBA')

reds <- c('#FAA0A0', '#FF5959', '#F75151', '#EF0029', '#D2042D')

purples <- c('#B200ED',  '#DCBEFF', '#6F2DA8')

# general
fibroblast <- c('SLC38A2', 'SLC4A10', 'SLC26A2', 'SLC47A1', 'FXYD5', 
                'ATP1B1', 'COL4A1', 'COL4A2', 'COL15A1', 'COLA1',
                'COL3A1')

endothelial <- c('TM4SF1', 'SLC38A5', 'CYTL1', 'BMX', 'MGP',
                 'FBLN5', 'ELN', 'IGFBP4', 'CLU')

pericytes <- c('VTN', 'KCNJ8', 'ABCC9', 'ART3')

bergmann <- c('NPY', 'TNC', 'LINC01727', 'FST', 'MT2A', 'PIFO', 'RSPH1')
kozareva <- c('PPP1R17', 'GABRA6', 'EOMES', 'LYPD6', 'PRKCD', 'SORC3', 
              'PTPRK', 'PRKCD', 'NXPH1', 'CDH22', 'KLHL1', 'ALDH1A3', 'SLC6A5', 'HTR2A', 'EDIL3',
              'DCN', 'KCNJ8', 'MRC1', 'FIT1', 'FOXJ1', 'SLC6A5', 'GRM2', 'SST', 'PTPRC')
leuko <- c("PTPRC", "SKAP1", "ARHGAP15", "PRKCH", "IKZF1", "STAT4", "DOCK8", 
           "CD247", "TC2N", "IQGAP2", "FYB1", "SAMD3", "BCL11B", "CARD11", 
           "EMB", "ETS1", "HLA-E", "LCP1", "CD96", "THEMIS", "STK17B", "APBB1IP", 
           "IKZF3", "TNFAIP8", "CLEC2D", "GNG2", "CCL5", "CD53", "FLI1", 
           "ZC3HAV1")

# dput(read_tsv('~/Desktop/dystonia_snRNAseq_2024/results/01R_objects/cer_marker_genes.tsv') %>%
#        filter(cluster == 'Cer-adult-Leuko?') %>%
#        slice_head(n = 30) %>%
#        pull(gene))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
