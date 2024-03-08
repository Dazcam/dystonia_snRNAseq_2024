#--------------------------------------------------------------------------------------
#
#    Download public data - Stiletti
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Stiletti data: https://cellxgene.cziscience.com/collections/283d65eb-dd53-496d-adb7-7570c7caa443
#  Note that the Stiletti data has counts stored in seurat_obj[["RNA"]]$data
#  Using Seurat 5 primarily, but also bioconductor packages for QC 

##  Set local to local or remote -----------------------------------------------------
locale <- 'local'

##  Load Packages, functions and variables  -------------------------------------------
source('~/Desktop/dystonia_snRNAseq_2024/workflow/scripts/Renvs.R')
source('~/Desktop/dystonia_snRNAseq_2024/workflow/scripts/dystonia_functions.R')

## Load Data --------------------------------------------------------------------------
# Download dissection data - run once
get_dissection_data(get(paste0(region, '_anns')), anns_table, R_dir, file_format = '.rds')

if (locale == 'remote') { file.create( snakemake@output ) }
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

                    