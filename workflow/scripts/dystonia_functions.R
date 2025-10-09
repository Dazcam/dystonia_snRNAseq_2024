# Load data from Stiletti data repository
#' 
#' Downloads Stiletti data for the region of interest, by dissection, as specified
#' in the anns_table. Note there is currently functionality for only 4 brain regions.
#' 
#' @param dissection_anns A vector of unique regional dissection annotation abbrieviations to download.
#' @param anns_table An excel sheet containing the regional dissection for download.
#' @param out_dir A string specifying the output directory for the downloaded data.
#' @param data_link A string specifying the url where the data are located.
#' @param file_format A string indicating the file format of the data being downloaded `.rds` or `.h5ad`
#' 
#' @examples
#' get_dissection_data(c('CaB', 'Pu'), Stiletti_downloads_table.xlsx, 'results/',
#'                     'https://datasets.cellxgene.cziscience.com/', .rds)
#'                     
get_dissection_data <- function(
    
  dissection_anns = NULL, 
  anns_table = NULL,
  out_dir = NULL, 
  data_link = 'https://datasets.cellxgene.cziscience.com/',
  file_format = NULL
  
) {
  
  options(timeout = max(1000, getOption("timeout"))) # Prevent timeout errors 
  
  dir.create(out_dir)
  stopifnot("dissection_anns must be a vector of Stiletti brain region abbrvs" = is.vector(dissection_anns))
  stopifnot("anns_table must be a Stiletti downloads tibble" = tibble::is_tibble(anns_table))
  
  # Load required libraries for SN handling
  if ('SN' %in% dissection_anns) {
    library(Matrix)
  }
  
  download_ids <- anns_table %>%
    dplyr::filter(abbr %in% dissection_anns) 
  
  message('Downloading data for the following:\n')
  message(paste0(capture.output(download_ids), collapse = "\n"))
  
  message('\nOutdir set to: ', out_dir,' \n')
  
  for (ann in 1:nrow(download_ids)) {
    
    download_id <- download_ids$download_link[ann]
    region_abbr <- download_ids$abbr[ann]
    
    message('Downloading data for: ', region_abbr, '... ')
    
    if (region_abbr == 'SN') {
      
      # Load pre-exported files for SN
      message('Access counts matrix and metadata ... ')
      counts_file <- paste0(out_dir, "SN_counts.mtx")
      cell_meta_file <- paste0(out_dir, "SN_cell_meta.csv")
      gene_meta_file <- paste0(out_dir, "SN_gene_meta.csv")
      
      counts <- readMM(counts_file)
      cellMeta <- read.csv(cell_meta_file)
      geneMeta <- read.csv(gene_meta_file)
      
      # Set the rownames and colnames of matrix
      rownames(counts) <- cellMeta$Barcode
      colnames(counts) <- geneMeta$ensembl_id
      
      message('Creating Seurat object with count matrix ... ')
      seurat_obj <- CreateSeuratObject(counts = t(counts))
      
      # Set the meta data
      message('Adding metadata to Seurat object ... ')
      seurat_obj@meta.data <- cbind(cellMeta, seurat_obj@meta.data)
      rownames(seurat_obj@meta.data) <- colnames(seurat_obj)
      
      # Add dissection for consistency - should already be there
      #seurat_obj@meta.data$dissection <- "SN"
      
      # Remove data frames no longer needed
      rm(counts, geneMeta, cellMeta)
      
      # Save the Seurat object
      message('Saving Seurat object ... ')
      saveRDS(seurat_obj, paste0(out_dir, region_abbr, file_format))
      
      message('Loaded and processed SN data from pre-exported files.')
      
    } else {
      
      download.file(paste0(data_link, download_id, file_format), 
                    destfile = paste0(out_dir, region_abbr, file_format), mode = "wb")
      
      # Repair CBV: too many cols
      if (region_abbr == 'CBV') {
        
        message('Repairing CBV. Too many cols. Metadata column no. for CB should be 32.')
        seurat_obj <- readRDS(paste0(out_dir, region_abbr, file_format))
        seurat_obj@meta.data[1:6] <- NULL
        saveRDS(seurat_obj, paste0(out_dir, region_abbr, file_format))
        
      }
      
    }
    
    # Check that download is what we expected
    seurat_obj <- readRDS(paste0(out_dir, region_abbr, file_format))
    
    meta_obj <- seurat_obj@meta.data %>%
      tibble::as_tibble() %>%
      dplyr::group_by(dissection) %>%
      dplyr::count()
    
    message('Dissections present in Seurat object:\n')
    message(paste0(capture.output(meta_obj), collapse = "\n"), '\n')
    
  }
  
}