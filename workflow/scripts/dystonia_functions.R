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
  
  download_ids <- anns_table %>%
    dplyr::filter(abbr %in% dissection_anns) 
  
  message('Downloading data for the following:\n')
  message(paste0(capture.output(download_ids), collapse = "\n"))
  
  message('\nOutdir set to: ', out_dir,' \n')
  
  for (ann in 1:nrow(download_ids)) {
    
    download_id <- download_ids$download_link[ann]
    region_abbr <- download_ids$abbr[ann]
    
    message('Downloading data for: ', region_abbr, '... ')
    
    download.file(paste0(data_link, download_id, file_format), 
                  destfile = paste0(out_dir, region_abbr, file_format), mode = "wb")
    
    # Repair CBV: too many cols
    if (region_abbr == 'CBV') {
      
      message('Repairing CBV. Too many cols. Metadata column no. for CB should be 32.')
      seurat_obj <- readRDS(paste0(out_dir, region_abbr, file_format))
      seurat_obj@meta.data[1:6] <- NULL
      saveRDS(seurat_obj, paste0(out_dir, region_abbr, file_format))
      
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

# Get a 2 column  the Stiletti cluster annotations
get_cluster_anns <- function(
    
  directory = NULL
  
) {
  
  readxl::read_xlsx(paste0(directory, 'science.add7046_table_s3.xlsx'),
                    sheet = 'Sheet1') %>%
    dplyr::select("Cluster ID", "Class auto-annotation", "Neurotransmitter auto-annotation") %>%
    tidyr::unite(clusterAnn, c("Class auto-annotation", "Neurotransmitter auto-annotation")) %>%
    dplyr::mutate(clusterAnn = str_replace_all(clusterAnn, " ", "_")) %>%
    dplyr::rename(cluster_id = "Cluster ID") %>%
    dplyr::mutate(clusterAnn = str_replace_all(clusterAnn, "_NA", "")) %>%
    dplyr::as_tibble() %>%
    na.omit()
  
}

create_BPCell_seurat_object <- function(
    
  annotations = NULL,
  file_type = '.rds',
  file_dir = NULL
  
){


  message('Creating BPCell object for ', region)
  
  file_set <- annotations
  data_list <- c()
  metadata_list <- c()
  
  for (i in 1:length(file_set)) {
    path <- paste0(file_dir, file_set[i], file_type)
    seurat_obj <- readRDS(path)
    write_matrix_dir(
      mat = seurat_obj[["RNA"]]$data, # Note Stiletti has raw counts in data layer
      dir = paste0(gsub(file_type, "", path), "_BP"),
      overwrite = TRUE)
    
    # Load in BP matrices
    mat <- open_matrix_dir(dir = paste0(gsub(file_type, "", path), "_BP"))
    mat <- Azimuth:::ConvertEnsembleToSymbol(mat = mat, species = "human")
    
    # Get metadata
    metadata_list[[i]] <- seurat_obj[[]]
    data_list[[i]] <- mat
  }
  
  # Name layers
  names(data_list) <- file_set
  
  # Add Metadata
  for (i in 1:length(metadata_list)) {
    metadata_list[[i]]$dataset <- names(data_list)[i]
  }
  
  metadata <- Reduce(rbind, metadata_list)
  
  seurat_merged <- CreateSeuratObject(counts = data_list, meta.data = metadata)
  return(seurat_merged)
  
}

get_meta_col_counts <- function(
    
  seurat_obj = NULL,
  meta_col = NULL
  
  
) {

  message('Creating cell counts for: ', meta_col, '...')  
  seurat_obj@meta.data %>%
    tibble::as_tibble(rownames = 'cell_id') %>%
    dplyr::select(any_of(c(meta_col))) %>%
    dplyr::group_by(.data[[meta_col]]) %>%
    tally() %>%
    mutate(freq = n / sum(n))
  
}

# Adds per cell QCs and identifies cell outliers in sce object
get_cell_outliers <- function(
    
  sce_obj = NULL,
  mad_thresh = 3,
  mad_range = NULL,
  mito_thresh = 5,
  ribo_thresh = 5
  
) {
  
  message('Identifying cell outliers in sce object ...')
  # Pull out Mito / Ribo gene names
  mt_genes <- rownames(sce_obj)[grepl("^MT-", rownames(sce_obj))]
  ribo_genes <- rownames(sce_obj)[grepl("^RP[LS]", rownames(sce_obj))]
  
  # Add per cell QC
  sce_obj <- addPerCellQC(sce_obj, subsets = list(Mito = mt_genes, Ribo = ribo_genes))
  
  # Need log to avoid negative numbers lower threshold
  sum_outlier <- scuttle::isOutlier(sce_obj$sum, nmads = mad_thresh, 
                                    type = mad_range, batch = sce_obj$sample_id, log = TRUE)
  detected_outlier <- scuttle::isOutlier(sce_obj$detected,  nmads = mad_thresh, 
                                         type = mad_range, batch = sce_obj$sample_id, log = TRUE)
  mito_outlier <- sce_obj$subsets_Mito_percent > mito_thresh
  ribo_outlier <- sce_obj$subsets_Ribo_percent > ribo_thresh

  cell_outliers <- sum_outlier | detected_outlier | mito_outlier | ribo_outlier 
  
  outlier_cnts_tbl <- tibble(
    measure = c('umi', 'genes', 'mito', 'ribo', 'total'), 
    count = c(sum(sum_outlier), sum(detected_outlier), sum(mito_outlier),
              sum(ribo_outlier), sum(cell_outliers))  
  ) 
  
  sce_obj$sum_outlier <- sum_outlier
  sce_obj$detected_outlier <- detected_outlier
  sce_obj$mito_outlier <- mito_outlier
  sce_obj$ribo_outlier <- ribo_outlier
  sce_obj$cell_outlier <- cell_outliers
  
  message('Cell numbers that will be excluded at specified thresholds:')
  message(paste0(capture.output(outlier_cnts_tbl), collapse = "\n"), '\n')
  
  # Plot outliers
  message('Plotting ...')
  create_outlier_plots(sce_obj, sum_outlier, detected_outlier, 
                       mito_outlier, ribo_outlier)
  
  return(sce_obj)
  
}

# Create outlier plots, run within `get_cell_outliers` function
# Need to fix the top 50 gene plot
# Should add choice of meta data column for x axis
# Booleans of whether cell is an outlier. T = outlier.
create_outlier_plots <- function(
    
  sce_obj = NULL,
  sum_outlier = NULL,
  detected_outlier = NULL,
  mito_outlier = NULL,
  ribo_outlier = NULL
  
) {
  
  message('Creating outlier plots ', region, ' ... ')

  # Generate main outlier plots
  umi_plot <- scater::plotColData(sce_obj, x = "sample_id", y = "sum",
                                  colour_by = I(sum_outlier)) + ggtitle('UMI per cell')
  gene_plot <- scater::plotColData(sce_obj, x = "sample_id", y = "detected",
                                   colour_by = I(detected_outlier)) + ggtitle('Genes per cell')
  mito_plot <- scater::plotColData(sce_obj, x = "sample_id", y = "subsets_Mito_percent",
                                   colour_by = I(mito_outlier)) + ggtitle('Mito genes per cell')
  ribo_plot <- scater::plotColData(sce_obj, x = "sample_id", y = "subsets_Ribo_percent",
                                   colour_by = I(ribo_outlier)) + ggtitle('Ribo genes per cell')
  
  # Get relative expression of each gene per cell and plot
  # C.sce = counts(sce_obj)
  # C.sce@x = C.sce@x/rep.int(colSums(C.sce), diff(C.sce@p))
  # most_expressed <- order(Matrix::rowSums(C.sce), decreasing = T)[50:1]
  # top_50_plot <- boxplot(as.matrix(t(C.sce[most_expressed, ])), cex = 0.05, 
  #                        las = 1, xlab = "% total count per cell", cex.axis=0.5,
  #                        col = (scales::hue_pal())(50)[50:1], horizontal = TRUE)
  # 
  outliers_plot <- cowplot::plot_grid(umi_plot, gene_plot, mito_plot, ribo_plot)
  plot(outliers_plot) 

}

# Cell_col is bool vector to append to seurat obj - to ID cells for removal
# Could streamline this by collating all genes for removal in one var
subset_seurat_object <- function(
    
  seurat_obj = NULL,
  cell_outliers = NULL,
  genesExpInCell_thresh = 3
  
) {

  message('Subsetting Seurat object: ', region, ' ...\n')
  seurat_obj$cell_outlier <- cell_outliers
  
  # Check seurat obj dimensions before filters
  message('Before filtering:')
  message(paste0(capture.output(seurat_obj), collapse = "\n"), '\n')
  
  # Keep cells that are not outliers
  message('Filtering ', sum(seurat_obj$cell_outlier), ' outlier cells ...')
  seurat_obj <- subset(x = seurat_obj, subset = cell_outlier == FALSE)
  
  # Remove MT genes and MALAT1
  genes_rm <- c(rownames(seurat_obj)[grepl("^MT-", rownames(seurat_obj))], "MALAT1")
  genes_in_many_cells <- rowSums(seurat_obj[["RNA"]]$counts) > genesExpInCell_thresh
  message('Filtering MT genes and MALAT1, ', length(genes_rm), ' genes will be removed ...\n')
  seurat_obj <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), genes_rm))
  
  # Remove genes expresssed in fewer than 3 cells
  genes_in_many_cells <- rowSums(seurat_obj[["RNA"]]$counts) > genesExpInCell_thresh
  message('Filtering genes expressed in fewer than 3 cells, removing a further ', 
          sum(!genes_in_many_cells),  ' genes...\n')
  seurat_obj <- subset(seurat_obj, features = names(genes_in_many_cells[genes_in_many_cells]))
  
  # Check seurat obj dimensions after filters
  message('After filtering:')
  message(paste0(capture.output(seurat_obj), collapse = "\n"), '\n')
  
  return(seurat_obj)
  
}

# Assess hard thresholds for seurat objects / combine with get_cell_outliers?
# Produces plots and counts for thresholds
# Need to get it working per batch
# Probs better doing batch specific thresholds maybe 3 X MAD??
get_hard_thresholds <- function(
    
  seurat_obj = NULL,
  exclusion_vec = c(5, 5, 1000, 5000, 3)
  
) {
  
  exclusion_tbl <- tibble(
    Criterion = c('Mito', 'Ribo', 'Low gene thresh', 'High gene thresh', 'No. genes exp in cells'),
    Threshold = c(paste0('> ', exclusion_vec[1], '% per cell'), 
                  paste0('> ', exclusion_vec[2], '% per cell'), 
                  paste('< ', exclusion_vec[3], ' per cell'), 
                  paste('< ', exclusion_vec[4], ' per cell'), 
                  paste('< ', exclusion_vec[5], ' per cell')),
    'Type excluded' = c('Cell', 'Cell', 'Cell', 'Cell', 'Gene'),
    'No Excluded' = c(sum(seurat_obj@meta.data$percent_mito > exclusion_vec[1]), 
                      sum(seurat_obj@meta.data$percent_ribo > exclusion_vec[2]),
                      sum(seurat_obj@meta.data$nFeature_RNA < exclusion_vec[3]),
                      sum(seurat_obj@meta.data$nFeature_RNA > exclusion_vec[4]),
                      sum(rowSums(seurat_obj[["RNA"]]$counts) < exclusion_vec[5]))
  )
  
  message('Cell / Gene numbers that will be excluded at specified thresholds:')
  message(paste0(capture.output(exclusion_tbl), collapse = "\n"), '\n')
  QC_Plots_Combined_Vln(seurat_object = seurat_obj, 
                        feature_cutoffs = c(exclusion_vec[3], exclusion_vec[4]), 
                        UMI_cutoffs = c(1200, 45000), 
                        mito_cutoffs = exclusion_vec[1], pt.size = 0.1)
}

#' Run the basic Seurat process
#' 
#' Run the basic Seurat pipeline using either log or SCT normalisation.
#' 
#' Note: The option 'log_sketch' runs the log normalisation Seurat pipeline 
#' on a sketch object, but omits the normalisation step, as log normalisation 
#' is run on the raw object just prior to the creation of the sketch object.
#' See create_sketch_object() and https://satijalab.org/seurat/articles/seurat5_sketch_analysis
#' 
#' @param seurat_obj An uncorrected Seurat object.
#' @param norm_method A method for mornalisation either 'log', 'log_sketch', or sct.
#' @param dims Number of principal components (input for run_seurat_process())
#' @param resolution A set of resolution params for clustering (input for run_seurat_process())
#' 
#' @returns A Seurat object
#' 
#' @examples
#' run_seurat_process(seurat_small, 'log', 30, c(0.3, 0.5))
#' 
run_seurat_process <- function( 
    
  seurat_obj = NULL,
  norm_method = 'log',
  dims = 30, 
  resolution = c(0.3, 0.5, 0.8)
  
){
  
  
  if (norm_method == 'log') {
    
    message('\nRunning Seurat pipeline norm method: ', norm_method, '\n')
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", 
                                scale.factor = 10000) %>%
      Seurat::FindVariableFeatures(seurat_obj, verbose = F) %>%
      Seurat::ScaleData(verbose = F) %>%
      Seurat::RunPCA(verbose = F) %>%
      Seurat::FindNeighbors(dims = 1:dims) %>%
      Seurat::FindClusters(resolution = resolution) %>%
      Seurat::RunUMAP(dims = 1:dims)}
  
  
  if (norm_method == 'log_sketch') {
    
    message('\nRunning Seurat pipeline norm method: ', norm_method, '\n')
    seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, verbose = F) %>%
      Seurat::ScaleData(verbose = F) %>%
      Seurat::RunPCA(verbose = F) %>%
      Seurat::FindNeighbors(dims = 1:dims) %>%
      Seurat::FindClusters(resolution = resolution) %>%
      Seurat::RunUMAP(dims = 1:dims)}
  
  if (norm_method == 'sct') {
    
    message('\nRunning Seurat pipeline norm method: ', norm_method, '\n')
    Seurat::SCTransform(seurat_obj, conserve.memory = T, verbose = F) %>%
      Seurat::RunPCA(verbose = F) %>%
      Seurat::RunUMAP(dims = 1:dims, verbose = F) %>%
      Seurat::FindNeighbors(dims = 1:dims, verbose = F) %>%
      Seurat::FindClusters(resolution = resolution, verbose = F)}
  
  return(seurat_obj)
  
}

#' Create a Seurat sketch object
#' 
#' Use this function on an uncorrected Seurat object. Note that log normalisation 
#' is carried out on the full object first, and then run_seurat_process() is run,
#' though normalisation is only run once (hence setting norm_method = 'log_sketch')
#' 
#' Note: that it is currently not recommended to Seurat::SketchData with SCT 
#' normalisation. See [issue 7336](https://github.com/satijalab/seurat/issues/7336).
#' 
#' @param seurat_obj An uncorrected Seurat object.
#' @param norm_method Set to log_sketch (input for run_seurat_process())
#' @param dims Number of principal components (input for run_seurat_process())
#' @param resolution A set of resolution params for clustering (input for run_seurat_process())
#' @param cell_num A positive integer indicating the no. of cells per sample to extract for sketch obj.
#' 
#' @examples
#' create_sketch_object(seurat_small, 'log_sketch', 30, c(0.3, 0.5))
create_sketch_object <- function(
    
  seurat_obj = NULL, 
  norm_method = 'log_sketch',
  dims = 30, 
  resolution = c(0.3, 0.5, 0.8),
  cell_num = 5000
  
) {
  
  
  message('Creating Seurat sketch object: ', region, ' ...')
  seurat_obj <- Seurat::NormalizeData(seurat_obj) %>%
    FindVariableFeatures(verbose = FALSE)
  
  seurat_sketch <- Seurat::SketchData(object = seurat_obj, 
                                      ncells = cell_num, 
                                      method = "LeverageScore", 
                                      sketched.assay = "sketch")
  
  # Take a small subset of data to 
  Seurat::DefaultAssay(seurat_sketch) <- "sketch"
  
  # Run Seurat pipeline
  seurat_sketch <- run_seurat_process(seurat_sketch, 
                                      norm_method = norm_method,
                                      dims = dims,
                                      resolution = resolution)
  
  return(seurat_sketch)
  
}

#' Create a plot list for a set of resolution params
#' 
#' Creates 3 plots for each resolution param. A cluster plot with default cluster IDs.
#' A stacked violin plot and a cluster plot with the braod stiletti annotations. 
#' Only has functionality for sketch object atm. 
#' 
#' @param seurat_obj An uncorrected Seurat object.
#' @param resolution A set of resolution params from a sketch object to plot
#' @param meta_id A vector of cluster ids for each cell, i.e. a col from Seurat object metadata.
#' 
#' @returns A list of plots 
#' 
#' @examples
#' create_resolution_plotlist(seurat_small, c(0.3, 0.5))
create_resolution_plotlist <- function(
    
  seurat_obj = NULL,
  resolution = c(0.3, 0.5, 0.8),
  meta_id = NULL
  
) {
  
  plot_list <- list()
  
  for (res_level in resolution) {
    
    message('Creating plots for res level:', res_level ,'...')
    res_plot <- DimPlot(seurat_obj, group.by = paste0('sketch_snn_res.', res_level), 
                        label = T) + NoLegend() 
    stiletti_plot <- DimPlot(seurat_obj, group.by = 'cell_type', pt.size = 1, raster = F) 
    meta_plot <- DimPlot(seurat_obj, group.by = meta_id,
                         label = T, repel = T, pt.size = 1, raster = F) +  
      Seurat::NoLegend()
    bar_plot <- create_proportion_barplot(seurat_obj, paste0('sketch_snn_res.', res_level), meta_id)
    plot_list[[paste0('res_', res_level)]]  <- res_plot
    plot_list[[paste0('vln_', res_level)]] <- create_stacked_vln_plot(seurat_obj, 
                                                                      paste0('sketch_snn_res.', res_level), 
                                                                      general_genes, 
                                                                      paste0(region, ' res. ', res_level))
    plot_list[[paste0('stiletti_', res_level)]] <- scCustomize::DimPlot_scCustom(seurat_obj, group.by = 'cell_type',
                                                                                 label = T, repel = T,
                                                                                 pt.size = 1, raster = F) + 
      Seurat::NoLegend()
    plot_list[[paste0('meta_id_', res_level)]]  <- meta_plot
    plot_list[[paste0('bar_', res_level)]]  <- bar_plot
    
  }
  
  return(plot_list)
  
}

# Basic QC plot
create_basic_qc_plots <- function(seurat_obj = NULL,
                                  point_size = 0,
                                  meta_id = NULL) {
  
  message('Creating basic qc plots: ', region, ' ...')
  Idents(object = seurat_obj) <- meta_id
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj$complexity <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  vln_plot <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", 
                                               "fraction_mitochondrial", "complexity"), 
                      ncol = 4, pt.size = 0)
  fs_plot_1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  fs_plot_2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  fs_plot <- cowplot::plot_grid(fs_plot_1, fs_plot_2)
  big_plot <- cowplot::plot_grid(vln_plot, fs_plot, ncol = 1)
  
  return(big_plot)
  
}

create_cluster_qc_plot <- function(
    
  seurat_obj = NULL, 
  dims = 30, 
  meta_id = 'sample_id'
  
) {
  
  message('Creating cluster qc plots: ', region, ' ...')
  cluster_plot <- Seurat::DimPlot(seurat_obj, reduction = "umap")
  elbow_plot <- Seurat::ElbowPlot(seurat_obj, ndims = dims)
  dataset_plot <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = 'dataset') 
  donor_plot <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = meta_id)
  bar_plot_dataset <- create_proportion_barplot(seurat_obj, paste0('seurat_clusters'), meta_id = 'dataset')
  bar_plot_donor <- create_proportion_barplot(seurat_obj, paste0('seurat_clusters'), meta_id = meta_id)
  
  qc_plot <- cowplot::plot_grid(cluster_plot, elbow_plot, 
                                dataset_plot, bar_plot_dataset,
                                donor_plot, bar_plot_donor, ncol = 2)
  
}

#' Run integration analysis on a Seurat object
#' 
#' This has functionality to run one or more batch correction algorithims. 
#' Note that CCA takes a very long time to run compared to the others. 
#' Atm there is only functionality for a single resolution param as it's 
#' not possible to set the resolution and cluster.name params in 
#' Seurat::FindClusters.
#' 
#' 
#' @param seurat_obj An uncorrected Seurat object.
#' @param reductions An vector of integration method to run eith 'harmony', 'cca', 'rpca', 'fastmnn'.
#' @param dimensions Number of principal components.
#' @param resolution A single resolution param to send to findClusters().
#' 
#' @returns A Seurat object.
#' 
#' @examples
#' run_integration(seurat_small, 'harmony', 30, 0.1)
#' 
run_integration <- function(
    
  seurat_obj = NULL,
  reductions = 'harmony',
  dimensions = 30, 
  resolution = 0.8
  
  ) {
  
  if ('harmony' %in% reductions) { 
  message('Running Harmony integration ...')
  seurat_obj <- IntegrateLayers(
    object = seurat_obj,
    orig.reduction = "pca", 
    new.reduction = "harmony",
    method = HarmonyIntegration,
    #normalization.method = "SCT",
    verbose = TRUE
  )
  
  message('Finding Harmony clusters and generating UMAP ...')
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dimensions, reduction = "harmony") %>%
    RunUMAP(reduction = "harmony", dims = 1:dimensions, reduction.name = "umap.harmony")
  
  for (res_level in resolution) {
    message('Finding harmony clusters for res level: ', res_level, '...')
    seurat_obj <- FindClusters(seurat_obj, 
                               cluster.name = paste0("harmony_clusters_", res_level), 
                               resolution = res_level)}}
  
  if ('cca' %in% reductions) { 
  message('Running CCA integration ...')
  seurat_obj <- IntegrateLayers(
    object = seurat_obj,
    orig.reduction = "pca", 
    new.reduction = "cca",
    method = CCAIntegration,
    #normalization.method = "SCT",
    verbose = TRUE
  )
  
  message('Finding CCA clusters and generating UMAP ...')
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dimensions, reduction = "cca") %>%
    FindClusters(cluster.name = "cca_clusters", resolution = resolution) %>%
    RunUMAP(reduction = "cca", dims = 1:dimensions, reduction.name = "umap.cca")
  
  }
  
  if ('rpca' %in% reductions) { 
  message('Running RPCA integration ...')
  seurat_obj <- IntegrateLayers(
    object = seurat_obj,
    orig.reduction = "pca", 
    new.reduction = "rpca",
    method = RPCAIntegration,
    #normalization.method = "SCT",
    verbose = TRUE
  )
  
  message('Finding RPCA clusters and generating UMAP ...')
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dimensions, reduction = "rpca") %>%
    FindClusters(cluster.name = "rpca_clusters", resolution = resolution) %>%
    RunUMAP(reduction = "rpca", dims = 1:dimensions, reduction.name = "umap.rpca")
  
  }
  
  if ('fastmnn' %in% reductions) { 
  message('Running FastMNN integration ...')
  seurat_obj <- IntegrateLayers(
    object = seurat_obj,
    orig.reduction = "pca", 
    new.reduction = "fastmnn",
    method = FastMNNIntegration,
    #normalization.method = "SCT",
    verbose = TRUE
  )
  
  message('Finding FastMNN clusters and generating UMAP ...')
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dimensions, reduction = "fastmnn") %>%
    FindClusters(cluster.name = "fastmnn_clusters", resolution = resolution) %>%
    RunUMAP(reduction = "fastmnn", dims = 1:dimensions, reduction.name = "umap.fastmnn")
  
  }
  
  return(seurat_obj)
  
}

#' Plot clusters after batch correction integration
#' 
#' This plots 4 plots, 
#' 
#' @param seurat_obj An uncorrected Seurat object.
#' @param reductions An vector of integration method to run eith 'harmony', 'cca', 'rpca', 'fastmnn'.
#' @param meta_id A vector of cluster ids for each cell, i.e. a col from Seurat object metadata.
#' @param dimensions Number of principal components.
#' @param resolution A single resolution param to send to findClusters().
#' 
#' @returns A ggplot plot list.
#' 
#' @examples
#' create_integration_plot(seurat_small, 'harmony', 'sample_id')
#' 
create_integration_plotlist <- function(
    
  seurat_obj = NULL, 
  algorithm = c('harmony'), 
  meta_id = 'sample_id',
  dims = 30,
  reduction = c(0.3, 0.5, 0.8)
  
) {

  message('Creating integration plot: ', region, ' ...')
  
  plot_list <- list()
  
  for (i in 1:length(algorithm)) {
    
    for (j in 1:length(reduction)) {
      
      message('Creating plots for res level:', j ,'...')
    
    cluster_plot <- DimPlot(seurat_obj, reduction = paste0("umap.", algorithm[i]), 
                            group.by = paste0(algorithm[i], '_clusters_', reduction[j]),
                            label = T, repel = T, pt.size = 1, raster = FALSE) +
      NoLegend()
    vln_plot <- create_stacked_vln_plot(seurat_object, 
                                        paste0(algorithm[i], '_clusters_', reduction[j]), 
                                        general_genes, 
                                        paste0(region, ' ', algorithm[i], ' clusters ', reduction[j]))
    meta_plot <- DimPlot(seurat_obj, reduction = paste0("umap.", algorithm[i]), group.by = meta_id, 
                         label = T, repel = T, pt.size = 1, raster = FALSE) +
      NoLegend()
    bar_plot <- create_proportion_barplot(seurat_obj, paste0(algorithm[i], '_clusters_', reduction[j]), meta_id)
    stiletti_plot <- scCustomize::DimPlot_scCustom(seurat_object, 
                                                   reduction = paste0("umap.", algorithm[i]), 
                                                   group.by = 'cell_type', label = T, repel = T,
                                                   pt.size = 1, raster = FALSE) +
      NoLegend()
    
    plot_list[[paste0(algorithm[i], '_cluster_', reduction[j])]]  <- cluster_plot
    plot_list[[paste0(algorithm[i], '_vln_', reduction[j])]]  <- vln_plot
    plot_list[[paste0(algorithm[i], '_stiletti_', reduction[j])]]  <- stiletti_plot
    plot_list[[paste0(algorithm[i], '_meta_', reduction[j])]]  <- meta_plot
    plot_list[[paste0(algorithm[i], '_bar_', reduction[j])]]  <- bar_plot
    
    }
    
  } 
  
  return(plot_list)
  
}


create_proportion_barplot <- function(seurat_obj = NULL, 
                                      cluster_id = NULL,
                                      meta_id = NULL) {
  
  message('Creating proportion barplot: ', region, ' ...')
  # Implementing this using dplyr: 
  # prop.table(table(seurat_sk_str$seurat_clusters, seurat_sk_str$dataset), margin = 2)
  
  bar_plot <- seurat_obj@meta.data %>%
    as_tibble(rownames = 'cell_id') %>%
    dplyr::select(any_of(c(cluster_id, meta_id))) %>%
    dplyr::count(.data[[cluster_id]], .data[[meta_id]]) %>%
    group_by(.data[[cluster_id]]) %>%          # now required with changes to dplyr::count()
    mutate(prop = prop.table(n)) %>%
    ggplot(aes(fill = .data[[meta_id]], y = prop, x = .data[[cluster_id]])) + 
    geom_bar(position = "fill", stat = "identity") +
    theme_minimal()
  
  return(bar_plot)
  
}

create_stacked_vln_plot <- function(
    
  seurat_obj = NULL,
  set_ident = 'seurat_clusters',
  genes = NULL,
  plot_title = NULL,
  col_pal = c("#FF7373", "#CC4E3D", "#993617", "#FF5C00", "#CC925C", "#996E2E", 
              "#FFC226", "#CCAB00", "#999645", "#F1FF4D", "#A9CC1F", "#689900", 
              "#C1FF73", "#7CCC3D", "#419917", "#33FF00", "#65CC5C", "#2E9932", 
              "#26FF49", "#00CC39", "#459967", "#4DFFA9", "#1FCC8E", "#009974", 
              "#73FFEE", "#3DCCCC", "#178999", "#00C2FF", "#5CA4CC", "#2E6699", 
              "#267DFF", "#0039CC", "#455299", "#4D54FF", "#2C1FCC", "#1F0099", 
              "#A073FF", "#7C3DCC", "#601799", "#AD00FF", "#B65CCC", "#902E99", 
              "#FF26F6", "#CC00AB", "#994581", "#FF4DB8", "#CC1F72", "#990037", 
              "#FF7394", "#CC3D4E")
  
) {
  
  message('Creating stacked vln plot: ', region, ' ...')
  Idents(seurat_obj) <- unname(unlist((seurat_obj[[set_ident]])))
  VlnPlot(seurat_obj, genes, stack = TRUE, flip = TRUE, 
          same.y.lims = TRUE, fill.by = 'ident', cols = col_pal) +
    theme(legend.position = "none",
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 18),
          axis.text.x  = element_text(colour = "#000000", size = 16),
          axis.text.y  = element_text(colour = "#000000", size = 16)) +
    xlab('Cell type') +
    ggtitle(plot_title) +
    facet_grid(feature ~ ., drop = F)
  
}


#' Recode cluster ids in a Seurat Object
#' 
#' @param seurat_obj An Seurat object.
#' @param region The string inidcating a brain region of interest. Currently 'str', 'cer', or 'fcx'.
#' @param meta_id A vector of cluster ids for each cell, i.e. a col from Seurat object metadata.
#' 
#' @returns A factor of cell IDs to attach to bind to Seurat metadata
#' 
#' @examples
#' recode_cluster_ids (seurat_small, 'str', 'harmony_clusters_0.1')
#' 
recode_cluster_ids <- function(
    
  seurat_obj = NULL,
  region = NULL,
  meta_id = NULL
  
) {
  
  message('\nRecoding Cluster IDs for, ', region, ' ...\n')
  
  if (region == 'str') {
    
    str_levels <- c("Str-adult-InN-1", "Str-adult-InN-2", "Str-adult-InN-3", 
                    "Str-adult-InN-4", "Str-adult-InN-5", "Str-adult-InN-6", 
                    "Str-adult-InN-7", "Str-adult-InN-8", "Str-adult-InN-9", 
                    "Str-adult-ExN", "Str-adult-Ast-1", "Str-adult-Ast-2", 
                    "Str-adult-Olig-1", "Str-adult-OPC", "Str-adult-MG", 
                    "Str-adult-Misc-1")
    
    clusters_recode <- seurat_obj@meta.data %>% 
      tibble::as_tibble() %>%
      dplyr::mutate(clust_recode = recode(.data[[meta_id]], 
                                          `0` = "Str-adult-InN-1", 
                                          `1` = "Str-adult-InN-2",
                                          `2` = "Str-adult-Olig-1",
                                          `3` = "Str-adult-InN-5",
                                          `4` = "Str-adult-InN-6",
                                          `5` = "Str-adult-InN-7",
                                          `6` = "Str-adult-Ast-2",
                                          `7` = "Str-adult-OPC",
                                          `8` = "Str-adult-ExN",
                                          `9` = "Str-adult-InN-8",
                                          `10` = "Str-adult-InN-3",
                                          `11` = "Str-adult-InN-4",
                                          `12` = "Str-adult-MG",
                                          `13` = "Str-adult-Misc-1",
                                          `14` = "Str-adult-Ast-1",
                                          `15` = "Str-adult-InN-9")) %>%
      pull(clust_recode)}
  
  if (region == 'cer') {
    
    cer_levels <- c("Cer-adult-ExN", "Cer-adult-UBC", "Cer-adult-InN-1", 
                    "Cer-adult-InN-2", "Cer-adult-InN-3", "Cer-adult-BGli?", 
                    "Cer-adult-Olig", "Cer-adult-OPC", "Cer-adult-Ast", 
                    "Cer-adult-MG", "Cer-adult-Endo?", "Cer-adult-Pericyte?",
                    "Cer-adult-Leuko?")
    
    clusters_recode <- seurat_obj@meta.data %>% 
      tibble::as_tibble() %>%
      dplyr::mutate(clust_recode = recode(.data[[meta_id]], 
                                          `0` = "Cer-adult-ExN", 
                                          `1` = "Cer-adult-Olig",
                                          `2` = "Cer-adult-BGli?",
                                          `3` = "Cer-adult-OPC",
                                          `4` = "Cer-adult-InN-1",
                                          `5` = "Cer-adult-MG",
                                          `6` = "Cer-adult-Ast",
                                          `7` = "Cer-adult-Endo?",
                                          `8` = "Cer-adult-UBC",
                                          `9` = "Cer-adult-Pericyte?",
                                          `10` = "Cer-adult-InN-2",
                                          `11` = "Cer-adult-Leuko?",
                                          `12` = "Cer-adult-InN-3")) %>%
      pull(clust_recode)}
  
  if (region == 'fcx') {
    
    fcx_levels <- c("FC-adult-ExN-1", "FC-adult-ExN-2", "FC-adult-ExN-3", 
                    "FC-adult-ExN-4", "FC-adult-ExN-5", "FC-adult-ExN-6", 
                    "FC-adult-ExN-7", "FC-adult-ExN-8", "FC-adult-InN-1", 
                    "FC-adult-InN-2", "FC-adult-InN-3", "FC-adult-InN-4", 
                    "FC-adult-InN-5", "FC-adult-InN-6", "FC-adult-InN-7", 
                    "FC-adult-InN-8", "FC-adult-Olig", "FC-adult-OPC", 
                    "FC-adult-Ast", "FC-adult-MG", "FC-adult-Undef")
    
    clusters_recode <- seurat_obj@meta.data %>% 
      tibble::as_tibble() %>% 
      dplyr::mutate(clust_recode = recode(.data[[meta_id]], 
                                          `0` = "FC-adult-ExN-1", 
                                          `1` = "FC-adult-InN-1",
                                          `2` = "FC-adult-InN-2",
                                          `3` = "FC-adult-InN-3",
                                          `4` = "FC-adult-ExN-2",
                                          `5` = "FC-adult-ExN-3",
                                          `6` = "FC-adult-ExN-4",
                                          `7` = "FC-adult-InN-4",
                                          `8` = "FC-adult-ExN-5",
                                          `9` = "FC-adult-Ast",
                                          `10` = "FC-adult-Olig",
                                          `11` = "FC-adult-ExN-6",
                                          `12` = "FC-adult-InN-5",
                                          `13` = "FC-adult-ExN-7",
                                          `14` = "FC-adult-InN-6",
                                          `15` = "FC-adult-OPC",
                                          `16` = "FC-adult-InN-7",
                                          `17` = "FC-adult-ExN-8",
                                          `18` = "FC-adult-MG",
                                          `19` = "FC-adult-InN-8",
                                          `20` = "FC-adult-Undef")) %>%
      pull(clust_recode)}
  
  clusters_recode <- factor(clusters_recode, get(paste0(region, '_levels')))
  
  return(clusters_recode)
  
}


#' Plot a paired plot, a umap and a violin plot
#' 
#' @param seurat_obj A Seurat object.
#' @param umap_reduct A string stating the umap reduction you would like to plot.
#' @param meta_id A vector of cluster ids for each cell, i.e. a col from Seurat object metadata.
#' @param genes A vector, or factor, of genes to plot in violin plot.
#' @param col_pal_umap A palette of colours for umap plot.
#' @param col_pal_vln A palette of colours for violin plot.
#' 
#' @returns A patchwork object of a umap and violin plot.
#' 
#' @examples
#' plot_paired_umap_vln(seurat_small, 'harmony', harmony_clusters_0.1', c('GAD1', 'GAD2'))
#' 
plot_paired_umap_vln <- function(
    
  seurat_obj = NULL,
  umap_reduct = NULL,
  meta_id = NULL,
  genes = NULL,
  col_pal_umap = c("#FF7373", "#CC4E3D", "#993617", "#FF5C00", "#CC925C", "#996E2E", 
                   "#FFC226", "#CCAB00", "#999645", "#F1FF4D", "#A9CC1F", "#689900", 
                   "#C1FF73", "#7CCC3D", "#419917", "#33FF00", "#65CC5C", "#2E9932", 
                   "#26FF49", "#00CC39", "#459967", "#4DFFA9", "#1FCC8E", "#009974", 
                   "#73FFEE", "#3DCCCC", "#178999", "#00C2FF", "#5CA4CC", "#2E6699", 
                   "#267DFF", "#0039CC", "#455299", "#4D54FF", "#2C1FCC", "#1F0099", 
                   "#A073FF", "#7C3DCC", "#601799", "#AD00FF", "#B65CCC", "#902E99", 
                   "#FF26F6", "#CC00AB", "#994581", "#FF4DB8", "#CC1F72", "#990037", 
                   "#FF7394", "#CC3D4E"),
  col_pal_vln = c("#FF7373", "#CC4E3D", "#993617", "#FF5C00", "#CC925C", "#996E2E", 
                  "#FFC226", "#CCAB00", "#999645", "#F1FF4D", "#A9CC1F", "#689900", 
                  "#C1FF73", "#7CCC3D", "#419917", "#33FF00", "#65CC5C", "#2E9932", 
                  "#26FF49", "#00CC39", "#459967", "#4DFFA9", "#1FCC8E", "#009974", 
                  "#73FFEE", "#3DCCCC", "#178999", "#00C2FF", "#5CA4CC", "#2E6699", 
                  "#267DFF", "#0039CC", "#455299", "#4D54FF", "#2C1FCC", "#1F0099", 
                  "#A073FF", "#7C3DCC", "#601799", "#AD00FF", "#B65CCC", "#902E99", 
                  "#FF26F6", "#CC00AB", "#994581", "#FF4DB8", "#CC1F72", "#990037", 
                  "#FF7394", "#CC3D4E")
  
  
) {
  
  message('Creating UMAP and Vln Paired plot:\n')
  umap <- Seurat::DimPlot(seurat_obj, group.by = meta_id, reduction = umap_reduct,
                          cols = col_pal_umap, pt.size = 0.1,
                          label = T, label.size = 3, label.box = T, repel = T,
                          raster = F) + 
    NoLegend()
  vln <- create_stacked_vln_plot(seurat_obj, meta_id, genes,
                                 toupper(region), col_pal_vln)
  
  
  umap | vln
  
} 

#' Plot 2 stacked violin plots to compare genes expressed
#' 
#' @param seurat_obj An Seurat object.
#' @param meta_id A vector of cluster ids for each cell, i.e. a col from Seurat object metadata.
#' @param genes_a A vector, or factor, of genes to plot in plot a.
#' @param genes_a A vector, or factor, of genes to plot in plot b.
#' @param col_pal A palette of colours for violin plot.
#' 
#' @returns A patchwork object of 2 violin plots.
#' 
#' @examples
#' plot_paired_vlns(seurat_small, 'harmony_clusters', c('GAD1', 'GAD2'), c('DRD1', 'DRD2'))
#' 
plot_paired_vlns <- function(
    
  seurat_obj = NULL,
  meta_id = NULL,
  genes_a = NULL,
  genes_b =NULL,
  col_pal = c("#FF7373", "#CC4E3D", "#993617", "#FF5C00", "#CC925C", "#996E2E", 
              "#FFC226", "#CCAB00", "#999645", "#F1FF4D", "#A9CC1F", "#689900", 
              "#C1FF73", "#7CCC3D", "#419917", "#33FF00", "#65CC5C", "#2E9932", 
              "#26FF49", "#00CC39", "#459967", "#4DFFA9", "#1FCC8E", "#009974", 
              "#73FFEE", "#3DCCCC", "#178999", "#00C2FF", "#5CA4CC", "#2E6699", 
              "#267DFF", "#0039CC", "#455299", "#4D54FF", "#2C1FCC", "#1F0099", 
              "#A073FF", "#7C3DCC", "#601799", "#AD00FF", "#B65CCC", "#902E99", 
              "#FF26F6", "#CC00AB", "#994581", "#FF4DB8", "#CC1F72", "#990037", 
              "#FF7394", "#CC3D4E")
  
) {
  
  message('Creating Vlns to compare markers:\n')
  vln_a <- create_stacked_vln_plot(seurat_obj, meta_id, genes_a,
                                   toupper(region), col_pal)
  vln_b <- create_stacked_vln_plot(seurat_obj, meta_id, genes_b,
                                   toupper(region), col_pal)
  
  vln_a | vln_b
  
} 

# Add functionality to add a prefix to colnames
# Add functionality to extract normalised values
calculate_average_expression <- function(
    
  seurat_obj = NULL,
  gene_list = NULL
  
) {
  
  message('Calculating average expression ...')
  av_exp_mat <- AverageExpression(seurat_obj, layer = 'counts', features = gene_list)
  av_exp_mat <- av_exp_mat$RNA
  #colnames(av_exp_mat) <- paste0(prefix, '_', colnames(av_exp_mat))
  
  return(av_exp_mat)
  
}

# Add functionality to add a prefix to colnames
# Add functionality to extract normalised values
calculate_aggregated_expression <- function(
    
  seurat_obj = NULL,
  gene_list = NULL
  
) {
  
  message('Calculating aggregated feature expression ...')
  agg_exp_mat <- AggregateExpression(seurat_obj, features = gene_list)
  agg_exp_mat <- agg_exp_mat$RNA
  #colnames(agg_exp_mat) <- paste0(prefix, '_', colnames(agg_exp_mat))
  
  return(agg_exp_mat)
  
}

#' Project Seurat sketch object data onto entire Seurat object 
#' 
#' @param seurat_obj An uncorrected Seurat object.
#' @param dimensions Number of principal components.
#' @param reduction A reduction object in the sketch object to use for projection.
#' @param umap_model The umap model in the sketch object to use for the projection.
#' @param cluster_model The cluster model in the sketch object to use for the projection.
#' 
#' @returns A Seurat object.
#' 
#' @examples
#' project_sketch_data(seurat_small, 30, 'harmony', 'umap.harmony', harmony_clusters_0.1')
#' 
project_sketch_data <- function(
    
  seurat_obj = NULL,
  dimensions = 30,
  reduction = NULL,
  umap_model = NULL,
  cluster_model = NULL
  
) {
  
  message('Projecting integration data from sketch to all cells: ', region, ' ...')
  seurat_obj <- ProjectIntegration(object =  seurat_obj, 
                                   sketched.assay = "sketch", 
                                   assay = "RNA", 
                                   reduction = reduction)
  
  message('Projecting cell labels from sketch to all cells: ', region, ' ...')
  seurat_obj <- ProjectData(
    object = seurat_obj,
    sketched.assay = "sketch",
    assay = "RNA",
    sketched.reduction = paste0(reduction, '.full'),
    full.reduction = paste0(reduction, '.full'),
    dims = 1:dimensions,
    umap.model = umap_model,
    refdata = list(cluster_full = cluster_model)
  )
  
  message('Running UMAP for full object:')
  seurat_obj <- RunUMAP(seurat_obj, reduction = paste0(reduction, '.full'), 
                        dims = 1:dimensions, reduction.name = "umap.full",
                        reduction.key = "UMAPfull_")
  
  message('Seurat object after project integration:')
  message(paste0(capture.output(seurat_obj), collapse = "\n"), '\n')
  
  return(seurat_obj)
  
}


# Function to create logfiles in snakemake
log_smk <- function() {
  if (exists("snakemake") & length(snakemake@log) != 0) {
    log <- file(snakemake@log[1][[1]], open = "wt")
    sink(log, append = TRUE)
    sink(log, append = TRUE, type = "message")
  }
}

## Functions for hdWGCNA

# Catch cell types removed by hdWGCNA
catch_clusters_rm_warning <- function() {
  
  cell_split <- str_split(names(last.warning), c(', ', ': '))
  cell_split <- cell_split[[2]][2]
  cell_split <- unlist(str_split(cell_split, c(', ')))
    
  return(cell_split )
  
}

#' Compare cell type proportions between Seurat cell types and hdWGCNA metacell types
#' 
#' @param seurat_obj An Seurat object.
#' @param meta_obj A Seurat object generated by hdWGCNA::GetMetacellObject(seurat_obj)
#' 
#' @returns A list. 2 dfs. 1st, cell type proportions in original Seurat object. 2nd,
#' meta-cell type proportions generate by hdWGCNA.
#' 
#' @examples
#' get_wgcna_cell_props(seurat_small, metacell_small)
#' 
get_wgcna_cell_props <- function(
    
  seurat_obj = NULL,
  meta_obj = NULL
  
) {
  
  # Compare cell- and meta cell proportions
  orig_prop_tbl <- seurat_obj$cellIDs %>%
    as_tibble() %>%
    group_by(value) %>% 
    summarise(Count = n()) %>%
    mutate(Percent = Count / sum(Count) * 100)
  
  meta_prop_tbl <- colnames(meta_obj) %>%
    as_tibble() %>%
    separate(value, c('value', NA), sep = '#') %>%
    group_by(value) %>% 
    summarise(Count = n()) %>%
    mutate(Percent = Count / sum(Count) * 100) %>%
    arrange(desc(Percent))
  
  prop_list <- list(
    'orig' = orig_prop_tbl,
    'meta' = meta_prop_tbl
  )
  
  return(prop_list)  
  
}

# Run WGCNA original 
# This follows the hdWGCNA pipeline in the vinegette 
# Saves a Seurat object to file for each cell type
run_wgcna_orig <- function(
    
  seurat_obj = NULL,
  cell_types = NULL,
  region = NULL,
  outdir = NULL
  
) {
  
  # Get cell
  
  for (cell_type in cell_types) {
    
    sink(paste0(outdir, region, '_', cell_type, '_hdWGCNA.log'))
    message("\n\n\nRunning hdWGCNA for: ", toupper(region), ', ', cell_type)
    
    # Required to prevent running completed runs after error
    if (file.exists(paste0(outdir, cell_type, '_TOM.rda'))) {
      message("Skipping ", cell_type, " as TOM exists already")
      sink()
      next
    }
    
    # Set up the expression matrix
    seurat_object <- SetDatExpr(
      seurat_obj,
      group_name = cell_type, 
      group.by = 'cellIDs', 
      assay = 'RNA', # using RNA assay
      slot = 'data' # using normalized data
    )
    
    # Select soft power threshold
    seurat_object <- TestSoftPowers(
      seurat_object,
      networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
    )
    
    power_val <- GetPowerTable(seurat_object) %>%
      select(Power, SFT.R.sq) %>%
      filter(SFT.R.sq > 0.8) %>%
      pull(Power) %>%
      dplyr::first()
    
    message("Soft Power threshold set to: ", power_val)
    
    if (power_val < 1 | power_val > 50 | is.na(power_val)) {
      message('Skipping ', cell_type, 'due to spurious Power value')
      print(cell_type, ': Stopped due to spurious power value')
      sink()
      next 
      }
    
    # Construct co-expression network
    seurat_object <- ConstructNetwork(
      seurat_object,
      tom_name = cell_type, # name of the topoligical overlap matrix written to disk
      soft_power = power_val,
      overwrite_tom = T,
      tom_outdir = outdir
    )
    
    # Compute Eigengenes and Connectivity
    # Compute all MEs in the full single-cell dataset
    seurat_object <- ModuleEigengenes(
      seurat_object,
      group.by.vars = "Sample"
    )
    
    # Compute module connectivity
    # compute eigengene-based connectivity (kME):
    seurat_object <- ModuleConnectivity(
      seurat_object,
      group.by = 'cellIDs', 
      group_name = cell_type
    )
    
    # rename the modules
    seurat_object <- ResetModuleNames(
      seurat_object,
      new_name = paste0(cell_type, "-M")
    )
    
    saveRDS(seurat_object, file = paste0(outdir, toupper(region), 
                                         '_', cell_type, '_hdWGCNA.rds'))
    rm(seurat_object)
    sink()
    
  }
  
}

# Run WGCNA
# This is my attempt so far to save all WGCNA data to the same 
# Seurat object 
run_wgcna <- function(
    
  seurat_obj = NULL,
  cell_types = NULL,
  region = NULL,
  outdir = NULL
  
) {
  
  for (i in 1:length(cell_types)) {
    
    sink(paste0(outdir, region, '_', cell_types[i], '_hdWGCNA.log'))
    message("\n\n\nRunning hdWGCNA for: ", toupper(region), ', ', cell_types[i])
    
    # Set up the expression matrix 1st seurat_obj the seurat_object
    if (i == 1)
      seurat_object <- SetDatExpr(
        seurat_obj,
        group_name = cell_types[i], 
        group.by = 'cellIDs', 
        assay = 'RNA', # using RNA assay
        slot = 'data', # using normalized data
        wgcna_name = paste0(cell_types[i], '_wgcna')
      )
    else
      seurat_object <- SetDatExpr(
        seurat_object,
        group_name = cell_types[i], 
        group.by = 'cellIDs', 
        assay = 'RNA', # using RNA assay
        slot = 'data', # using normalized data
        wgcna_name = paste0(cell_types[i], '_wgcna')
      )
      
    
    # Select soft power threshold
    seurat_object <- TestSoftPowers(
      seurat_object,
      networkType = 'signed', # you can also use "unsigned" or "signed hybrid"
      wgcna_name = paste0(cell_types[i], '_wgcna'))
    
    power_val <- GetPowerTable(seurat_object, paste0(cell_types[i], '_wgcna')) %>%
      select(Power, SFT.R.sq) %>%
      filter(SFT.R.sq > 0.8) %>%
      pull(Power) %>%
      dplyr::first()
    
    message("Soft Power threshold set to: ", power_val)
    
    message("Construct co-expression network ...")
    # Construct co-expression network
    seurat_object <- ConstructNetwork(
      seurat_object,
      tom_name = cell_types[i], # name of the topoligical overlap matrix written to disk
      soft_power = power_val,
      overwrite_tom = T,
      wgcna_name = paste0(cell_types[i], '_wgcna'),
      tom_outdir = outdir
    )
    message("Compute Eigengenes ...")
    # Compute Eigengenes and Connectivity
    # Compute all MEs in the full single-cell dataset
    seurat_object <- ModuleEigengenes(
      seurat_object,
      modules = GetModules(seurat_object, wgcna_name = paste0(cell_types[i], '_wgcna')),
      group.by.vars = "Sample",
      wgcna_name = paste0(cell_types[i], '_wgcna')
    )
    
    message("Module Connect ...")
    # Compute module connectivity
    # compute eigengene-based connectivity (kME):
    seurat_object <- ModuleConnectivity(
      seurat_object,
      group.by = 'cellIDs', 
      group_name = cell_types[i],
      wgcna_name = paste0(cell_types[i], '_wgcna')
    )
    the 
    message("Rename modules ...")
    # rename the modules
    seurat_object <- ResetModuleNames(
      seurat_object,
      new_name = paste0(cell_types[i], "-M"),
      wgcna_name = paste0(cell_types[i], '_wgcna')
    )
    
  }
  
  message("Saving RDS file ...")
  saveRDS(seurat_object, file = paste0(outdir, toupper(region), 
                                       '_', cell_types[i], '_hdWGCNA.rds'))
  rm(seurat_object)
  sink()
  
}

#' Set up object for hdWGCNA and construct metacells
#' Note this has to be run twice. Once to set up the general WGCNA metadata slot
#' in the Seurat object. Then once the meta cells have been established, run again
#' to add a seprate metadata slot for each meta cell. 
#' 
#' @param seurat_object An Seurat object.
#' @param gene_select A string. Either 'fraction', or 'variable' specifying method
#' by which genes are selected to create GeX networks.
#' @param wgcna_name A string. Name for the WGCNA object stored in seurat_object@misc. 
#' @param metacell_location A string. Name of the WGCNA experiment to copy the metacell object from. 
#' Set this to null when first setting up object.
#' @param meta_cell_types A vector of the cell types that survived min_cells cull 
#' after first running this function. Set this to null when first setting up object.
#' @param set_k A string. No of cell types to merge for kNN.
#' 
#' @returns A Seurat object.
#' 
#' @examples
#' create_wgcna_metacells(seurat_small, 'fraction', 'wgcna')
#' create_wgcna_metacells(seurat_small, 'fraction', 'wgcna')
#' 
create_wgcna_metacells <- function(
    
  seurat_object = NULL,
  gene_select = 'fraction',
  wgcna_name = 'wgcna',
  metacell_location = NULL,
  meta_cell_types = NULL,
  set_k = 25
  
) {
  

  
  if (is.null(metacell_location)) {
    
    message('Setting up Seurat object for WGCNA ...')
    seurat_obj <- hdWGCNA::SetupForWGCNA(
        seurat_object,
        gene_select = gene_select, # Default fraction is 0.05
        wgcna_name = wgcna_name,
        metacell_location = metacell_location)
    
    message('Constructing metacells for WGCNA ...')
    seurat_obj <- hdWGCNA::MetacellsByGroups(
      seurat_obj = seurat_obj,
      group.by = c("cellIDs", "Sample"), # specify the columns in seurat_obj@meta.data to group by
      reduction = 'umap', # select the dimensionality reduction to perform KNN on
      k = set_k, # nearest-neighbors parameter
      max_shared = 10, # maximum number of shared cells between two metacells
      ident.group = 'cellIDs', # set the Idents of the metacell seurat object
    )
    
    # Normalize metacell expression matrix (only required on general misc object)
    seurat_obj <- NormalizeMetacells(seurat_obj)
  
    return(seurat_obj)
    
  } else {
    
    for (i in 1:length(meta_cell_types)) {
      
      if (i == 1) {

        message('Setting up Seurat object for ', meta_cell_types[i], ' ...')
        seurat_all <- hdWGCNA::SetupForWGCNA(
          seurat_object,
          gene_select = gene_select, # Default fraction is 0.05
          wgcna_name = paste0(meta_cell_types[i], '_wgcna'),
          metacell_location = metacell_location)
      
      } else {

        message('Setting up Seurat object for ', meta_cell_types[i], ' ...')
        seurat_all <- hdWGCNA::SetupForWGCNA(
          seurat_all,
          gene_select = gene_select, # Default fraction is 0.05
          wgcna_name = paste0(meta_cell_types[i], '_wgcna'),
          metacell_location = metacell_location) }}
    
  return(seurat_all)}
  
}

#' Test for overlap of dystonia genes and wgcna coexpression modules
#' 
#' @param seurat_object An Seurat object.
#' @param meta_cells A string. The cellIDs retained after generating meta cells.
#' @param region A string spcifying brain region of interest.
#' @param wgcna_name A string. Name for the WGCNA object stored in seurat_object@misc. 
#' @param obj_dir A string where the suerat object is stored. 
#' 
#' @returns A list of overlapping genes for each metacell type in the surat object.
#' 
#' @examples
#' run_dyst_gene_overlap(seurat_obj, meta_cell_types, region, paste0(region, '_wgcna'), wgcna_dir)
#' 
run_dyst_gene_overlap <- function(
    
  seurat_object = NULL,  
  meta_cells = NULL,
  region = NULL,
  wgcna_name = NULL,
  obj_dir = NULL
  
) {
  
  overlap_list <- list()
  
  for (cell_type in meta_cells) {
    
    seurat_object <- readRDS(paste0(obj_dir, toupper(region), 
                                    '_', cell_type, '_hdWGCNA.rds'))
    
    # Get the module assignment table
    module_names <- GetModules(seurat_object, wgcna_name) %>% 
      subset(module != 'grey') %>%
      select(-gene_name, -module, -color, -kME_grey) %>%
      rename_with(~gsub("kME_", "", .x, fixed = TRUE)) %>%
      colnames()
    
    for (mod in module_names) {
      
      print(mod)
      genes <- GetHubGenes(seurat_object, n_hubs = 50) %>%
        filter(module == mod) %>%
        filter(gene_name %in% dystonia_genes) 
      print(genes)
      
      overlap_list[[mod]] <- genes
      
    }
    
  }
  
  return(overlap_list)
  
}


#' Get wgcna stats from seurat_object@misc
#' 
#' @param seurat_object An Seurat object.
#' @param meta_cells A string. The cellIDs retained after generating meta cells.
#' @param region A string spcifying brain region of interest.
#' @param wgcna_name A string. Name for the WGCNA object stored in seurat_object@misc. 
#' @param obj_dir A string where the suerat object is stored. 
#' 
#' @returns A tibble.
#' 
#' @examples
#' get_wgcna_stats(seurat_obj, meta_cell_types, region, paste0(region, '_wgcna'), wgcna_dir)
#' 
get_wgcna_stats <- function(
    
  seurat_object = NULL,  
  meta_cells = NULL,
  region = NULL,
  wgcna_name = NULL,
  obj_dir = NULL
  
) {
  
  for (cell_type in meta_cells) {
    
    message("Getting stats for ", toupper(region), ', ', cell_type, '...')
    seurat_object <- readRDS(paste0(obj_dir, toupper(region), 
                                    '_', cell_type, '_hdWGCNA.rds'))
    
    if (!exists('stats_tbl')) {
      
      stats_tbl <- tibble(
        'Cell_type' = cell_type,
        'Metacells' = ncol(seurat_object@misc[[wgcna_name]][["wgcna_metacell_obj"]]),
        'Power' = seurat_object@misc[[wgcna_name]][["wgcna_params"]][["power"]],
        'Modules' = ncol(seurat_object@misc[[wgcna_name]][["wgcna_modules"]]) - 4
      )
      
    } else {
      
      stats_tbl <- rbind(stats_tbl,
                         tibble(
                           'Cell_type' = cell_type,
                           'Metacells' = ncol(seurat_object@misc[[wgcna_name]][["wgcna_metacell_obj"]]),
                           'Power' = seurat_object@misc[[wgcna_name]][["wgcna_params"]][["power"]],
                           'Modules' = ncol(seurat_object@misc[[wgcna_name]][["wgcna_modules"]]) - 4
                         ))
      
    }
    
  }
  
  return(stats_tbl)
  
}

# Aggregate cell types for wgcna
recode_wgcna_clusters <- function(
    
  seurat_object = NULL,
  region = NULL
  
) {
  
  message('Aggregating cells for ', region, '...')
  
  if (region == 'fcx_fetal') {
    
    seurat_obj$cellIDs <- seurat_obj@meta.data %>% 
      tibble::as_tibble() %>%
      dplyr::mutate(cluster_recode = recode(.data[['cellIDs']], 
                                            "FC-ExN-1" = "FC-ExN", 
                                            "FC-ExN-2" = "FC-ExN",
                                            "FC-ExN-3" = "FC-ExN",
                                            "FC-ExN-4" = "FC-ExN",
                                            "FC-ExN-5" = "FC-ExN",
                                            "FC-InN-1" = "FC-InN",
                                            "FC-InN-2" = "FC-InN",
                                            "FC-InN-3" = "FC-InN",
                                            "FC-InN-4" = "FC-InN",
                                            "FC-RG-1" = "FC-RG",
                                            "FC-RG-2" = "FC-RG")) %>%
      pull(cluster_recode)
    
  }
  
  if (region == 'ge_fetal') {
    
    seurat_obj$cellIDs <- seurat_obj@meta.data %>% 
      tibble::as_tibble() %>%
      dplyr::mutate(cluster_recode = recode(.data[['cellIDs']], 
                                            "FC-InN-1" = "FC-InN",
                                            "FC-InN-2" = "FC-InN",
                                            "FC-InN-3" = "FC-InN",
                                            "FC-InN-4" = "FC-InN",
                                            "FC-InN-5" = "FC-InN",
                                            "FC-InN-6" = "FC-InN",
                                            "FC-InN-7" = "FC-InN",
                                            "FC-RG-1" = "FC-RG",
                                            "FC-RG-2" = "FC-RG",
                                            "FC-RG-3" = "FC-RG")) %>%
      pull(cluster_recode)
    
  }
  
  if (region == 'cer_fetal') {
    
    seurat_obj$cellIDs <- seurat_obj@meta.data %>% 
      tibble::as_tibble() %>%
      dplyr::mutate(cluster_recode = recode(.data[['cellIDs']], 
                                            "Cer-ExN-1" = "Cer-ExN", 
                                            "Cer-ExN-2" = "Cer-ExN",
                                            "Cer-ExN-3" = "Cer-ExN",
                                            "Cer-ExN-4" = "Cer-ExN",
                                            "Cer-ExN-5" = "Cer-ExN",
                                            "Cer-InN-1" = "Cer-InN",
                                            "Cer-InN-2" = "Cer-InN",
                                            "Cer-InN-3" = "Cer-InN",
                                            "Cer-InN-4" = "Cer-InN",
                                            "Cer-InN-5" = "Cer-InN",
                                            "Cer-InN-6" = "Cer-InN",
                                            "Cer-InN-7" = "Cer-InN",
                                            "Cer-InN-8" = "Cer-InN",
                                            "Cer-RG-1" = "Cer-RG",
                                            "Cer-RG-2" = "Cer-RG",
                                            "Cer-RG-3" = "Cer-RG")) %>%
      pull(cluster_recode)
    
  }
  
  message("New cell IDs are: ", seurat_obj$cellIDs)
  return(seurat_object)
  
} 