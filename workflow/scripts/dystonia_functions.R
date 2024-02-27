# Load data from Stiletti data repository
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
  stopifnot("anns_table must be a Stiletti downloads tibble" = is_tibble(anns_table))
  
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
  root_directory = NULL
  
){
  
  file_set <- annotations
  data_list <- c()
  metadata_list <- c()
  
  for (i in 1:length(file_set)) {
    path <- paste0(root_directory, file_set[i], file_type)
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
  
  seurat_obj@meta.data %>%
    as_tibble(rownames = 'cell_id') %>%
    dplyr::select(any_of(c(meta_col))) %>%
    group_by(.data[[meta_col]]) %>%
    count()
  
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
  sum_outlier <- scuttle::isOutlier(sce_obj$sum,  nmads = mad_thresh, 
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
  
  sce_obj$cell_outlier <- cell_outliers
  
  message('Cell numbers that will be excluded at specified thresholds:')
  message(paste0(capture.output(outlier_cnts_tbl), collapse = "\n"), '\n')
  
  # Plot outliers
  message('Plotting ...')
  create_outlier_plots(sce_obj)
  
  return(sce_obj)
  
  
}

# Create outlier plots, run within `get_cell_outliers` function
# Need to fix the top 50 gene plot
# Should add choice of meta data column for x axis
create_outlier_plots <- function(
    
  sce_obj = NULL
  
) {
  
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


create_sketch_object <- function(
    
  seurat_obj = NULL, 
  dims = 30
  
) {
  
  seurat_obj <- Seurat::NormalizeData(seurat_obj) %>%
    FindVariableFeatures(verbose = FALSE)
  
  seurat_sketch <- Seurat::SketchData(object = seurat_obj, 
                                      ncells = 5000, 
                                      method = "LeverageScore", 
                                      sketched.assay = "sketch")
  
  # Take a small subset of data to 
  Seurat::DefaultAssay(seurat_sketch) <- "sketch"
  
  # Basic commands
  seurat_sketch <- Seurat::FindVariableFeatures(seurat_sketch, verbose = F) %>%
    Seurat::ScaleData(verbose = F) %>%
    Seurat::RunPCA(verbose = F) %>%
    Seurat::FindNeighbors(dims = 1:dims) %>%
    Seurat::FindClusters() %>%
    Seurat::RunUMAP(dims = 1:dims)
  
  return(seurat_sketch)
  
}

# Basic QC plot
create_basic_qc_plots <- function(seurat_obj = NULL,
                                  point_size = 0) {
  
  
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
  dims = 30
  
) {
  
  cluster_plot <- Seurat::DimPlot(seurat_obj, reduction = "umap")
  elbow_plot <- Seurat::ElbowPlot(seurat_obj, ndims = dims)
  dataset_plot <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = 'dataset') 
  donor_plot <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = 'sample_id')
  bar_plot_dataset <- create_proportion_barplot(seurat_obj, paste0('seurat_clusters'), meta_id = 'dataset')
  bar_plot_donor <- create_proportion_barplot(seurat_obj, paste0('seurat_clusters'), meta_id = 'sample_id')
  
  qc_plot <- cowplot::plot_grid(cluster_plot, elbow_plot, 
                                dataset_plot, bar_plot_dataset,
                                donor_plot, bar_plot_donor, ncol = 2)
  
}

run_integration <- function(
    
  seurat_obj = NULL,
  reductions = 'harmony',
  dimensions = 30
  
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
    FindClusters(cluster.name = "harmony_clusters") %>%
    RunUMAP(reduction = "harmony", dims = 1:dimensions, reduction.name = "umap.harmony")
  
  }
  
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
    FindClusters(cluster.name = "cca_clusters") %>%
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
    FindClusters(cluster.name = "rpca_clusters") %>%
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
    FindClusters(cluster.name = "fastmnn_clusters") %>%
    RunUMAP(reduction = "fastmnn", dims = 1:dimensions, reduction.name = "umap.fastmnn")
  
  }
  
  return(seurat_obj)
  
}

create_proportion_barplot <- function(seurat_obj = NULL, 
                                      cluster_id = NULL,
                                      meta_id = NULL) {
  
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

create_integration_plot <- function(
    
  seurat_obj = NULL, 
  reductions = c('harmony'), 
  meta_id = 'sample_id',
  dims = 30
  
) {
  
  plot_list <- list()
  
  for (i in 1:length(reductions)) {
    
    cluster_plot <- DimPlot(seurat_obj, reduction = paste0("umap.", reductions[i]), group.by = paste0(reductions[i], '_clusters'))
    meta_plot <- DimPlot(seurat_obj, reduction = paste0("umap.", reductions[i]), group.by = meta_id)
    bar_plot <- create_proportion_barplot(seurat_obj, paste0(reductions[i], '_clusters'), meta_id)
    
    plot_list[[paste0(reductions[i], '_cluster')]]  <- cluster_plot
    plot_list[[paste0(reductions[i], '_meta')]]  <- meta_plot
    plot_list[[paste0(reductions[i], '_bar')]]  <- bar_plot
    
  } 
  
  group_plot <- plot_grid(plotlist = plot_list, ncol = 2)
  
  return(group_plot)
  
}

create_stacked_vln_plot <- function(
    
  seurat_obj = NULL,
  set_ident = 'seurat_clusters',
  genes = NULL,
  plot_title = NULL,
  col_pal = NULL
  
) {
  
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
    ggtitle(plot_title)
  
}

calculate_average_expression <- function(
    
  seurat_obj = NULL,
  region = NULL,
  gene_list = NULL
  
) {
  
  message('Calculating average expression ...')
  av_exp_mat <- AverageExpression(seurat_obj, layer = 'counts', features = gene_list)
  av_exp_mat <- av_exp_mat$RNA
  colnames(av_exp_mat) <- paste0(region, '_', seq(0, ncol(av_exp_mat) - 1, 1))
  
  return(av_exp_mat)
  
}

# Project info generated on Seurat sketch object onto entire dataset
project_sketch_data <- function(
    
  seurat_obj = NULL,
  dimensions = 30
  
) {
  
  seurat_obj <- ProjectData(
    object = seurat_obj,
    assay = "RNA",
    full.reduction = "pca.full",
    sketched.assay = "sketch",
    sketched.reduction = "pca",
    umap.model = "umap",
    dims = dimensions,
    refdata = list(cluster_full = "seurat_clusters")
    
  )
  
  seurat_obj <- ProjectIntegration(object = seurat_obj,
                                   sketched.assay = "sketch",
                                   assay = "RNA",
                                   reduction = "harmony")
  
  return(seurat_obj)
  
}