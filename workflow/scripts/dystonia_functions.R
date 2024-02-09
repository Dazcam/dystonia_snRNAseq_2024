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
      as_tibble() %>%
      group_by(dissection) %>%
      count()

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
  donor_plot <- Seurat::DimPlot(seurat_obj, reduction = "umap", group.by = 'orig.ident')
  bar_plot_dataset <- create_proportion_barplot(seurat_obj, paste0('seurat_clusters'), meta_id = 'dataset')
  bar_plot_donor <- create_proportion_barplot(seurat_obj, paste0('seurat_clusters'), meta_id = 'orig.ident')
  
  qc_plot <- cowplot::plot_grid(cluster_plot, elbow_plot, 
                                dataset_plot, bar_plot_dataset,
                                donor_plot, bar_plot_donor, ncol = 2)
  
}

run_integration_all <- function(seurat_obj = NULL) {
  
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
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction = "harmony") %>%
    FindClusters(cluster.name = "harmony_clusters") %>%
    RunUMAP(reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
  
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
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction = "cca") %>%
    FindClusters(cluster.name = "cca_clusters") %>%
    RunUMAP(reduction = "cca", dims = 1:30, reduction.name = "umap.cca")
  
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
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction = "rpca") %>%
    FindClusters(cluster.name = "rpca_clusters") %>%
    RunUMAP(reduction = "rpca", dims = 1:30, reduction.name = "umap.rpca")
  
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
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction = "fastmnn") %>%
    FindClusters(cluster.name = "fastmnn_clusters") %>%
    RunUMAP(reduction = "fastmnn", dims = 1:30, reduction.name = "umap.fastmnn")
  
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

create_integration_compare_plot <- function(
    
  seurat_obj = seurat_sk_str, 
  reductions = c('harmony', 'cca', 'rpca', 'fastmnn'), 
  meta_id = 'dataset',
  dims = 30
  
) {
  
  plot_list <- list()
  
  for (i in 1:length(reductions)) {
    
    cluster_plot <- scCustomize::DimPlot_scCustom(seurat_obj, reduction = paste0("umap.", reductions[i]), group.by = paste0(reductions[i], '_clusters'))
    meta_plot <- scCustomize::DimPlot_scCustom(seurat_obj, reduction = paste0("umap.", reductions[i]), group.by = meta_id)
    bar_plot <- create_proportion_barplot(seurat_obj, paste0(reductions[i], '_clusters'), meta_id)
    
    plot_list[[paste0(reductions[i], '_cluster')]]  <- cluster_plot
    plot_list[[paste0(reductions[i], '_meta')]]  <- meta_plot
    plot_list[[paste0(reductions[i], '_bar')]]  <- bar_plot
    
  } 
  
  group_plot <- plot_grid(plotlist = plot_list, ncol = 3)
  
  return(group_plot)
  
}

create_stacked_vln_plot <- function(
    
  seurat_obj = NULL,
  genes = NULL,
  plot_title = NULL
  
) {
  
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
  VlnPlot(seurat_obj, genes, stack = TRUE, flip = TRUE, 
          same.y.lims = TRUE, fill.by = 'ident') +
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