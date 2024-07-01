#--------------------------------------------------------------------------------------
#
#    Dystonia - hdWGCNA
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# Vingette: https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html

##  Load Packages, functions and variables  -------------------------------------------
message('Setting environment variables ...')
if (Sys.info()[["nodename"]] == "Darrens-iMac-2.local") {
  
  library(yaml)
  root_dir <- '~/Desktop/dystonia_snRNAseq_2024/'
  yaml_file <- yaml.load_file(paste0(root_dir, 'config/config.yaml'))
  region <- yaml.load(yaml_file$region)
  
  source(paste0(root_dir, 'workflow/scripts/dystonia_functions.R'))
  source(paste0(root_dir, 'workflow/scripts/dystonia_Renvs.R'))
  source(paste0(root_dir, 'workflow/scripts/dystonia_gene_lists.R'))
  
} else {
  
  source('scripts/dystonia_functions.R')
  source('scripts/dystonia_Renvs.R')
  source('scripts/dystonia_gene_lists.R')
  
}

### Need to configure these for hawk ####

# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 4)

# load the Zhou et al snRNA-seq dataset
#seurat_obj <- readRDS('Zhou_2020.rds')

DATA_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/resources/R_objects/"
FIG_DIR <- "~/Desktop/fetal_brain_snRNAseq_110122/results/figures/"
REGIONS <- c('hip', 'pfc', 'wge')

### Need to configure these for hawk ####

## Load Data --------------------------------------------------------------------------
# Seurat objects  ----
for (REGION in REGIONS) { 
  
  seurat.obj <- readRDS(paste0(DATA_DIR, 'seurat.', REGION, '.final.rds'))
  assign(paste0('seurat.', REGION), seurat.obj, .GlobalEnv)
  
}

# Set up object for WGCNA
seurat_obj <- SetupForWGCNA(
  seurat.pfc,
  gene_select = "fraction", # genes exp in frac of cells of Seurat var genes
  fraction = 0.05, 
  wgcna_name = "pfc_wcgna",
  min_cells = 100
)

# Construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("cellIDs", "Sample"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cellIDs', # set the Idents of the metacell seurat object
)



# Catch cell types that are removed due to low cell counts
cell_split <- str_split(names(last.warning), c(', ', ': '))
cell_split <- cell_split[[2]][2]
cell_split <- unlist(str_split(cell_split, c(', ')))

# Compare cell- and meta cell proportions
seurat_obj$cellIDs %>%
  as_tibble() %>%
  group_by(value) %>% 
  summarise(Count = n()) %>%
  mutate(Percent = Count / sum(Count) * 100)

meta_obj <- GetMetacellObject(seurat_obj)
meta_cells <- colnames(meta_obj) %>%
  as_tibble() %>%
  separate(value, c('value', NA), sep = '#') %>%
  group_by(value) %>% 
  summarise(Count = n()) %>%
  mutate(Percent = Count / sum(Count) * 100) %>%
  arrange(desc(Percent))

# Normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

for (cell_type in meta_cells %>% pull(value)) {
  
  message("Running hdWGCNA for: ", cell_type)
  
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
  
  # Construct co-expression network
  seurat_object <- ConstructNetwork(
    seurat_object,
    tom_name = cell_type, # name of the topoligical overlap matrix written to disk
    soft_power = power_val,
    overwrite_tom = T
  )
  
  # Compute Eigengenes and Connectivity
  # Compute all MEs in the full single-cell dataset
  seurat_object <- ModuleEigengenes(
    seurat_object,
    group.by.vars = "Sample"
  )
  
  # harmonized module eigengenes:
  hMEs <- GetMEs(seurat_object)
  # module eigengenes:
  MEs <- GetMEs(seurat_object, harmonized = FALSE)
  
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
  
  saveRDS(seurat_obj, file = paste0('~/Desktop/hdWCGNA/', cell_type, '_hdWGCNA.rds'))
  rm(seurat_object)
  
  }
  
}


PlotDendrogram(seurat_object, main = paste0(cell_type, ' hdWGCNA Dendrogram'))
TOM <- GetTOM(seurat_object)
TOM[1:10, 1:10]
dim(TOM)

# plot genes ranked by kME for each module
# kME is the correlation of a gene's expression with a module eigengene.
p <- PlotKMEs(seurat_object, ncol = 5)
p

for (cell_type in meta_cells %>% pull(value)) {
  
  seurat_object <- readRDS(file = paste0('~/Desktop/hdWCGNA/', 
                                  cell_type, '_hdWGCNA.rds'))

  # Get the module assignment table
  module_names <- GetModules(seurat_object) %>% 
    subset(module != 'grey') %>%
    select(-gene_name, -module, -color, -kME_grey) %>%
    rename_with(~gsub("kME_", "", .x, fixed = TRUE)) %>%
    colnames()
  
  test <- GetHubGenes(seurat_object, n_hubs = 10)
  boxplot(test$kME)
  
  for (mod in module_names) {
    
    print(mod)
    genes <- GetHubGenes(seurat_object, n_hubs = 50) %>%
      filter(module == mod) %>%
      filter(gene_name %in% dystonia_genes) 
    print(genes)

  }
  
}
  
# saveRDS(seurat_obj, file='hdWGCNA_object.rds')

## ------------------------------------------------------------------------------------
## ------------------------------------------------------------------------------------

# Compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# Visualization
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order = TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

# Hub genes
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order = TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

seurat_obj$cluster <- do.call(rbind, strsplit(as.character(seurat_obj$cellIDs), ' '))[,1]

ModuleRadarPlot(
  seurat_obj,
  group.by = 'cellIDs',
  barcodes = seurat_obj@meta.data %>% subset(cell_type == 'FC-ExN-1') %>% rownames(),
  axis.label.size = 4,
  grid.label.size = 4
)


# Testing joining vln plots
FC_plot <- VlnPlot(seurat.pfc, dystonia_genes, stack = TRUE, flip = TRUE, 
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
  ggtitle("Frontal Cortex")

GE_plot <- VlnPlot(seurat.wge, dystonia_genes, stack = TRUE, flip = TRUE, 
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
  ggtitle("Ganglionic Eminence")

Hip_plot <- VlnPlot(seurat.hip, dystonia_genes, stack = TRUE, flip = TRUE, 
                    same.y.lims = TRUE, fill.by = 'ident')  +
  theme(legend.position = "none",
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 18),
        axis.text.x  = element_text(colour = "#000000", size = 16),
        axis.text.y  = element_text(colour = "#000000", size = 16)) +
  xlab('Cell type') +
  ggtitle("Hippocampus")

cowplot::plot_grid(FC_plot, GE_plot, Hip_plot, ncol = 3)
# -------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------
