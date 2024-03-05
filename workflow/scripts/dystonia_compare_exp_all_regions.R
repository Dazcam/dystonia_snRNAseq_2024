#--------------------------------------------------------------------------------------
#
#    Dystonia - Compare average and aggregate expression
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Script to compare expression of dystonia genes across all regions / ages
#  Compare average and aggregate expression of 25 dystonia genes across all regions

##  Load Packages  --------------------------------------------------------------------
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(SeuratWrappers)
library(Azimuth) 
library(tidyverse)
library(scCustomize)
library(readxl)
library(cowplot)

## Set variables  ---------------------------------------------------------------------
fetal_dir <- '~/Desktop/fetal_brain_snRNAseq_110122/resources/R_objects/'
R_dir <- '~/Desktop/dystonia_snRNAseq_2024/results/01R/'

## Load Data --------------------------------------------------------------------------
# Fetal data
for (region in c('cer', 'hip', 'pfc', 'wge', 'tha')) {
  
  seurat_obj <- readRDS(paste0(fetal_dir, 'seurat.', region, '.final.rds'))
  aver_exp_mat <- calculate_average_expression(seurat_obj, paste0(region, '_fetal'), dystonia_genes)
  aggr_exp_mat <- calculate_aggregated_expression(seurat_obj, paste0(region, '_fetal'), dystonia_genes)
  assign(paste0('fetal_aver_', region), aver_exp_mat, envir = .GlobalEnv)
  assign(paste0('fetal_aggr_', region), aggr_exp_mat, envir = .GlobalEnv)
  
}

# Adult data
for (region in c('fcx', 'str')) {
  
  seurat_obj_av <- readRDS(paste0(R_dir, "seurat_aver_exp_", region, ".Rds"))
  seurat_obj_ag <- readRDS(paste0(R_dir, "seurat_aggr_exp_", region, ".Rds"))

  assign(paste0('adult_aver_', region), seurat_obj_av, envir = .GlobalEnv)
  assign(paste0('adult_aggr_', region), seurat_obj_ag, envir = .GlobalEnv)
  
}

##  Join data 
av_exp_mat <- cbind(fetal_aver_cer, fetal_aver_pfc, fetal_aver_wge, fetal_aver_hip, fetal_aver_tha,
                    adult_aver_fcx, adult_aver_str)
ag_exp_mat <- cbind(fetal_aggr_cer, fetal_aggr_pfc, fetal_aggr_wge, fetal_aggr_hip, fetal_aggr_tha,
                    adult_aggr_fcx, adult_aggr_str)

# Scale data
av_exp_mat <- normalizeCounts(av_exp_mat)
av_exp_mat_scaled <- t(apply(av_exp_mat, 1, scale)) # centre and scale each col (Z-Score) then t
colnames(av_exp_mat_scaled) <- colnames(av_exp_mat) 

# Normalise??
ag_exp_mat <- normalizeCounts(ag_exp_mat)
ag_exp_mat_scaled <- t(apply(ag_exp_mat, 1, scale)) # centre and scale each col (Z-Score) then t
colnames(ag_exp_mat_scaled) <- colnames(ag_exp_mat) 

# Generate Heatmaps - huge in size
jpeg(paste0('~/Desktop/dystonia_av_exp.jpg'), width = 960, height = 2500, 
     units = "px", pointsize = 12, quality = 150)
ComplexHeatmap::Heatmap(t(av_exp_mat_scaled), 
                       column_labels = rownames(av_exp_mat_scaled), 
                       cluster_columns = F, 
                       cluster_rows = T,
                       width = ncol(t(av_exp_mat_scaled))*unit(5, "mm"), 
                       height = nrow(t(av_exp_mat_scaled))*unit(5, "mm"))

dev.off()

jpeg(paste0('~/Desktop/dystonia_ag_exp.jpg'), width = 960, height = 2500, 
     units = "px", pointsize = 12, quality = 150)
ComplexHeatmap::Heatmap(t(ag_exp_mat_scaled), 
                        column_labels = rownames(ag_exp_mat_scaled), 
                        cluster_columns = F, 
                        cluster_rows = T,
                        width = ncol(t(av_exp_mat_scaled))*unit(5, "mm"), 
                        height = nrow(t(av_exp_mat_scaled))*unit(5, "mm"))
dev.off()


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

# PCA

av_exp_pca <- prcomp(t(av_exp_mat_scaled))

### test code 
gc_vst <- read.table("~/Downloads/counts_vst.txt", header = T, row.names = 1, sep = "\t")
vst_pca <- prcomp(t(gc_vst))

# Calculate variance explained
frac_var <- function(x) x^2 / sum(x^2)

av_exp_pca$sdev %>% 
  as_tibble() %>% 
  frac_var() %>% 
  mutate(Comp = colnames(av_exp_pca$x)) %>% 
  slice(1:9) %>% 
  ggplot(aes(x=Comp, y = value)) + 
  geom_bar(stat = "identity", fill = "#4DC5F9") +
  geom_hline(yintercept = 0.03, linetype = 2) +
  xlab("Principal Components") +
  scale_y_continuous(name = "Variance Explained", breaks = seq(0,0.8,0.1), labels = scales::percent_format(accuracy = 5L)) +
  theme_classic(base_size = 14)

# Prep PCA plot
genes.selected <- av_exp_pca$rotation[c(which.max(av_exp_pca$rotation[,"PC1"]), 
                                        which.min(av_exp_pca$rotation[,"PC1"]), 
                                        which.max(av_exp_pca$rotation[,"PC2"]), 
                                        which.min(av_exp_pca$rotation[,"PC2"])),
                                      c("PC1", "PC2")]
genes.selected <- genes.selected %>%
  as.data.frame() %>%
  rownames_to_column(var = "genes")
genes.selected

pca_all <- av_exp_pca$x %>%
  as_tibble(rownames = 'cell_type') %>%
  mutate(region = str_split_i(cell_type, pattern = '_', 1)) %>%
  ggplot(aes(x= PC1, y= PC2, color = region)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype=2) +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

ggplot(genes_tbl, aes(x=PC1, y=PC2)) +
  geom_point() +
  geom_segment(aes(xend=PC1, yend=PC2), x=0, y=0, color="Grey") +
  geom_label(aes(x=PC1, y=PC2, label=genes), size=2, vjust="outward") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),  panel.grid.minor = element_blank()) 


# Bioconductor method - https://www.bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html
pca_test <- PCAtools::pca(av_exp_mat_scaled, scale = F, 
                          metadata = data.frame(region = str_split_i(colnames(av_exp_mat_scaled), 'ult', 1), 
                                                row.names = colnames(av_exp_mat_scaled)))
scree_plot <- PCAtools::screeplot(pca_test, axisLabSize = 18, titleLabSize = 22)
pca_plot <- PCAtools::biplot(pca_test, showLoadings = FALSE, ntopLoadings = 25, colby = 'region',
                             hline = 0, vline = 0)
pca_load_plot <- PCAtools::biplot(pca_test, showLoadings = TRUE, ntopLoadings = 25, colby = 'region',
                                  hline = 0, vline = 0)
pairs_plot <- PCAtools::pairsplot(pca_test, colby = 'region')
plotloadings(pca_test,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1, PC2, PC3, PC4, PC5', 
             caption = 'Top 1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)

plot_grid(pca_all, pca_bio, align = 'hvlr')


