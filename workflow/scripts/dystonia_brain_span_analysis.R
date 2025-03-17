#--------------------------------------------------------------------------------------
#
#    Dystonia - Brain span linear regression analysis
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

# We have gene expression data and metadata for:
  
# 39 individuals
# 447 independent dissections (see `dystonia_brain_span_prep_data.r`)
# RPKM correct expression values
# 7 dev stages
# Covariates: `sex`, `rin`, `ethnicity` 

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

# library(pheatmap)
# library(cowplot)
# library(reshape2)
# library(performance) # LM diagnostic plots
# library(dlookr)
# library(explore)
# library(car)
# library(lmtest) 
# library(preprocessCore) # Quantile normalisation
 library(bestNormalize) # 
# library(MASS) # LAD model
# library(flextable)
# library(tidyverse)

##  Load Data  ------------------------------------------------------------------------
expr_tbl <- readRDS(paste0(bulk_dir, 'dystonia_brain_span_clean_tbl.rds')) %>%
  dplyr::select(donor, region, dev_stage, ethnicity, rin, sex, sort(names(.)[7:dim(.)[2]]))

##  Step 1  ---------------------------------------------------------------------------
all_yeojohn_lm_lst <- list()
all_yeojohn_lm_res_lst <- list()

# yeojohn transform 
exp_yeojohn_tbl <- expr_tbl |> 
  mutate(across(ACTB:YY1, ~ yeojohnson(.x)$x.t))

for (gene in c(dystonia_genes)) {
  
  # Run linear model for dystonia genes
  linear_model <- lm(paste0(gene, '~', 'rin + ethnicity + sex'), data = exp_yeojohn_tbl)
  lm_residuals <- linear_model$residuals
  linear_model <- broom::tidy(linear_model)
  
  # Add gene summary (with name) to list
  all_yeojohn_lm_lst[[paste0(gene, '_lm_all')]] <- linear_model
  all_yeojohn_lm_res_lst[[paste0(gene, '_res')]] <- lm_residuals
  
} 

# Combine residuals
step1_yeoJohn_tbl <- expr_tbl |>
  dplyr::select(dev_stage, region, donor) |>
  bind_cols(all_yeojohn_lm_res_lst |> as_tibble()) |> # Add step 1 residuals to tbl
  pivot_longer(ADCY5_res:VPS16_res, names_to = "gene", values_to = "step1_residuals") 
glimpse(step1_yeoJohn_tbl)

step1_yeoJohn_tbl <- model.matrix(~ dev_stage - 1, data = step1_yeoJohn_tbl) |>
  as_tibble() |>
  rename_with(~ str_replace(., 'dev_stage', '')) |>
  bind_cols(step1_yeoJohn_tbl) 

# Run model for each Dev_stage
dev_stage_lm_lst <- list()
dev_stage_lm_p <- list()

for (dev_stage in unique(step1_yeoJohn_tbl$dev_stage)) {
  
  linear_model <- lm(paste0('step1_residuals', '~', dev_stage, '+ region'), step1_yeoJohn_tbl)
  
  lm_pval <- summary(linear_model)$coefficients[2,]
  linear_model <- broom::tidy(linear_model)
  
  # Add gene summary (with name) to list
  dev_stage_lm_lst[[paste0(dev_stage, '_lm')]] <- linear_model
  dev_stage_lm_p[[paste0(dev_stage, '_P')]] <- lm_pval
  
}

# Pull out region line of lms and retain item name
dev_stage_lms <- lapply(names(dev_stage_lm_lst), function(model_name) {
  
  dev_stage <- str_split_i(model_name, '_', 1)
  
  df <- dev_stage_lm_lst[[model_name]]
  intercept_data <- df |> 
    filter(term == dev_stage)
  
  # Replace '(Intercept)' with the model name
  intercept_data <- intercept_data %>%
    mutate(term = model_name) 
  
  return(intercept_data)})

plt_data <- bind_rows(dev_stage_lms) |>
  mutate(dev_stage = str_replace(term, '_lm', '')) |>
  mutate(P_fdr = p.adjust(p.value, method = "fdr")) |>
  mutate(P_fdr_bool = ifelse(P_fdr < 0.05, TRUE, FALSE)) |>
  mutate(neglog10p = -log10(P_fdr)) |>
  mutate(neglog10p = ifelse(estimate < 0, -neglog10p, neglog10p)) |>
  mutate(dev_stage = factor(dev_stage, levels = rev(dev_levels))) |>
  dplyr::select(-term) |>
  relocate(dev_stage) |>
  mutate(dev_stage = factor(if_else(dev_stage %in% c("EarlyFetal", "MidFetal", "LateFetal"), 
                                    str_replace(dev_stage, "([A-Za-z]+)(Fetal)", "\\1-Fetal"), 
                                    dev_stage), 
                            levels = c("Adulthood", "Adolescence", "Childhood", "Infancy", 
                                       "Late-Fetal", "Mid-Fetal", "Early-Fetal")))

dev_stage_plt <- plt_data |>
  ggplot(aes(x = dev_stage, y = neglog10p)) +
  geom_col(aes(fill = P_fdr_bool), color = 'black') +
  coord_flip() +  # Flip the coordinates for better readability
  labs(x = "", y = expression("-log"[10] * "(P"[adj] * ") sgn(" * beta * ")")) +
  
  theme_minimal() +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  theme(panel.grid.minor.x = element_blank(),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 13),
        axis.text = element_text(size = 13, colour = 'black'),
        axis.title.x = element_text(size = 14, vjust = -1),
        ) +
  # scale_fill_manual(values = "#00BFC4") + # Need to add this as everything is sig.
  ylim(-140, 140) +
  Seurat::NoLegend()


print(dev_stage_plt)

bind_rows(plt_data) |>
  write_csv(paste0(table_dir, 'brain_span_lm.csv'))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

