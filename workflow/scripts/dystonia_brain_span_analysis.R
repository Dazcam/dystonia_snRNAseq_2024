#BrainSpan Analysis 08_02_2023
#Data downloaded from website
#Calculating adjusted expression levels based on Clifton et al Translational Psychiatry 2019
#Tables generated manually from downloaded data
#Using limma

#Install
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
#install.packages("dplyr") # Loaded with tidyverse
install.packages("openxlsx", dependencies = TRUE) 
install.packages("magrittr")
install.packages("broom", dependencies = TRUE)
install.packages("tidyverse", dependencies = TRUE)
#install.packages("ggplot") # Loaded with tidyverse
install.packages("ggpubr")
install.packages("vctrs")


library(limma)
library(dplyr)
library(openxlsx)
library(writexl)
library(magrittr)
library(tidyverse)
library(ggpubr)
library(broom)
library(ggplot2)
library(gridExtra)
library(psych)
library(irr)
library(plotrix)
library(ggsignif)
library(ggpubr)
library(viridisLite)
library(rlang)
library(RColorBrewer)
library(agricolae)

#Step 1 - Initial QC and file set up
#Manual selection of Dystonia Genes and manual set up file with expression data where RIN, Sex and Ethnicity also available
#Resultant file = RPKM_Adjustment_2.0
#First remove samples with RIN <7.0
RPKM_Adjustment_3.1 <- RPKM_Adjustment_2.1 %>% filter(RIN>=7)
write.xlsx(RPKM_Adjustment_3.1,'~/Desktop/RPKM_Adjustment_3.1.xlsx')
#Re-import RPKM_Adjustment_3.0 file having checked it

#Step 2 - Sort out Age & Developmental Stage
#Manually convert age to days
# block ages into developmental stages
RPKM_Adjustment_4.0$Dev_stage <- sapply(RPKM_Adjustment_4.0$Age.days, function(x) {
  if (x < -196) {
    "EarlyFetal"
  } else if (x >= -196 & x < -154) {
    "EarlyMidfetal"
  } else if (x >= -154 & x < -147) {
    "Midfetal"
  } else if (x >= -147 & x < -98) {
    "LateMidfetal"
  } else if (x >= -98 & x < 0) {
    "LateFetal"
  } else if (x >= 0 & x < 180) {
    "EarlyInfancy"
  } else if (x >= 180 & x < 366) {
    "LateInfancy"
  } else if (x >= 366 & x < 2190) {
    "EarlyChildhood"
  } else if (x >= 2190 & x < 4745) {
    "LateChildhood"
  } else if (x >= 4745 & x < 7300) {
    "Adolescence"
  } else if (x >= 7300 & x < 10950) {
    "YoungAdulthood"
  } else if (x >= 10950 & x < 21900) {
    "MidAdulthood"
  } else {
    print(x)
    print("No developmental stage")
    "NoStage"
  }
})

RPKM_Adjustment_5.0 <- data.frame(RPKM_Adjustment_4.0, RPKM_Adjustment_4.0$Dev_stage)
write.xlsx(RPKM_Adjustment_5.0,'~/Desktop/RPKM_Adjustment_5.0.xlsx')

#Step 3 - RPKM expression values adjusted to control for RIN, ethnicity, gender and age
#linear regression model of expression with above factors as independent variables
#set levels for individual variables
#Difference between 5.0 and 6.0 is one gene added at end, and multiple others removed due to clinical phenotype
#Difference between 6.0 and 7.0 is baseline developmental stage added and MLL2 genes removed as duplicate (with MLL4) for KMT2B

#Whole Brain Analysis
RPKM_Adjustment_7.1$Ethnicity <- factor(RPKM_Adjustment_7.1$Ethnicity, 
                                        levels=c("E", "As", "A", "A_E", "H"))
RPKM_Adjustment_7.1$Sex <- factor(RPKM_Adjustment_7.1$Sex, 
                                  levels=c("M", "F"))
RPKM_Adjustment_7.1$Dev_stage <- factor(RPKM_Adjustment_7.1$Dev_stage,
                                        levels = c("Early Fetal",
                                                   "Early Midfetal", 
                                                   "Midfetal",
                                                   "Late Midfetal", 
                                                   "Late Fetal",
                                                   "Early Infancy", 
                                                   "Late Infancy",
                                                   "Early Childhood", 
                                                   "Late Childhood", 
                                                   "Adolescence", 
                                                   "Early Adulthood", 
                                                   "Mid Adulthood"))
table(RPKM_Adjustment_7.1$Dev_stage)
table(RPKM_Adjustment_7.1$Region.2)

#Pre-frontal Cortex
Prefrontal.data <-RPKM_Adjustment_7.0 %>% filter(RPKM_Adjustment_7.0$Region.2=="Prefrontal Cortex")
write.xlsx(Prefrontal.data,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal.xlsx')
#Primary Motor-Sensory Cortex
MS.data <-RPKM_Adjustment_7.0 %>% filter(RPKM_Adjustment_7.0$Region.2=="Primary Motor-Sensory Cortex")
write.xlsx(MS.data,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS.xlsx')
#Striatum
Striatum.data <-RPKM_Adjustment_7.0 %>% filter(RPKM_Adjustment_7.0$Region.2=="Striatum")
write.xlsx(Striatum.data,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum.xlsx')
#Thalamus
Thalamus.data <-RPKM_Adjustment_7.0 %>% filter(RPKM_Adjustment_7.0$Region.2=="Thalamus")
write.xlsx(Thalamus.data,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus.xlsx')
#Cerebellum
Cerebellum.data <- RPKM_Adjustment_7.0 %>% filter(RPKM_Adjustment_7.0$Region.2=="Cerebellum")
write.xlsx(Cerebellum.data,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum.xlsx')
#Non-Prefrontal Cortex
NPF.data <-RPKM_Adjustment_7.0 %>% filter(RPKM_Adjustment_7.0$Region.2=="Non-Prefrontal Cortex")
write.xlsx(NPF.data,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NPF.xlsx')
#Hippocamppus
Hippocampus.data <- RPKM_Adjustment_7.0 %>% filter(RPKM_Adjustment_7.0$Region.2=="Hippocampus")
write.xlsx(Hippocampus.data,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus.xlsx')
#Other Subcortical Regions
OSR.data <- RPKM_Adjustment_7.0 %>% filter(RPKM_Adjustment_7.0$Region.2=="Other Subcortical Regions")
write.xlsx(OSR.data,'~/Desktop/Embryonic.Expression.Work/Brain_Span/OSR.xlsx')

#Manually add in baseline '0' data for re-leveling in regression model. File created 'Prefrontal.csv'
#Manually import 'Prefrontal.csv'
Prefrontal$Dev_stage <- factor(Prefrontal$Dev_stage,
                                        levels = c("Baseline",
                                                  "Early Fetal",
                                                   "Early Midfetal", 
                                                   "Midfetal",
                                                   "Late Midfetal", 
                                                   "Late Fetal",
                                                   "Early Infancy", 
                                                   "Late Infancy",
                                                   "Early Childhood", 
                                                   "Late Childhood", 
                                                   "Adolescence", 
                                                   "Early Adulthood", 
                                                   "Mid Adulthood"))
Prefrontal$Dev_stage <- relevel(factor(Prefrontal$Dev_stage), ref ="Baseline")

### -------  SUGGESTIONS TO SIMPLIFY CODE  ---------
# 1. Create a named list to collect all the lm output in a single object
# 2. Run all you analyses in a loop
# 3. Output the named list to a single xlsx file, one tab per gene
# 4. For now we'll just build a seprate loop for each region, later
#    we can try and nest two loops to pick up all genes across all regions.

# List genes of intreset taken from line 1757
genes <- c("EIF2AK2","SLC2A1","HPCA","MECR","SPR","PRKRA","PNKD",
           "COL6A3","ADCY5","SQSTM1","SLC6A3","SGCE","RELN","THAP1",
           "TOR1A", "CIZ1", "KCNMA1","DNAJC12","SLC18A2","ANO3","TH",
           "NKX2.1","GCH1","PRRT2","GNAO1","GNAL","TUBB4A","MLL4","ATP1A3",
           "VPS16", "KCTD17","TAF1") 

# Initialise list for lm output
gene_list <- list()

for (gene in c(genes)) {
  

  # Run linear model for dystonia genes
  linear_model <- lm(gene ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
  summary(linear_model)
  linear_model <- broom::tidy(linear_model)

  # Add gene summary (with name) to list
  gene_list[[paste0(gene, '.PF')]] <- linear_model

} 

# Check list
gene_list

# Export gene list to file
openxlsx::write.xlsx(gene_list, "~/Desktop/prefrontal.xlsx")


##. -------------------------

#Each gene in turn - completely failed to generate a loop!
EIF2AK2.PF <- lm(EIF2AK2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(EIF2AK2.PF)
EIF2AK2.PF <- tidy(EIF2AK2.PF)
write.xlsx(EIF2AK2.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/EIF2AK2.PF.xlsx')

SLC2A1.PF <- lm(SLC2A1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(SLC2A1.PF)
SLC2A1.PF <- tidy(SLC2A1.PF)
write.xlsx(SLC2A1.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/SLC2A1.PF.xlsx')

HPCA.PF <- lm(HPCA ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(HPCA.PF)
HPCA.PF <- tidy(HPCA.PF)
write.xlsx(HPCA.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/HPCA.PF.xlsx')

MECR.PF <- lm(MECR ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(MECR.PF)
MECR.PF <- tidy(MECR.PF)
write.xlsx(MECR.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/MECR.PF.xlsx')

SPR.PF <- lm(SPR ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(SPR.PF)
SPR.PF <- tidy(SPR.PF)
write.xlsx(SPR.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/SPR.PF.xlsx')

PRKRA.PF <- lm(PRKRA ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(PRKRA.PF)
PRKRA.PF <- tidy(PRKRA.PF)
write.xlsx(PRKRA.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/PRKRA.PF.xlsx')

PNKD.PF <- lm(PNKD ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(PNKD.PF)
PNKD.PF <- tidy(PNKD.PF)
write.xlsx(PNKD.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/PNKD.PF.xlsx')

COL6A3.PF <- lm(COL6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(COL6A3.PF)
COL6A3.PF <- tidy(COL6A3.PF)
write.xlsx(COL6A3.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/COL6A3.PF.xlsx')

ADCY5.PF <- lm(ADCY5 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(ADCY5.PF)
ADCY5.PF <- tidy(ADCY5.PF)
write.xlsx(ADCY5.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/ADCY5.PF.xlsx')

SQSTM1.PF <- lm(SQSTM1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(SQSTM1.PF)
SQSTM1.PF <- tidy(SQSTM1.PF)
write.xlsx(SQSTM1.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/SQSTM1.PF.xlsx')

SLC6A3.PF <- lm(SLC6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(SLC6A3.PF)
SLC6A3.PF <- tidy(SLC6A3.PF)
write.xlsx(SLC6A3.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/SLC6A3.PF.xlsx')

SGCE.PF <- lm(SGCE ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(SGCE.PF)
SGCE.PF <- tidy(SGCE.PF)
write.xlsx(SGCE.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/SGCE.PF.xlsx')

RELN.PF <- lm(RELN ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(RELN.PF)
RELN.PF <- tidy(RELN.PF)
write.xlsx(RELN.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/RELN.PF.xlsx')

THAP1.PF <- lm(THAP1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(THAP1.PF)
THAP1.PF <- tidy(THAP1.PF)
write.xlsx(THAP1.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/THAP1.PF.xlsx')

TOR1A.PF <- lm(TOR1A ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(TOR1A.PF)
TOR1A.PF <- tidy(TOR1A.PF)
write.xlsx(TOR1A.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/TOR1A.PF.xlsx')

CIZ1.PF <- lm(CIZ1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(CIZ1.PF)
CIZ1.PF <- tidy(CIZ1.PF)
write.xlsx(CIZ1.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/CIZ1.PF.xlsx')

KCNMA1.PF <- lm(KCNMA1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(KCNMA1.PF)
KCNMA1.PF <- tidy(KCNMA1.PF)
write.xlsx(KCNMA1.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/KCNMA1.PF.xlsx')

DNAJC12.PF <- lm(DNAJC12 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(DNAJC12.PF)
DNAJC12.PF <- tidy(DNAJC12.PF)
write.xlsx(DNAJC12.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/DNAJC12.PF.xlsx')

SLC18A2.PF <- lm(SLC18A2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(SLC18A2.PF)
SLC18A2.PF <- tidy(SLC18A2.PF)
write.xlsx(SLC18A2.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/SLC18A2.PF.xlsx')

ANO3.PF <- lm(ANO3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(ANO3.PF)
ANO3.PF <- tidy(ANO3.PF)
write.xlsx(ANO3.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/ANO3.PF.xlsx')

TH.PF <- lm(TH ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(TH.PF)
TH.PF <- tidy(TH.PF)
write.xlsx(TH.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/TH.PF.xlsx')

NKX2.1.PF <- lm(NKX2.1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(NKX2.1.PF)
NKX2.1.PF <- tidy(NKX2.1.PF)
write.xlsx(NKX2.1.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/NKX2.1.PF.xlsx')

GCH1.PF <- lm(GCH1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(GCH1.PF)
GCH1.PF <- tidy(GCH1.PF)
write.xlsx(GCH1.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/GCH1.PF.xlsx')

PRRT2.PF <- lm(PRRT2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(PRRT2.PF)
PRRT2.PF <- tidy(PRRT2.PF)
write.xlsx(PRRT2.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/PRRT2.PF.xlsx')

GNAO1.PF <- lm(GNAO1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(GNAO1.PF)
GNAO1.PF <- tidy(GNAO1.PF)
write.xlsx(GNAO1.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/GNAO1.PF.xlsx')

GNAL.PF <- lm(GNAL ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(GNAL.PF)
GNAL.PF <- tidy(GNAL.PF)
write.xlsx(GNAL.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/GNAL.PF.xlsx')

TUBB4A.PF <- lm(TUBB4A ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(TUBB4A.PF)
TUBB4A.PF <- tidy(TUBB4A.PF)
write.xlsx(TUBB4A.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/TUBB4A.PF.xlsx')

MLL4.PF <- lm(MLL4 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(MLL4.PF)
MLL4.PF <- tidy(MLL4.PF)
write.xlsx(MLL4.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/MLL4.PF.xlsx')

ATP1A3.PF <- lm(ATP1A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(ATP1A3.PF)
ATP1A3.PF <- tidy(ATP1A3.PF)
write.xlsx(ATP1A3.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/ATP1A3.PF.xlsx')

VPS16.PF <- lm(VPS16 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(VPS16.PF)
VPS16.PF <- tidy(VPS16.PF)
write.xlsx(VPS16.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/VPS16.PF.xlsx')

KCTD17.PF <- lm(KCTD17 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(KCTD17.PF)
KCTD17.PF <- tidy(KCTD17.PF)
write.xlsx(KCTD17.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/KCTD17.PF.xlsx')

TAF1.PF <- lm(TAF1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Prefrontal)
summary(TAF1.PF)
TAF1.PF <- tidy(TAF1.PF)
write.xlsx(TAF1.PF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Prefrontal/TAF1.PF.xlsx')


#################################
#Manually import 'MS.csv'
MS$Dev_stage <- factor(MS$Dev_stage,
                               levels = c("Baseline",
                                          "Early Fetal",
                                          "Early Midfetal", 
                                          "Midfetal",
                                          "Late Midfetal", 
                                          "Late Fetal",
                                          "Early Infancy", 
                                          "Late Infancy",
                                          "Early Childhood", 
                                          "Late Childhood", 
                                          "Adolescence", 
                                          "Early Adulthood", 
                                          "Mid Adulthood"))
MS$Dev_stage <- relevel(factor(MS$Dev_stage), ref ="Baseline")
#Each gene in turn - completely failed to generate a loop!
EIF2AK2.MS <- lm(EIF2AK2 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(EIF2AK2.MS)
EIF2AK2.MS <- tidy(EIF2AK2.MS)
write.xlsx(EIF2AK2.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/EIF2AK2.MS.xlsx')

SLC2A1.MS <- lm(SLC2A1 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(SLC2A1.MS)
SLC2A1.MS <- tidy(SLC2A1.MS)
write.xlsx(SLC2A1.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/SLC2A1.MS.xlsx')

HPCA.MS <- lm(HPCA ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(HPCA.MS)
HPCA.MS <- tidy(HPCA.MS)
write.xlsx(HPCA.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/HPCA.MS.xlsx')

MECR.MS <- lm(MECR ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(MECR.MS)
MECR.MS <- tidy(MECR.MS)
write.xlsx(MECR.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/MECR.MS.xlsx')

SPR.MS <- lm(SPR ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(SPR.MS)
SPR.MS <- tidy(SPR.MS)
write.xlsx(SPR.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/SPR.MS.xlsx')

PRKRA.MS <- lm(PRKRA ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(PRKRA.MS)
PRKRA.MS <- tidy(PRKRA.MS)
write.xlsx(PRKRA.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/PRKRA.MS.xlsx')

PNKD.MS <- lm(PNKD ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(PNKD.MS)
PNKD.MS <- tidy(PNKD.MS)
write.xlsx(PNKD.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/PNKD.MS.xlsx')

COL6A3.MS <- lm(COL6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(COL6A3.MS)
COL6A3.MS <- tidy(COL6A3.MS)
write.xlsx(COL6A3.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/COL6A3.MS.xlsx')

ADCY5.MS <- lm(ADCY5 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(ADCY5.MS)
ADCY5.MS <- tidy(ADCY5.MS)
write.xlsx(ADCY5.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/ADCY5.MS.xlsx')

SQSTM1.MS <- lm(SQSTM1 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(SQSTM1.MS)
SQSTM1.MS <- tidy(SQSTM1.MS)
write.xlsx(SQSTM1.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/SQSTM1.MS.xlsx')

SLC6A3.MS <- lm(SLC6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(SLC6A3.MS)
SLC6A3.MS <- tidy(SLC6A3.MS)
write.xlsx(SLC6A3.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/SLC6A3.MS.xlsx')

SGCE.MS <- lm(SGCE ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(SGCE.MS)
SGCE.MS <- tidy(SGCE.MS)
write.xlsx(SGCE.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/SGCE.MS.xlsx')

RELN.MS <- lm(RELN ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(RELN.MS)
RELN.MS <- tidy(RELN.MS)
write.xlsx(RELN.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/RELN.MS.xlsx')

THAP1.MS <- lm(THAP1 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(THAP1.MS)
THAP1.MS <- tidy(THAP1.MS)
write.xlsx(THAP1.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/THAP1.MS.xlsx')

TOR1A.MS <- lm(TOR1A ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(TOR1A.MS)
TOR1A.MS <- tidy(TOR1A.MS)
write.xlsx(TOR1A.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/TOR1A.MS.xlsx')

CIZ1.MS <- lm(CIZ1 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(CIZ1.MS)
CIZ1.MS <- tidy(CIZ1.MS)
write.xlsx(CIZ1.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/CIZ1.MS.xlsx')

KCNMA1.MS <- lm(KCNMA1 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(KCNMA1.MS)
KCNMA1.MS <- tidy(KCNMA1.MS)
write.xlsx(KCNMA1.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/KCNMA1.MS.xlsx')

DNAJC12.MS <- lm(DNAJC12 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(DNAJC12.MS)
DNAJC12.MS <- tidy(DNAJC12.MS)
write.xlsx(DNAJC12.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/DNAJC12.MS.xlsx')

SLC18A2.MS <- lm(SLC18A2 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(SLC18A2.MS)
SLC18A2.MS <- tidy(SLC18A2.MS)
write.xlsx(SLC18A2.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/SLC18A2.MS.xlsx')

ANO3.MS <- lm(ANO3 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(ANO3.MS)
ANO3.MS <- tidy(ANO3.MS)
write.xlsx(ANO3.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/ANO3.MS.xlsx')

TH.MS <- lm(TH ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(TH.MS)
TH.MS <- tidy(TH.MS)
write.xlsx(TH.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/TH.MS.xlsx')

NKX2.1.MS <- lm(NKX2.1 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(NKX2.1.MS)
NKX2.1.MS <- tidy(NKX2.1.MS)
write.xlsx(NKX2.1.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/NKX2.1.MS.xlsx')

GCH1.MS <- lm(GCH1 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(GCH1.MS)
GCH1.MS <- tidy(GCH1.MS)
write.xlsx(GCH1.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/GCH1.MS.xlsx')

PRRT2.MS <- lm(PRRT2 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(PRRT2.MS)
PRRT2.MS <- tidy(PRRT2.MS)
write.xlsx(PRRT2.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/PRRT2.MS.xlsx')

GNAO1.MS <- lm(GNAO1 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(GNAO1.MS)
GNAO1.MS <- tidy(GNAO1.MS)
write.xlsx(GNAO1.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/GNAO1.MS.xlsx')

GNAL.MS <- lm(GNAL ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(GNAL.MS)
GNAL.MS <- tidy(GNAL.MS)
write.xlsx(GNAL.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/GNAL.MS.xlsx')

TUBB4A.MS <- lm(TUBB4A ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(TUBB4A.MS)
TUBB4A.MS <- tidy(TUBB4A.MS)
write.xlsx(TUBB4A.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/TUBB4A.MS.xlsx')

MLL4.MS <- lm(MLL4 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(MLL4.MS)
MLL4.MS <- tidy(MLL4.MS)
write.xlsx(MLL4.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/MLL4.MS.xlsx')

ATP1A3.MS <- lm(ATP1A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(ATP1A3.MS)
ATP1A3.MS <- tidy(ATP1A3.MS)
write.xlsx(ATP1A3.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/ATP1A3.MS.xlsx')

VPS16.MS <- lm(VPS16 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(VPS16.MS)
VPS16.MS <- tidy(VPS16.MS)
write.xlsx(VPS16.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/VPS16.MS.xlsx')

KCTD17.MS <- lm(KCTD17 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(KCTD17.MS)
KCTD17.MS <- tidy(KCTD17.MS)
write.xlsx(KCTD17.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/KCTD17.MS.xlsx')

TAF1.MS <- lm(TAF1 ~ RIN + Ethnicity + Sex + Dev_stage, data = MS)
summary(TAF1.MS)
TAF1.MS <- tidy(TAF1.MS)
write.xlsx(TAF1.MS,'~/Desktop/Embryonic.Expression.Work/Brain_Span/MS/TAF1.MS.xlsx')

########################################
#Manually import 'Thalamus.csv'
Thalamus$Dev_stage <- factor(Thalamus$Dev_stage,
                       levels = c("Baseline",
                                  "Early Fetal",
                                  "Early Midfetal", 
                                  "Midfetal",
                                  "Late Midfetal", 
                                  "Late Fetal",
                                  "Early Infancy", 
                                  "Late Infancy",
                                  "Early Childhood", 
                                  "Late Childhood", 
                                  "Adolescence", 
                                  "Early Adulthood", 
                                  "Mid Adulthood"))
Thalamus$Dev_stage <- relevel(factor(Thalamus$Dev_stage), ref ="Baseline")
#Each gene in turn - completely failed to generate a loop!
EIF2AK2.Thalamus <- lm(EIF2AK2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(EIF2AK2.Thalamus)
EIF2AK2.Thalamus <- tidy(EIF2AK2.Thalamus)
write.xlsx(EIF2AK2.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/EIF2AK2.Thalamus.xlsx')

SLC2A1.Thalamus <- lm(SLC2A1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(SLC2A1.Thalamus)
SLC2A1.Thalamus <- tidy(SLC2A1.Thalamus)
write.xlsx(SLC2A1.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/SLC2A1.Thalamus.xlsx')

HPCA.Thalamus <- lm(HPCA ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(HPCA.Thalamus)
HPCA.Thalamus <- tidy(HPCA.Thalamus)
write.xlsx(HPCA.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/HPCA.Thalamus.xlsx')

MECR.Thalamus <- lm(MECR ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(MECR.Thalamus)
MECR.Thalamus <- tidy(MECR.Thalamus)
write.xlsx(MECR.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/MECR.Thalamus.xlsx')

SPR.Thalamus <- lm(SPR ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(SPR.Thalamus)
SPR.Thalamus <- tidy(SPR.Thalamus)
write.xlsx(SPR.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/SPR.Thalamus.xlsx')

PRKRA.Thalamus <- lm(PRKRA ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(PRKRA.Thalamus)
PRKRA.Thalamus <- tidy(PRKRA.Thalamus)
write.xlsx(PRKRA.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/PRKRA.Thalamus.xlsx')

PNKD.Thalamus <- lm(PNKD ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(PNKD.Thalamus)
PNKD.Thalamus <- tidy(PNKD.Thalamus)
write.xlsx(PNKD.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/PNKD.Thalamus.xlsx')

COL6A3.Thalamus <- lm(COL6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(COL6A3.Thalamus)
COL6A3.Thalamus <- tidy(COL6A3.Thalamus)
write.xlsx(COL6A3.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/COL6A3.Thalamus.xlsx')

ADCY5.Thalamus <- lm(ADCY5 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(ADCY5.Thalamus)
ADCY5.Thalamus <- tidy(ADCY5.Thalamus)
write.xlsx(ADCY5.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/ADCY5.Thalamus.xlsx')

SQSTM1.Thalamus <- lm(SQSTM1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(SQSTM1.Thalamus)
SQSTM1.Thalamus <- tidy(SQSTM1.Thalamus)
write.xlsx(SQSTM1.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/SQSTM1.Thalamus.xlsx')

SLC6A3.Thalamus <- lm(SLC6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(SLC6A3.Thalamus)
SLC6A3.Thalamus <- tidy(SLC6A3.Thalamus)
write.xlsx(SLC6A3.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/SLC6A3.Thalamus.xlsx')

SGCE.Thalamus <- lm(SGCE ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(SGCE.Thalamus)
SGCE.Thalamus <- tidy(SGCE.Thalamus)
write.xlsx(SGCE.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/SGCE.Thalamus.xlsx')

RELN.Thalamus <- lm(RELN ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(RELN.Thalamus)
RELN.Thalamus <- tidy(RELN.Thalamus)
write.xlsx(RELN.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/RELN.Thalamus.xlsx')

THAP1.Thalamus <- lm(THAP1 ~ RIN + Ethnicity + Sex + Dev_stage, data =Thalamus)
summary(THAP1.Thalamus)
THAP1.Thalamus <- tidy(THAP1.Thalamus)
write.xlsx(THAP1.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/THAP1.Thalamus.xlsx')

TOR1A.Thalamus <- lm(TOR1A ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(TOR1A.Thalamus)
TOR1A.Thalamus <- tidy(TOR1A.Thalamus)
write.xlsx(TOR1A.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/TOR1A.Thalamus.xlsx')

CIZ1.Thalamus <- lm(CIZ1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(CIZ1.Thalamus)
CIZ1.Thalamus <- tidy(CIZ1.Thalamus)
write.xlsx(CIZ1.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/CIZ1.Thalamus.xlsx')

KCNMA1.Thalamus <- lm(KCNMA1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(KCNMA1.Thalamus)
KCNMA1.Thalamus <- tidy(KCNMA1.Thalamus)
write.xlsx(KCNMA1.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/KCNMA1.Thalamus.xlsx')

DNAJC12.Thalamus <- lm(DNAJC12 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(DNAJC12.Thalamus)
DNAJC12.Thalamus <- tidy(DNAJC12.Thalamus)
write.xlsx(DNAJC12.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/DNAJC12.Thalamus.xlsx')

SLC18A2.Thalamus <- lm(SLC18A2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(SLC18A2.Thalamus)
SLC18A2.Thalamus <- tidy(SLC18A2.Thalamus)
write.xlsx(SLC18A2.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/SLC18A2.Thalamus.xlsx')

ANO3.Thalamus <- lm(ANO3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(ANO3.Thalamus)
ANO3.Thalamus <- tidy(ANO3.Thalamus)
write.xlsx(ANO3.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/ANO3.Thalamus.xlsx')

TH.Thalamus <- lm(TH ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(TH.Thalamus)
TH.Thalamus <- tidy(TH.Thalamus)
write.xlsx(TH.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/TH.Thalamus.xlsx')

NKX2.1.Thalamus <- lm(NKX2.1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(NKX2.1.Thalamus)
NKX2.1.Thalamus <- tidy(NKX2.1.Thalamus)
write.xlsx(NKX2.1.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/NKX2.1.Thalamus.xlsx')

GCH1.Thalamus <- lm(GCH1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(GCH1.Thalamus)
GCH1.Thalamus <- tidy(GCH1.Thalamus)
write.xlsx(GCH1.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/GCH1.Thalamus.xlsx')

PRRT2.Thalamus <- lm(PRRT2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(PRRT2.Thalamus)
PRRT2.Thalamus <- tidy(PRRT2.Thalamus)
write.xlsx(PRRT2.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/PRRT2.Thalamus.xlsx')

GNAO1.Thalamus <- lm(GNAO1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(GNAO1.Thalamus)
GNAO1.Thalamus <- tidy(GNAO1.Thalamus)
write.xlsx(GNAO1.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/GNAO1.Thalamus.xlsx')

GNAL.Thalamus <- lm(GNAL ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(GNAL.Thalamus)
GNAL.Thalamus <- tidy(GNAL.Thalamus)
write.xlsx(GNAL.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/GNAL.Thalamus.xlsx')

TUBB4A.Thalamus <- lm(TUBB4A ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(TUBB4A.Thalamus)
TUBB4A.Thalamus <- tidy(TUBB4A.Thalamus)
write.xlsx(TUBB4A.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/TUBB4A.Thalamus.xlsx')

MLL4.Thalamus <- lm(MLL4 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(MLL4.Thalamus)
MLL4.Thalamus <- tidy(MLL4.Thalamus)
write.xlsx(MLL4.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/MLL4.Thalamus.xlsx')

ATP1A3.Thalamus <- lm(ATP1A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(ATP1A3.Thalamus)
ATP1A3.Thalamus <- tidy(ATP1A3.Thalamus)
write.xlsx(ATP1A3.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/ATP1A3.Thalamus.xlsx')

VPS16.Thalamus <- lm(VPS16 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(VPS16.Thalamus)
VPS16.Thalamus <- tidy(VPS16.Thalamus)
write.xlsx(VPS16.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/VPS16.Thalamus.xlsx')

KCTD17.Thalamus <- lm(KCTD17 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(KCTD17.Thalamus)
KCTD17.Thalamus <- tidy(KCTD17.Thalamus)
write.xlsx(KCTD17.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/KCTD17.Thalamus.xlsx')

TAF1.Thalamus <- lm(TAF1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Thalamus)
summary(TAF1.Thalamus)
TAF1.Thalamus <- tidy(TAF1.Thalamus)
write.xlsx(TAF1.Thalamus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Thalamus/TAF1.Thalamus.xlsx')

#####################
#Manually import 'Striatum.csv'
Striatum$Dev_stage <- factor(Striatum$Dev_stage,
                             levels = c("Baseline",
                                        "Early Fetal",
                                        "Early Midfetal", 
                                        "Midfetal",
                                        "Late Midfetal", 
                                        "Late Fetal",
                                        "Early Infancy", 
                                        "Late Infancy",
                                        "Early Childhood", 
                                        "Late Childhood", 
                                        "Adolescence", 
                                        "Early Adulthood", 
                                        "Mid Adulthood"))
Striatum$Dev_stage <- relevel(factor(Striatum$Dev_stage), ref ="Baseline")
#Each gene in turn - completely failed to generate a loop!
EIF2AK2.Striatum <- lm(EIF2AK2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(EIF2AK2.Striatum)
EIF2AK2.Striatum <- tidy(EIF2AK2.Striatum)
write.xlsx(EIF2AK2.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/EIF2AK2.Striatum.xlsx')

SLC2A1.Striatum <- lm(SLC2A1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(SLC2A1.Striatum)
SLC2A1.Striatum <- tidy(SLC2A1.Striatum)
write.xlsx(SLC2A1.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/SLC2A1.Striatum.xlsx')

HPCA.Striatum <- lm(HPCA ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(HPCA.Striatum)
HPCA.Striatum <- tidy(HPCA.Striatum)
write.xlsx(HPCA.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/HPCA.Striatum.xlsx')

MECR.Striatum <- lm(MECR ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(MECR.Striatum)
MECR.Striatum <- tidy(MECR.Striatum)
write.xlsx(MECR.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/MECR.Striatum.xlsx')

SPR.Striatum <- lm(SPR ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(SPR.Striatum)
SPR.Striatum <- tidy(SPR.Striatum)
write.xlsx(SPR.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/SPR.Striatum.xlsx')

PRKRA.Striatum <- lm(PRKRA ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(PRKRA.Striatum)
PRKRA.TStriatum <- tidy(PRKRA.Striatum)
write.xlsx(PRKRA.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/PRKRA.Striatum.xlsx')

PNKD.Striatum <- lm(PNKD ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(PNKD.Striatum)
PNKD.Striatum <- tidy(PNKD.Striatum)
write.xlsx(PNKD.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/PNKD.Striatum.xlsx')

COL6A3.Striatum <- lm(COL6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(COL6A3.Striatum)
COL6A3.Striatum <- tidy(COL6A3.Striatum)
write.xlsx(COL6A3.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/COL6A3.Striatum.xlsx')

ADCY5.Striatum <- lm(ADCY5 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(ADCY5.Striatum)
ADCY5.Striatum <- tidy(ADCY5.Striatum)
write.xlsx(ADCY5.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/ADCY5.Striatum.xlsx')

SQSTM1.Striatum <- lm(SQSTM1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(SQSTM1.Striatum)
SQSTM1.Striatum <- tidy(SQSTM1.Striatum)
write.xlsx(SQSTM1.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/SQSTM1.Striatum.xlsx')

SLC6A3.Striatum <- lm(SLC6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(SLC6A3.Striatum)
SLC6A3.Striatum <- tidy(SLC6A3.Striatum)
write.xlsx(SLC6A3.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/SLC6A3.Striatum.xlsx')

SGCE.Striatum <- lm(SGCE ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(SGCE.Striatum)
SGCE.Striatum <- tidy(SGCE.Striatum)
write.xlsx(SGCE.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/SGCE.Striatum.xlsx')

RELN.Striatum <- lm(RELN ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(RELN.Striatum)
RELN.Striatum <- tidy(RELN.Striatum)
write.xlsx(RELN.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/RELN.Striatum.xlsx')

THAP1.Striatum <- lm(THAP1 ~ RIN + Ethnicity + Sex + Dev_stage, data =Striatum)
summary(THAP1.Striatum)
THAP1.Striatum <- tidy(THAP1.Striatum)
write.xlsx(THAP1.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/THAP1.Striatum.xlsx')

TOR1A.Striatum <- lm(TOR1A ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(TOR1A.Striatum)
TOR1A.Striatum <- tidy(TOR1A.Striatum)
write.xlsx(TOR1A.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/TOR1A.Striatum.xlsx')

CIZ1.Striatum <- lm(CIZ1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(CIZ1.Striatum)
CIZ1.Striatum <- tidy(CIZ1.Striatum)
write.xlsx(CIZ1.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/CIZ1.Striatum.xlsx')

KCNMA1.Striatum <- lm(KCNMA1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(KCNMA1.Striatum)
KCNMA1.Striatum <- tidy(KCNMA1.Striatum)
write.xlsx(KCNMA1.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/KCNMA1.Striatum.xlsx')

DNAJC12.Striatum <- lm(DNAJC12 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(DNAJC12.Striatum)
DNAJC12.Striatum <- tidy(DNAJC12.Striatum)
write.xlsx(DNAJC12.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/DNAJC12.Striatum.xlsx')

SLC18A2.Striatum <- lm(SLC18A2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(SLC18A2.Striatum)
SLC18A2.Striatum <- tidy(SLC18A2.Striatum)
write.xlsx(SLC18A2.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/SLC18A2.Striatum.xlsx')

ANO3.Striatum <- lm(ANO3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(ANO3.Striatum)
ANO3.Striatum <- tidy(ANO3.Striatum)
write.xlsx(ANO3.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/ANO3.Striatum.xlsx')

TH.Striatum <- lm(TH ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(TH.Striatum)
TH.Striatum <- tidy(TH.Striatum)
write.xlsx(TH.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/TH.Striatum.xlsx')

NKX2.1.Striatum <- lm(NKX2.1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(NKX2.1.Striatum)
NKX2.1.Striatum <- tidy(NKX2.1.Striatum)
write.xlsx(NKX2.1.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/NKX2.1.Striatum.xlsx')

GCH1.Striatum <- lm(GCH1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(GCH1.Striatum)
GCH1.Striatum <- tidy(GCH1.Striatum)
write.xlsx(GCH1.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/GCH1.Striatum.xlsx')

PRRT2.Striatum <- lm(PRRT2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(PRRT2.Striatum)
PRRT2.Striatum<- tidy(PRRT2.Striatum)
write.xlsx(PRRT2.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/PRRT2.Striatum.xlsx')

GNAO1.Striatum <- lm(GNAO1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(GNAO1.Striatum)
GNAO1.Striatum <- tidy(GNAO1.Striatum)
write.xlsx(GNAO1.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/GNAO1.Striatum.xlsx')

GNAL.Striatum <- lm(GNAL ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(GNAL.Striatum)
GNAL.Striatum <- tidy(GNAL.Striatum)
write.xlsx(GNAL.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/GNAL.Striatum.xlsx')

TUBB4A.Striatum <- lm(TUBB4A ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(TUBB4A.Striatum)
TUBB4A.Striatum <- tidy(TUBB4A.Striatum)
write.xlsx(TUBB4A.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/TUBB4A.Striatum.xlsx')

MLL4.Striatum <- lm(MLL4 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(MLL4.Striatum)
MLL4.Striatum <- tidy(MLL4.Striatum)
write.xlsx(MLL4.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/MLL4.Striatum.xlsx')

ATP1A3.Striatum <- lm(ATP1A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(ATP1A3.Striatum)
ATP1A3.Striatum <- tidy(ATP1A3.Striatum)
write.xlsx(ATP1A3.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/ATP1A3.Striatum.xlsx')

VPS16.Striatum <- lm(VPS16 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(VPS16.Striatum)
VPS16.Striatum <- tidy(VPS16.Striatum)
write.xlsx(VPS16.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/VPS16.Striatum.xlsx')

KCTD17.Striatum <- lm(KCTD17 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(KCTD17.Striatum)
KCTD17.Striatum <- tidy(KCTD17.Striatum)
write.xlsx(KCTD17.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/KCTD17.Striatum.xlsx')

TAF1.Striatum <- lm(TAF1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Striatum)
summary(TAF1.Striatum)
TAF1.Striatum <- tidy(TAF1.Striatum)
write.xlsx(TAF1.Striatum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Striatum/TAF1.Striatum.xlsx')


################
#Manually import 'Cerebellum.csv'
Cerebellum$Dev_stage <- factor(Cerebellum$Dev_stage,
                             levels = c("Baseline",
                                        "Early Fetal",
                                        "Early Midfetal", 
                                        "Midfetal",
                                        "Late Midfetal", 
                                        "Late Fetal",
                                        "Early Infancy", 
                                        "Late Infancy",
                                        "Early Childhood", 
                                        "Late Childhood", 
                                        "Adolescence", 
                                        "Early Adulthood", 
                                        "Mid Adulthood"))
Cerebellum$Dev_stage <- relevel(factor(Cerebellum$Dev_stage), ref ="Baseline")
#Each gene in turn - completely failed to generate a loop!
EIF2AK2.Cerebellum <- lm(EIF2AK2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(EIF2AK2.Cerebellum)
EIF2AK2.Cerebellum <- tidy(EIF2AK2.Cerebellum)
write.xlsx(EIF2AK2.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/EIF2AK2.Cerebellum.xlsx')

SLC2A1.Cerebellum <- lm(SLC2A1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(SLC2A1.Cerebellum)
SLC2A1.Cerebellum <- tidy(SLC2A1.Cerebellum)
write.xlsx(SLC2A1.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/SLC2A1.Cerebellum.xlsx')

HPCA.Cerebellum <- lm(HPCA ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(HPCA.Cerebellum)
HPCA.Cerebellum <- tidy(HPCA.Cerebellum)
write.xlsx(HPCA.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/HPCA.Cerebellum.xlsx')

MECR.Cerebellum <- lm(MECR ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(MECR.Cerebellum)
MECR.Cerebellum <- tidy(MECR.Cerebellum)
write.xlsx(MECR.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/MECR.Cerebellum.xlsx')

SPR.Cerebellum <- lm(SPR ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(SPR.Cerebellum)
SPR.Cerebellum <- tidy(SPR.Cerebellum)
write.xlsx(SPR.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/SPR.Cerebellum.xlsx')

PRKRA.Cerebellum <- lm(PRKRA ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(PRKRA.Cerebellum)
PRKRA.Cerebellum <- tidy(PRKRA.Cerebellum)
write.xlsx(PRKRA.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/PRKRA.Cerebellum.xlsx')

PNKD.Cerebellum <- lm(PNKD ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(PNKD.Cerebellum)
PNKD.Cerebellum <- tidy(PNKD.Cerebellum)
write.xlsx(PNKD.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/PNKD.Cerebellum.xlsx')

COL6A3.Cerebellum <- lm(COL6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(COL6A3.Cerebellum)
COL6A3.Cerebellum <- tidy(COL6A3.Cerebellum)
write.xlsx(COL6A3.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/COL6A3.Cerebellum.xlsx')

ADCY5.Cerebellum <- lm(ADCY5 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(ADCY5.Cerebellum)
ADCY5.Cerebellum <- tidy(ADCY5.Cerebellum)
write.xlsx(ADCY5.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/ADCY5.Cerebellum.xlsx')

SQSTM1.Cerebellum <- lm(SQSTM1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(SQSTM1.Cerebellum)
SQSTM1.Cerebellum <- tidy(SQSTM1.Cerebellum)
write.xlsx(SQSTM1.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/SQSTM1.Cerebellum.xlsx')

SLC6A3.Cerebellum <- lm(SLC6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(SLC6A3.Cerebellum)
SLC6A3.Cerebellum <- tidy(SLC6A3.Cerebellum)
write.xlsx(SLC6A3.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/SLC6A3.Cerebellum.xlsx')

SGCE.Cerebellum <- lm(SGCE ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(SGCE.Cerebellum)
SGCE.Cerebellum <- tidy(SGCE.Cerebellum)
write.xlsx(SGCE.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/SGCE.Cerebellum.xlsx')

RELN.Cerebellum <- lm(RELN ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(RELN.Cerebellum)
RELN.Cerebellum <- tidy(RELN.Cerebellum)
write.xlsx(RELN.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/RELN.Cerebellum.xlsx')

THAP1.Cerebellum <- lm(THAP1 ~ RIN + Ethnicity + Sex + Dev_stage, data =Cerebellum)
summary(THAP1.Cerebellum)
THAP1.Cerebellum <- tidy(THAP1.Cerebellum)
write.xlsx(THAP1.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/THAP1.Cerebellum.xlsx')

TOR1A.Cerebellum <- lm(TOR1A ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(TOR1A.Cerebellum)
TOR1A.Cerebellum <- tidy(TOR1A.Cerebellum)
write.xlsx(TOR1A.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/TOR1A.Cerebellum.xlsx')

CIZ1.Cerebellum <- lm(CIZ1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(CIZ1.Cerebellum)
CIZ1.Cerebellum <- tidy(CIZ1.Cerebellum)
write.xlsx(CIZ1.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/CIZ1.Cerebellum.xlsx')

KCNMA1.Cerebellum <- lm(KCNMA1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(KCNMA1.Cerebellum)
KCNMA1.Cerebellum <- tidy(KCNMA1.Cerebellum)
write.xlsx(KCNMA1.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/KCNMA1.Cerebellum.xlsx')

DNAJC12.Cerebellum <- lm(DNAJC12 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(DNAJC12.Cerebellum)
DNAJC12.Cerebellum <- tidy(DNAJC12.Cerebellum)
write.xlsx(DNAJC12.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/DNAJC12.Cerebellum.xlsx')

SLC18A2.Cerebellum <- lm(SLC18A2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(SLC18A2.Cerebellum)
SLC18A2.Cerebellum <- tidy(SLC18A2.Cerebellum)
write.xlsx(SLC18A2.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/SLC18A2.Cerebellum.xlsx')

ANO3.Cerebellum <- lm(ANO3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(ANO3.Cerebellum)
ANO3.Cerebellum <- tidy(ANO3.Cerebellum)
write.xlsx(ANO3.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/ANO3.Cerebellum.xlsx')

TH.Cerebellum <- lm(TH ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(TH.Cerebellum)
TH.Cerebellum <- tidy(TH.Cerebellum)
write.xlsx(TH.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/TH.Cerebellum.xlsx')

NKX2.1.Cerebellum <- lm(NKX2.1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(NKX2.1.Cerebellum)
NKX2.1.Cerebellum <- tidy(NKX2.1.Cerebellum)
write.xlsx(NKX2.1.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/NKX2.1.Cerebellum.xlsx')

GCH1.Cerebellum <- lm(GCH1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(GCH1.Cerebellum)
GCH1.Cerebellum <- tidy(GCH1.Cerebellum)
write.xlsx(GCH1.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/GCH1.Cerebellum.xlsx')

PRRT2.Cerebellum <- lm(PRRT2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(PRRT2.Cerebellum)
PRRT2.Cerebellum<- tidy(PRRT2.Cerebellum)
write.xlsx(PRRT2.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/PRRT2.Cerebellum.xlsx')

GNAO1.Cerebellum <- lm(GNAO1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(GNAO1.Cerebellum)
GNAO1.Cerebellum <- tidy(GNAO1.Cerebellum)
write.xlsx(GNAO1.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/GNAO1.Cerebellum.xlsx')

GNAL.Cerebellum <- lm(GNAL ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(GNAL.Cerebellum)
GNAL.Cerebellum <- tidy(GNAL.Cerebellum)
write.xlsx(GNAL.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/GNAL.Cerebellum.xlsx')

TUBB4A.Cerebellum <- lm(TUBB4A ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(TUBB4A.Cerebellum)
TUBB4A.Cerebellum <- tidy(TUBB4A.Cerebellum)
write.xlsx(TUBB4A.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/TUBB4A.Cerebellum.xlsx')

MLL4.Cerebellum <- lm(MLL4 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(MLL4.Cerebellum)
MLL4.Cerebellum <- tidy(MLL4.Cerebellum)
write.xlsx(MLL4.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/MLL4.Cerebellum.xlsx')

ATP1A3.Cerebellum <- lm(ATP1A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(ATP1A3.Cerebellum)
ATP1A3.Cerebellum <- tidy(ATP1A3.Cerebellum)
write.xlsx(ATP1A3.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/ATP1A3.Cerebellum.xlsx')

VPS16.Cerebellum <- lm(VPS16 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(VPS16.Cerebellum)
VPS16.Cerebellum <- tidy(VPS16.Cerebellum)
write.xlsx(VPS16.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/VPS16.Cerebellum.xlsx')

KCTD17.Cerebellum <- lm(KCTD17 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(KCTD17.Cerebellum)
KCTD17.Cerebellum <- tidy(KCTD17.Cerebellum)
write.xlsx(KCTD17.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/KCTD17.Cerebellum.xlsx')

TAF1.Cerebellum <- lm(TAF1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Cerebellum)
summary(TAF1.Cerebellum)
TAF1.Cerebellum <- tidy(TAF1.Cerebellum)
write.xlsx(TAF1.Cerebellum,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Cerebellum/TAF1.Cerebellum.xlsx')

#Manually import NPR.csv
NPF$Dev_stage <- factor(NPF$Dev_stage,
                               levels = c("Baseline",
                                          "Early Fetal",
                                          "Early Midfetal", 
                                          "Midfetal",
                                          "Late Midfetal", 
                                          "Late Fetal",
                                          "Early Infancy", 
                                          "Late Infancy",
                                          "Early Childhood", 
                                          "Late Childhood", 
                                          "Adolescence", 
                                          "Early Adulthood", 
                                          "Mid Adulthood"))
NPF$Dev_stage <- relevel(factor(NPF$Dev_stage), ref ="Baseline")
#Each gene in turn - completely failed to generate a loop!
EIF2AK2.NPF <- lm(EIF2AK2 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(EIF2AK2.NPF)
EIF2AK2.NPF <- tidy(EIF2AK2.NPF)
write.xlsx(EIF2AK2.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/EIF2AK2.NPF.xlsx')

SLC2A1.NPF <- lm(SLC2A1 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(SLC2A1.NPF)
SLC2A1.NPF <- tidy(SLC2A1.NPF)
write.xlsx(SLC2A1.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/SLC2A1.NPF.xlsx')

HPCA.NPF <- lm(HPCA ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(HPCA.NPF)
HPCA.NPF <- tidy(HPCA.NPF)
write.xlsx(HPCA.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/HPCA.NPF.xlsx')

MECR.NPF <- lm(MECR ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(MECR.NPF)
MECR.NPF <- tidy(MECR.NPF)
write.xlsx(MECR.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/MECR.NPF.xlsx')

SPR.NPF <- lm(SPR ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(SPR.NPF)
SPR.NPF <- tidy(SPR.NPF)
write.xlsx(SPR.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/SPR.NPF.xlsx')

PRKRA.NPF <- lm(PRKRA ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(PRKRA.NPF)
PRKRA.NPF <- tidy(PRKRA.NPF)
write.xlsx(PRKRA.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/PRKRA.NPF.xlsx')

PNKD.NPF <- lm(PNKD ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(PNKD.NPF)
PNKD.NPF <- tidy(PNKD.NPF)
write.xlsx(PNKD.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/PNKD.NPF.xlsx')

COL6A3.NPF <- lm(COL6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(COL6A3.NPF)
COL6A3.NPF <- tidy(COL6A3.NPF)
write.xlsx(COL6A3.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/COL6A3.NPF.xlsx')

ADCY5.NPF <- lm(ADCY5 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(ADCY5.NPF)
ADCY5.NPF <- tidy(ADCY5.NPF)
write.xlsx(ADCY5.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/ADCY5.NPF.xlsx')

SQSTM1.NPF <- lm(SQSTM1 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(SQSTM1.NPF)
SQSTM1.NPF <- tidy(SQSTM1.NPF)
write.xlsx(SQSTM1.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/SQSTM1.NPF.xlsx')

SLC6A3.NPF <- lm(SLC6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(SLC6A3.NPF)
SLC6A3.NPF <- tidy(SLC6A3.NPF)
write.xlsx(SLC6A3.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/SLC6A3.NPF.xlsx')

SGCE.NPF <- lm(SGCE ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(SGCE.NPF)
SGCE.NPF <- tidy(SGCE.NPF)
write.xlsx(SGCE.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/SGCE.NPF.xlsx')

RELN.NPF <- lm(RELN ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(RELN.NPF)
RELN.NPF <- tidy(RELN.NPF)
write.xlsx(RELN.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/RELN.NPF.xlsx')

THAP1.NPF <- lm(THAP1 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(THAP1.NPF)
THAP1.NPF <- tidy(THAP1.NPF)
write.xlsx(THAP1.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/THAP1.NPF.xlsx')

TOR1A.NPF <- lm(TOR1A ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(TOR1A.NPF)
TOR1A.NPF <- tidy(TOR1A.NPF)
write.xlsx(TOR1A.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/TOR1A.NPF.xlsx')

CIZ1.NPF <- lm(CIZ1 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(CIZ1.NPF)
CIZ1.NPF <- tidy(CIZ1.NPF)
write.xlsx(CIZ1.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/CIZ1.NPF.xlsx')

KCNMA1.NPF <- lm(KCNMA1 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(KCNMA1.NPF)
KCNMA1.NPF <- tidy(KCNMA1.NPF)
write.xlsx(KCNMA1.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/KCNMA1.NPF.xlsx')

DNAJC12.NPF <- lm(DNAJC12 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(DNAJC12.NPF)
DNAJC12.NPF <- tidy(DNAJC12.NPF)
write.xlsx(DNAJC12.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/DNAJC12.NPF.xlsx')

SLC18A2.NPF <- lm(SLC18A2 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(SLC18A2.NPF)
SLC18A2.NPF <- tidy(SLC18A2.NPF)
write.xlsx(SLC18A2.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/SLC18A2.NPF.xlsx')

ANO3.NPF <- lm(ANO3 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(ANO3.NPF)
ANO3.NPF <- tidy(ANO3.NPF)
write.xlsx(ANO3.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/ANO3.NPF.xlsx')

TH.NPF <- lm(TH ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(TH.NPF)
TH.NPF <- tidy(TH.NPF)
write.xlsx(TH.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/TH.NPF.xlsx')

NKX2.1.NPF <- lm(NKX2.1 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(NKX2.1.NPF)
NKX2.1.NPF <- tidy(NKX2.1.NPF)
write.xlsx(NKX2.1.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/NKX2.1.NPF.xlsx')

GCH1.NPF <- lm(GCH1 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(GCH1.NPF)
GCH1.NPF <- tidy(GCH1.NPF)
write.xlsx(GCH1.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/GCH1.NPF.xlsx')

PRRT2.NPF <- lm(PRRT2 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(PRRT2.NPF)
PRRT2.NPF <- tidy(PRRT2.NPF)
write.xlsx(PRRT2.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/PRRT2.NPF.xlsx')

GNAO1.NPF <- lm(GNAO1 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(GNAO1.NPF)
GNAO1.NPF <- tidy(GNAO1.NPF)
write.xlsx(GNAO1.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/GNAO1.NPF.xlsx')

GNAL.NPF <- lm(GNAL ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(GNAL.NPF)
GNAL.NPF <- tidy(GNAL.NPF)
write.xlsx(GNAL.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/GNAL.NPF.xlsx')

TUBB4A.NPF <- lm(TUBB4A ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(TUBB4A.NPF)
TUBB4A.NPF <- tidy(TUBB4A.NPF)
write.xlsx(TUBB4A.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/TUBB4A.NPF.xlsx')

MLL4.NPF <- lm(MLL4 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(MLL4.NPF)
MLL4.NPF <- tidy(MLL4.NPF)
write.xlsx(MLL4.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/MLL4.NPF.xlsx')

ATP1A3.NPF <- lm(ATP1A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(ATP1A3.NPF)
ATP1A3.NPF <- tidy(ATP1A3.NPF)
write.xlsx(ATP1A3.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/ATP1A3.NPF.xlsx')

VPS16.NPF <- lm(VPS16 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(VPS16.NPF)
VPS16.NPF <- tidy(VPS16.NPF)
write.xlsx(VPS16.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/VPS16.NPF.xlsx')

KCTD17.NPF <- lm(KCTD17 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(KCTD17.NPF)
KCTD17.NPF <- tidy(KCTD17.NPF)
write.xlsx(KCTD17.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/KCTD17.NPF.xlsx')

TAF1.NPF <- lm(TAF1 ~ RIN + Ethnicity + Sex + Dev_stage, data = NPF)
summary(TAF1.NPF)
TAF1.NPF <- tidy(TAF1.NPF)
write.xlsx(TAF1.NPF,'~/Desktop/Embryonic.Expression.Work/Brain_Span/NonPrefrontal/TAF1.NPF.xlsx')

#####
#Manually import OSR.csv = Amygdala
OSR$Dev_stage <- factor(OSR$Dev_stage,
                        levels = c("Baseline",
                                   "Early Fetal",
                                   "Early Midfetal", 
                                   "Midfetal",
                                   "Late Midfetal", 
                                   "Late Fetal",
                                   "Early Infancy", 
                                   "Late Infancy",
                                   "Early Childhood", 
                                   "Late Childhood", 
                                   "Adolescence", 
                                   "Early Adulthood", 
                                   "Mid Adulthood"))
OSR$Dev_stage <- relevel(factor(OSR$Dev_stage), ref ="Baseline")
#Each gene in turn - completely failed to generate a loop!
EIF2AK2.OSR <- lm(EIF2AK2 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(EIF2AK2.OSR)
EIF2AK2.OSR <- tidy(EIF2AK2.OSR)
write.xlsx(EIF2AK2.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/EIF2AK2.OSR.xlsx')

SLC2A1.OSR <- lm(SLC2A1 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(SLC2A1.OSR)
SLC2A1.OSR <- tidy(SLC2A1.OSR)
write.xlsx(SLC2A1.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/SLC2A1.OSR.xlsx')

HPCA.OSR <- lm(HPCA ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(HPCA.OSR)
HPCA.OSR <- tidy(HPCA.OSR)
write.xlsx(HPCA.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/HPCA.OSR.xlsx')

MECR.OSR <- lm(MECR ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(MECR.OSR)
MECR.OSR <- tidy(MECR.OSR)
write.xlsx(MECR.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/MECR.OSR.xlsx')

SPR.OSR <- lm(SPR ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(SPR.OSR)
SPR.OSR <- tidy(SPR.OSR)
write.xlsx(SPR.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/SPR.OSR.xlsx')

PRKRA.OSR <- lm(PRKRA ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(PRKRA.OSR)
PRKRA.OSR <- tidy(PRKRA.OSR)
write.xlsx(PRKRA.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/PRKRA.OSR.xlsx')

PNKD.OSR <- lm(PNKD ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(PNKD.OSR)
PNKD.OSR <- tidy(PNKD.OSR)
write.xlsx(PNKD.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/PNKD.OSR.xlsx')

COL6A3.OSR <- lm(COL6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(COL6A3.OSR)
COL6A3.OSR <- tidy(COL6A3.OSR)
write.xlsx(COL6A3.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/COL6A3.OSR.xlsx')

ADCY5.OSR <- lm(ADCY5 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(ADCY5.OSR)
ADCY5.OSR <- tidy(ADCY5.OSR)
write.xlsx(ADCY5.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/ADCY5.OSR.xlsx')

SQSTM1.OSR <- lm(SQSTM1 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(SQSTM1.OSR)
SQSTM1.OSR <- tidy(SQSTM1.OSR)
write.xlsx(SQSTM1.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/SQSTM1.OSR.xlsx')

SLC6A3.OSR <- lm(SLC6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(SLC6A3.OSR)
SLC6A3.OSR <- tidy(SLC6A3.OSR)
write.xlsx(SLC6A3.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/SLC6A3.OSR.xlsx')

SGCE.OSR <- lm(SGCE ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(SGCE.OSR)
SGCE.OSR <- tidy(SGCE.OSR)
write.xlsx(SGCE.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/SGCE.OSR.xlsx')

RELN.OSR <- lm(RELN ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(RELN.OSR)
RELN.OSR <- tidy(RELN.OSR)
write.xlsx(RELN.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/RELN.OSR.xlsx')

THAP1.OSR <- lm(THAP1 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(THAP1.OSR)
THAP1.OSR <- tidy(THAP1.OSR)
write.xlsx(THAP1.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/THAP1.OSR.xlsx')

TOR1A.OSR <- lm(TOR1A ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(TOR1A.OSR)
TOR1A.OSR <- tidy(TOR1A.OSR)
write.xlsx(TOR1A.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/TOR1A.OSR.xlsx')

CIZ1.OSR <- lm(CIZ1 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(CIZ1.OSR)
CIZ1.OSR <- tidy(CIZ1.OSR)
write.xlsx(CIZ1.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/CIZ1.OSR.xlsx')

KCNMA1.OSR <- lm(KCNMA1 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(KCNMA1.OSR)
KCNMA1.OSR <- tidy(KCNMA1.OSR)
write.xlsx(KCNMA1.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/KCNMA1.OSR.xlsx')

DNAJC12.OSR <- lm(DNAJC12 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(DNAJC12.OSR)
DNAJC12.OSR <- tidy(DNAJC12.OSR)
write.xlsx(DNAJC12.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/DNAJC12.OSR.xlsx')

SLC18A2.OSR <- lm(SLC18A2 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(SLC18A2.OSR)
SLC18A2.OSR <- tidy(SLC18A2.OSR)
write.xlsx(SLC18A2.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/SLC18A2.OSR.xlsx')

ANO3.OSR <- lm(ANO3 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(ANO3.OSR)
ANO3.OSR <- tidy(ANO3.OSR)
write.xlsx(ANO3.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/ANO3.OSR.xlsx')

TH.OSR <- lm(TH ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(TH.OSR)
TH.OSR <- tidy(TH.OSR)
write.xlsx(TH.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/TH.OSR.xlsx')

NKX2.1.OSR <- lm(NKX2.1 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(NKX2.1.OSR)
NKX2.1.OSR <- tidy(NKX2.1.OSR)
write.xlsx(NKX2.1.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/NKX2.1.OSR.xlsx')

GCH1.OSR <- lm(GCH1 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(GCH1.OSR)
GCH1.OSR <- tidy(GCH1.OSR)
write.xlsx(GCH1.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/GCH1.OSR.xlsx')

PRRT2.OSR <- lm(PRRT2 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(PRRT2.OSR)
PRRT2.OSR <- tidy(PRRT2.OSR)
write.xlsx(PRRT2.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/PRRT2.OSR.xlsx')

GNAO1.OSR <- lm(GNAO1 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(GNAO1.OSR)
GNAO1.OSR <- tidy(GNAO1.OSR)
write.xlsx(GNAO1.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/GNAO1.OSR.xlsx')

GNAL.OSR <- lm(GNAL ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(GNAL.OSR)
GNAL.OSR <- tidy(GNAL.OSR)
write.xlsx(GNAL.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/GNAL.OSR.xlsx')

TUBB4A.OSR <- lm(TUBB4A ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(TUBB4A.OSR)
TUBB4A.OSR <- tidy(TUBB4A.OSR)
write.xlsx(TUBB4A.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/TUBB4A.OSR.xlsx')

MLL4.OSR <- lm(MLL4 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(MLL4.OSR)
MLL4.OSR <- tidy(MLL4.OSR)
write.xlsx(MLL4.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/MLL4.OSR.xlsx')

ATP1A3.OSR <- lm(ATP1A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(ATP1A3.OSR)
ATP1A3.OSR <- tidy(ATP1A3.OSR)
write.xlsx(ATP1A3.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/ATP1A3.OSR.xlsx')

VPS16.OSR <- lm(VPS16 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(VPS16.OSR)
VPS16.OSR <- tidy(VPS16.OSR)
write.xlsx(VPS16.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/VPS16.OSR.xlsx')

KCTD17.OSR <- lm(KCTD17 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(KCTD17.OSR)
KCTD17.OSR <- tidy(KCTD17.OSR)
write.xlsx(KCTD17.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/KCTD17.OSR.xlsx')

TAF1.OSR <- lm(TAF1 ~ RIN + Ethnicity + Sex + Dev_stage, data = OSR)
summary(TAF1.OSR)
TAF1.OSR <- tidy(TAF1.OSR)
write.xlsx(TAF1.OSR,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Amygdala/TAF1.OSR.xlsx')

#####
#Manually import Hippocampus.csv 
Hippocampus$Dev_stage <- factor(Hippocampus$Dev_stage,
                        levels = c("Baseline",
                                   "Early Fetal",
                                   "Early Midfetal", 
                                   "Midfetal",
                                   "Late Midfetal", 
                                   "Late Fetal",
                                   "Early Infancy", 
                                   "Late Infancy",
                                   "Early Childhood", 
                                   "Late Childhood", 
                                   "Adolescence", 
                                   "Early Adulthood", 
                                   "Mid Adulthood"))
Hippocampus$Dev_stage <- relevel(factor(Hippocampus$Dev_stage), ref ="Baseline")
#Each gene in turn - completely failed to generate a loop!
EIF2AK2.Hippocampus <- lm(EIF2AK2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(EIF2AK2.Hippocampus)
EIF2AK2.Hippocampus <- tidy(EIF2AK2.Hippocampus)
write.xlsx(EIF2AK2.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/EIF2AK2.Hippocampus.xlsx')

SLC2A1.Hippocampus <- lm(SLC2A1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(SLC2A1.Hippocampus)
SLC2A1.Hippocampus <- tidy(SLC2A1.Hippocampus)
write.xlsx(SLC2A1.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/SLC2A1.Hippocampus.xlsx')

HPCA.Hippocampus <- lm(HPCA ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(HPCA.Hippocampus)
HPCA.Hippocampus <- tidy(HPCA.Hippocampus)
write.xlsx(HPCA.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/HPCA.Hippocampus.xlsx')

MECR.Hippocampus <- lm(MECR ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(MECR.Hippocampus)
MECR.Hippocampus <- tidy(MECR.Hippocampus)
write.xlsx(MECR.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/MECR.Hippocampus.xlsx')

SPR.Hippocampus <- lm(SPR ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(SPR.Hippocampus)
SPR.Hippocampus <- tidy(SPR.Hippocampus)
write.xlsx(SPR.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/SPR.Hippocampus.xlsx')

PRKRA.Hippocampus <- lm(PRKRA ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(PRKRA.Hippocampus)
PRKRA.Hippocampus <- tidy(PRKRA.Hippocampus)
write.xlsx(PRKRA.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/PRKRA.Hippocampus.xlsx')

PNKD.Hippocampus <- lm(PNKD ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(PNKD.Hippocampus)
PNKD.Hippocampus <- tidy(PNKD.Hippocampus)
write.xlsx(PNKD.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/PNKD.Hippocampus.xlsx')

COL6A3.Hippocampus <- lm(COL6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(COL6A3.Hippocampus)
COL6A3.Hippocampus <- tidy(COL6A3.Hippocampus)
write.xlsx(COL6A3.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/COL6A3.Hippocampus.xlsx')

ADCY5.Hippocampus <- lm(ADCY5 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(ADCY5.Hippocampus)
ADCY5.Hippocampus <- tidy(ADCY5.Hippocampus)
write.xlsx(ADCY5.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/ADCY5.Hippocampus.xlsx')

SQSTM1.Hippocampus <- lm(SQSTM1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(SQSTM1.Hippocampus)
SQSTM1.Hippocampus <- tidy(SQSTM1.Hippocampus)
write.xlsx(SQSTM1.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/SQSTM1.Hippocampus.xlsx')

SLC6A3.Hippocampus <- lm(SLC6A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(SLC6A3.Hippocampus)
SLC6A3.Hippocampus <- tidy(SLC6A3.Hippocampus)
write.xlsx(SLC6A3.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/SLC6A3.Hippocampus.xlsx')

SGCE.Hippocampus <- lm(SGCE ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(SGCE.Hippocampus)
SGCE.Hippocampus <- tidy(SGCE.Hippocampus)
write.xlsx(SGCE.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/SGCE.Hippocampus.xlsx')

RELN.Hippocampus <- lm(RELN ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(RELN.Hippocampus)
RELN.Hippocampus <- tidy(RELN.Hippocampus)
write.xlsx(RELN.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/RELN.Hippocampus.xlsx')

THAP1.Hippocampus <- lm(THAP1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(THAP1.Hippocampus)
THAP1.Hippocampus <- tidy(THAP1.Hippocampus)
write.xlsx(THAP1.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/THAP1.Hippocampus.xlsx')

TOR1A.Hippocampus <- lm(TOR1A ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(TOR1A.Hippocampus)
TOR1A.Hippocampus <- tidy(TOR1A.Hippocampus)
write.xlsx(TOR1A.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/TOR1A.Hippocampus.xlsx')

CIZ1.Hippocampus <- lm(CIZ1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(CIZ1.Hippocampus)
CIZ1.Hippocampus <- tidy(CIZ1.Hippocampus)
write.xlsx(CIZ1.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/CIZ1.Hippocampus.xlsx')

KCNMA1.Hippocampus <- lm(KCNMA1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(KCNMA1.Hippocampus)
KCNMA1.Hippocampus <- tidy(KCNMA1.Hippocampus)
write.xlsx(KCNMA1.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/KCNMA1.Hippocampus.xlsx')

DNAJC12.Hippocampus <- lm(DNAJC12 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(DNAJC12.Hippocampus)
DNAJC12.Hippocampus <- tidy(DNAJC12.Hippocampus)
write.xlsx(DNAJC12.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/DNAJC12.Hippocampus.xlsx')

SLC18A2.Hippocampus <- lm(SLC18A2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(SLC18A2.Hippocampus)
SLC18A2.Hippocampus <- tidy(SLC18A2.Hippocampus)
write.xlsx(SLC18A2.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/SLC18A2.Hippocampus.xlsx')

ANO3.Hippocampus <- lm(ANO3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(ANO3.Hippocampus)
ANO3.Hippocampus <- tidy(ANO3.Hippocampus)
write.xlsx(ANO3.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/ANO3.Hippocampus.xlsx')

TH.Hippocampus <- lm(TH ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(TH.Hippocampus)
TH.Hippocampus <- tidy(TH.Hippocampus)
write.xlsx(TH.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/TH.Hippocampus.xlsx')

NKX2.1.Hippocampus <- lm(NKX2.1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(NKX2.1.Hippocampus)
NKX2.1.Hippocampus <- tidy(NKX2.1.Hippocampus)
write.xlsx(NKX2.1.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/NKX2.1.Hippocampus.xlsx')

GCH1.Hippocampus <- lm(GCH1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(GCH1.Hippocampus)
GCH1.Hippocampus <- tidy(GCH1.Hippocampus)
write.xlsx(GCH1.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/GCH1.Hippocampus.xlsx')

PRRT2.Hippocampus <- lm(PRRT2 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(PRRT2.Hippocampus)
PRRT2.Hippocampus <- tidy(PRRT2.Hippocampus)
write.xlsx(PRRT2.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/PRRT2.Hippocampus.xlsx')

GNAO1.Hippocampus <- lm(GNAO1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(GNAO1.Hippocampus)
GNAO1.Hippocampus <- tidy(GNAO1.Hippocampus)
write.xlsx(GNAO1.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/GNAO1.Hippocampus.xlsx')

GNAL.Hippocampus <- lm(GNAL ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(GNAL.Hippocampus)
GNAL.Hippocampus <- tidy(GNAL.Hippocampus)
write.xlsx(GNAL.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/GNAL.Hippocampus.xlsx')

TUBB4A.Hippocampus <- lm(TUBB4A ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(TUBB4A.Hippocampus)
TUBB4A.Hippocampus <- tidy(TUBB4A.Hippocampus)
write.xlsx(TUBB4A.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/TUBB4A.Hippocampus.xlsx')

MLL4.Hippocampus <- lm(MLL4 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(MLL4.Hippocampus)
MLL4.Hippocampus <- tidy(MLL4.Hippocampus)
write.xlsx(MLL4.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/MLL4.Hippocampus.xlsx')

ATP1A3.Hippocampus <- lm(ATP1A3 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(ATP1A3.Hippocampus)
ATP1A3.Hippocampus <- tidy(ATP1A3.Hippocampus)
write.xlsx(ATP1A3.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/ATP1A3.Hippocampus.xlsx')

VPS16.Hippocampus <- lm(VPS16 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(VPS16.Hippocampus)
VPS16.Hippocampus <- tidy(VPS16.Hippocampus)
write.xlsx(VPS16.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/VPS16.Hippocampus.xlsx')

KCTD17.Hippocampus <- lm(KCTD17 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(KCTD17.Hippocampus)
KCTD17.Hippocampus <- tidy(KCTD17.Hippocampus)
write.xlsx(KCTD17.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/KCTD17.Hippocampus.xlsx')

TAF1.Hippocampus <- lm(TAF1 ~ RIN + Ethnicity + Sex + Dev_stage, data = Hippocampus)
summary(TAF1.Hippocampus)
TAF1.Hippocampus <- tidy(TAF1.Hippocampus)
write.xlsx(TAF1.Hippocampus,'~/Desktop/Embryonic.Expression.Work/Brain_Span/Hippocampus/TAF1.Hippocampus.xlsx')

########## Plots #########
#Testing Heatmaps
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

#Cerebellum
set.seed(123)
cerebellum <- data.matrix(Heatmap_Cerebellum[,c(2:13)])
rownames(cerebellum) = paste0(Heatmap_Cerebellum$Genes)
column_split = rep("Prenatal", 12)
column_split[6:12] = "Postnatal"
column_split <- factor(column_split, levels =c("Prenatal", "Postnatal"))

Heatmap(cerebellum, cluster_columns = FALSE,cluster_rows = TRUE,
        show_column_dend = FALSE,show_row_dend = FALSE,
        name = "Cerebellum",
        column_split = column_split,
        column_names_gp = gpar(fontsize = 20),
        row_names_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 40),
        column_title_gp = gpar(fontsize = 40),
        row_dend_width = unit(3, "cm"), show_row_names = TRUE,
        heatmap_legend_param = list(title="T Statistic", at = seq(-10, 50, 5)))

#MotorSensory Cortex
set.seed(123)
MS <- data.matrix(Heatmap_MS[,c(2:13)])
rownames(MS) = paste0(Heatmap_MS$Genes)
column_split = rep("Prenatal", 12)
column_split[6:12] = "Postnatal"
column_split <- factor(column_split, levels =c("Prenatal", "Postnatal"))
Heatmap(MS, cluster_columns = FALSE,cluster_rows = TRUE,
        show_column_dend = FALSE,show_row_dend = FALSE,
        name = "Motor-Sensory Cortex",
        column_split = column_split,
        column_names_gp = gpar(fontsize = 20),
        row_names_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 40),
        column_title_gp = gpar(fontsize = 40),
        row_dend_width = unit(3, "cm"), show_row_names = TRUE,
        heatmap_legend_param = list(title="T Statistic", at = seq(-10, 50, 5)))

#Prefrontal
set.seed(123)
Prefrontal <- data.matrix(Heatmap_Prefrontal[,c(2:13)])
rownames(Prefrontal) = paste0(Heatmap_Prefrontal$Genes)
column_split = rep("Prenatal", 12)
column_split[6:12] = "Postnatal"
column_split <- factor(column_split, levels =c("Prenatal", "Postnatal"))
Heatmap(Prefrontal, cluster_columns = FALSE,cluster_rows = TRUE,
        show_column_dend = FALSE,show_row_dend = FALSE,
        name = "Prefrontal Cortex",
        column_split = column_split,
        column_names_gp = gpar(fontsize = 20),
        row_names_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 40),
        column_title_gp = gpar(fontsize = 40),
        row_dend_width = unit(3, "cm"), show_row_names = TRUE,
        heatmap_legend_param = list(title="T Statistic", at = seq(-10, 50, 5)))

#Striatum
set.seed(123)
Striatum <- data.matrix(Heatmap_Striatum[,c(2:12)])
rownames(Striatum) = paste0(Heatmap_Striatum$Genes)
column_split = rep("Prenatal", 11)
column_split[6:11] = "Postnatal"
column_split <- factor(column_split, levels =c("Prenatal", "Postnatal"))
Heatmap(Striatum, cluster_columns = FALSE,cluster_rows = TRUE,
        show_column_dend = FALSE,show_row_dend = FALSE,
        name = "Motor-Sensory Cortex",
        column_split = column_split,
        column_names_gp = gpar(fontsize = 20),
        row_names_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 40),
        column_title_gp = gpar(fontsize = 40),
        row_dend_width = unit(3, "cm"), show_row_names = TRUE,
        heatmap_legend_param = list(title="T Statistic", at = seq(-10, 50, 5)))
#Thalamus
set.seed(123)
Thalamus <- data.matrix(Heatmap_Thalamus[,c(2:13)])
rownames(Thalamus) = paste0(Heatmap_Thalamus$Genes)
column_split = rep("Prenatal", 12)
column_split[6:12] = "Postnatal"
column_split <- factor(column_split, levels =c("Prenatal", "Postnatal"))
Heatmap(Thalamus, cluster_columns = FALSE,cluster_rows = TRUE,
        show_column_dend = FALSE,show_row_dend = FALSE,
        name = "Motor-Sensory Cortex",
        column_split = column_split,
        column_names_gp = gpar(fontsize = 20),
        row_names_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 40),
        column_title_gp = gpar(fontsize = 40),
        row_dend_width = unit(3, "cm"), show_row_names = TRUE,
        heatmap_legend_param = list(title="T Statistic", at = seq(-10, 50, 5)))

#Non-Prefrontal Cortex
set.seed(123)
NPFC <- data.matrix(Heatmap_NPFC[,c(2:13)])
rownames(NPFC) = paste0(Heatmap_NPFC$Genes)
column_split = rep("Prenatal", 12)
column_split[6:12] = "Postnatal"
column_split <- factor(column_split, levels =c("Prenatal", "Postnatal"))
Heatmap(NPFC, cluster_columns = FALSE,cluster_rows = TRUE,
        show_column_dend = FALSE,show_row_dend = FALSE,
        name = "Motor-Sensory Cortex",
        column_split = column_split,
        column_names_gp = gpar(fontsize = 20),
        row_names_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 40),
        column_title_gp = gpar(fontsize = 40),
        row_dend_width = unit(3, "cm"), show_row_names = TRUE,
        heatmap_legend_param = list(title="T Statistic", at = seq(-10, 50, 5)))

#Amygdala
set.seed(123)
Amygdala <- data.matrix(Heatmap_Amygdala[,c(2:12)])
rownames(Amygdala) = paste0(Heatmap_Amygdala$Genes)
column_split = rep("Prenatal", 11)
column_split[6:11] = "Postnatal"
column_split <- factor(column_split, levels =c("Prenatal", "Postnatal"))
Heatmap(Amygdala, cluster_columns = FALSE,cluster_rows = TRUE,
        show_column_dend = FALSE,show_row_dend = FALSE,
        name = "Motor-Sensory Cortex",
        column_split = column_split,
        column_names_gp = gpar(fontsize = 20),
        row_names_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 40),
        column_title_gp = gpar(fontsize = 40),
        row_dend_width = unit(3, "cm"), show_row_names = TRUE,
        heatmap_legend_param = list(title="T Statistic", at = seq(-10, 50, 5)))

#Hippocampus
set.seed(123)
Hippocampus <- data.matrix(Heatmap_Hippocampus[,c(2:12)])
rownames(Hippocampus) = paste0(Heatmap_Hippocampus$Genes)
column_split = rep("Prenatal", 11)
column_split[6:11] = "Postnatal"
column_split <- factor(column_split, levels =c("Prenatal", "Postnatal"))
Heatmap(Hippocampus, cluster_columns = FALSE,cluster_rows = TRUE,
        show_column_dend = FALSE,show_row_dend = FALSE,
        name = "Motor-Sensory Cortex",
        column_split = column_split,
        column_names_gp = gpar(fontsize = 20),
        row_names_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 40),
        column_title_gp = gpar(fontsize = 40),
        row_dend_width = unit(3, "cm"), show_row_names = TRUE,
        heatmap_legend_param = list(title="T Statistic", at = seq(-10, 50, 5)))

#Individual Genes Overall
set.seed(123)
Genes_Overall <- data.matrix(Heatmap_All_Gene[,c(2:13)])
rownames(Genes_Overall) = paste0(Heatmap_All_Gene$Genes)
column_split = rep("Prenatal", 12)
column_split[6:12] = "Postnatal"
column_split <- factor(column_split, levels =c("Prenatal", "Postnatal"))
row.subsections <- c(8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,
                     8,8,8,8,8,8)
row_split = data.frame(rep(c("EIF2AK2","SLC2A1","HPCA","MECR","SPR","PRKRA","PNKD",
                             "COL6A3","ADCY5","SQSTM1","SLC6A3","SGCE","RELN","THAP1",
                             "TOR1A", "CIZ1", "KCNMA1","DNAJC12","SLC18A2","ANO3","TH",
                             "NKX2.1","GCH1","PRRT2","GNAO1","GNAL","TUBB4A","MLL4","ATP1A3",
                             "VPS16", "KCTD17","TAF1"), row.subsections))

Heatmap(Genes_Overall, cluster_columns = FALSE, cluster_rows = FALSE,
        show_column_dend = FALSE,show_row_dend = FALSE,
        name = "Overall",
        column_split = column_split,
        row_split = row_split,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        column_names_gp = gpar(fontsize = 20),
        row_names_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 40),
        column_title_gp = gpar(fontsize = 40),
        row_dend_width = unit(3, "cm"), show_row_names = TRUE,
        heatmap_legend_param = list(title="T Statistic", at = seq(-50, 50, 10)))

#Individual Regions Overall
set.seed(123)
Regions_Overall <- data.matrix(Heatmap_All_Region[,c(2:13)])
rownames(Regions_Overall) = paste0(Heatmap_All_Region$Genes)
column_split = rep("Prenatal", 12)
column_split[6:12] = "Postnatal"
column_split <- factor(column_split, levels =c("Prenatal", "Postnatal"))
row.subsections <- c(32,32,32,32,32,32,32,32)
row_split = data.frame(rep(c("Prefrontal","MS","NPFC","Striatum","Thalamus",
                             "Cerebellum","Hippocampus","Amygdala"), row.subsections))

Heatmap(Regions_Overall, cluster_columns = FALSE, cluster_rows = FALSE,
        show_column_dend = FALSE,show_row_dend = FALSE,
        name = "Regions",
        column_split = column_split,
        row_split = row_split,
        cluster_row_slices = FALSE,
        cluster_column_slices = FALSE,
        column_names_gp = gpar(fontsize = 20),
        row_names_gp = gpar(fontsize = 20),
        row_title_gp = gpar(fontsize = 40),
        column_title_gp = gpar(fontsize = 40),
        row_dend_width = unit(3, "cm"), show_row_names = TRUE,
        heatmap_legend_param = list(title="T Statistic", at = seq(-50, 50, 10)))



#BrainSeq
#extracting data
library ("SummarizedExperiment")
load("~/Desktop/rse_gene_BrainSeq_Phase1_hg19_TopHat2_EnsemblV75.rda")
rse_gene

View(colData(rse_gene))
View(rowData(rse_gene))
rowData(rse_gene)[,c("Symbol")]

#Filtering
genes_to_filter <- c("SLC2A1","HPCA","MECR","SPR","EIF2AK2","PRKRA","PNKD","COL6A3","ADCY5","SQSTM1","SLC6A3","SGCE",
                     "RELN","THAP1","TOR1A","CIZ1","KCNMA1","DNAJC12","SLC18A2","ANO3","TH","NKX2-1","GCH1","PRRT2","GNAO1","GNAL","TUBB4A","KMT2B",
                     "ATP1A3","VPS16","KCTD17","TAF1")

#subset rowData for genes of interest
interesting_genes <- rowData(rse_gene) %>% 
  as.data.frame() %>% 
  rownames_to_column("ensembl_id") %>% 
  filter("Symbol" %in% genes_to_filter)

#subset countData for genes of interest
expression_data <- rse_gene@assays$data@listData$counts %>%  as.data.frame() %>% 
  rownames_to_column("ensembl_id") %>% 
  filter(ensembl_id %in% Filtered$ensembl_id)

install.packages('openxlsx') 
library(openxlsx)
library(broom)

expression_data <- tidy(expression_data)
WriteXLS(expression_data,'~/Desktop/expression_data.xlsx')
install.packages("writexl")

#other demographic data needed for analysis
RNum<- as.data.frame(rse_gene$RNum)
BrNum <- as.data.frame(rse_gene$BrNum)
RIN<- as.data.frame(rse_gene$RIN)
Age<- as.data.frame(rse_gene$Age)
Sex<- as.data.frame(rse_gene$Sex)
Race<- as.data.frame(rse_gene$Race)
Region<- as.data.frame(rse_gene$Region)
Dx<- as.data.frame(rse_gene$Dx)

demographics <- c(RNum, BrNum, RIN, Age, Sex, Race, Region, Dx)
demographics <- as.data.frame(demographics)
demographics_2.0 <- demographics %>% filter(RIN>=7)
WriteXLS(demographics_2.0,'~/Desktop/demographics_2.0.xlsx')
demographics_3.0 <- demographics_2.0 %>% filter(demographics_2.0$rse_gene.Region=="DLPFC")
WriteXLS(demographics_3.0,'~/Desktop/demographics_3.0.xlsx')
demographics_4.0 <- demographics_3.0 %>% filter(demographics_3.0$rse_gene.Dx=="Control")
WriteXLS(demographics_4.0,'~/Desktop/demographics_4.0.1.xlsx')

#Demographics_5.0 edited by hand + added in late fetal stage to early infancy stage
#Demographics_6.0 added in BrNum to enable cross-checking with BrainSpan dataset
#Demographics_7.0 remove cases with overlap to BrainSpan
#Demographics_8.0 reference each sample to baseline of 0
#Demographics_9.0 updated as age likely incorrect in Demographics_8.0 (raw data probably pcw, needed changing to days to match with BrainSpan analysis)
demographics_9.0$Ethnicity <- factor(demographics_9.0$Ethnicity, 
                                     levels=c("CAUC", "AA"))
demographics_9.0$Sex <- factor(demographics_9.0$Sex, 
                               levels=c("M", "F"))
demographics_9.0$Stage <- factor(demographics_9.0$Stage,
                                 levels=c("Baseline", "Late Fetal", "Early Infancy", "Late Infancy", "Early Childhood"))
demographics_9.0$Stage <- relevel(factor(demographics_9.0$Stage), ref ="Baseline")

EIF2AK2 <- lm(EIF2AK2 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
EIF2AK2 <- tidy(EIF2AK2)
WriteXLS(EIF2AK2,'~/Desktop/BrainSeq/EIF2AK2.xlsx')

SLC2A1 <- lm(SLC2A1 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
SLC2A1 <- tidy(SLC2A1)
WriteXLS(SLC2A1,'~/Desktop/BrainSeq/SLC2A1.xlsx')

HPCA <- lm(HPCA ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
HPCA <- tidy(HPCA)
WriteXLS(HPCA,'~/Desktop/BrainSeq/HPCA.xlsx')

MECR <- lm(MECR ~ RIN + Ethnicity + Sex + Stage , data = demographics_9.0)
MECR <- tidy(MECR)
WriteXLS(MECR,'~/Desktop/BrainSeq/MECR.xlsx')

SPR <- lm(SPR ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
SPR <- tidy(SPR)
WriteXLS(SPR,'~/Desktop/BrainSeq/SPR.xlsx')

PRKRA <- lm(PRKRA ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
PRKRA <- tidy(PRKRA)
WriteXLS(PRKRA,'~/Desktop/BrainSeq/PRKRA.xlsx')

PNKD <- lm(PNKD ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
PNKD <- tidy(PNKD)
WriteXLS(PNKD,'~/Desktop/BrainSeq/PNKD.xlsx')

COL6A3 <- lm(COL6A3 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
COL6A3 <- tidy(COL6A3)
WriteXLS(COL6A3,'~/Desktop/BrainSeq/COL6A3.xlsx')

ADCY5 <- lm(ADCY5 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
ADCY5 <- tidy(ADCY5)
WriteXLS(ADCY5,'~/Desktop/BrainSeq/ADCY5.xlsx')

SQSTM1 <- lm(SQSTM1 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
SQSTM1 <- tidy(SQSTM1)
WriteXLS(SQSTM1,'~/Desktop/BrainSeq/SQSTM1.xlsx')

SLC6A3 <- lm(SLC6A3 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
SLC6A3 <- tidy(SLC6A3)
WriteXLS(SLC6A3,'~/Desktop/BrainSeq/SLC6A3.xlsx')

SGCE <- lm(SGCE ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
SGCE <- tidy(SGCE)
WriteXLS(SGCE,'~/Desktop/BrainSeq/SGCE.xlsx')

RELN <- lm(RELN ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
RELN <- tidy(RELN)
WriteXLS(RELN,'~/Desktop/BrainSeq/RELN.xlsx')

THAP1 <- lm(THAP1 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
THAP1 <- tidy(THAP1)
WriteXLS(THAP1,'~/Desktop/BrainSeq/THAP1.xlsx')

TOR1A <- lm(TOR1A ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
TOR1A <- tidy(TOR1A)
WriteXLS(TOR1A,'~/Desktop/BrainSeq/TOR1A.xlsx')

CIZ1 <- lm(CIZ1 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
CIZ1 <- tidy(CIZ1)
WriteXLS(CIZ1,'~/Desktop/BrainSeq/CIZ1.xlsx')

KCNMA1 <- lm(KCNMA1 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
KCNMA1 <- tidy(KCNMA1)
WriteXLS(KCNMA1,'~/Desktop/BrainSeq/KCNMA1.xlsx')

DNAJC12 <- lm(DNAJC12 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
DNAJC12 <- tidy(DNAJC12)
WriteXLS(DNAJC12,'~/Desktop/BrainSeq/DNAJC12.xlsx')

SLC18A2 <- lm(SLC18A2 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
SLC18A2 <- tidy(SLC18A2)
WriteXLS(SLC18A2,'~/Desktop/BrainSeq/SLC18A2.xlsx')

ANO3 <- lm(ANO3 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
ANO3 <- tidy(ANO3)
WriteXLS(ANO3,'~/Desktop/BrainSeq/ANO3.xlsx')

TH <- lm(TH ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
TH <- tidy(TH)
WriteXLS(TH,'~/Desktop/BrainSeq/TH.xlsx')

NKX2.1 <- lm(NKX2.1 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
NKX2.1 <- tidy(NKX2.1)
WriteXLS(NKX2.1,'~/Desktop/BrainSeq/NKX2.1.xlsx')

GCH1 <- lm(GCH1 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
GCH1 <- tidy(GCH1)
WriteXLS(GCH1,'~/Desktop/BrainSeq/GCH1.xlsx')

PRRT2 <- lm(PRRT2 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
PRRT2 <- tidy(PRRT2)
WriteXLS(PRRT2,'~/Desktop/BrainSeq/PRRT2.xlsx')

GNAO1 <- lm(GNAO1 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
GNAO1 <- tidy(GNAO1)
WriteXLS(GNAO1,'~/Desktop/BrainSeq/GNAO1.xlsx')

GNAL <- lm(GNAL ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
GNAL <- tidy(GNAL)
WriteXLS(GNAL,'~/Desktop/BrainSeq/GNAL.xlsx')

TUBB4A <- lm(TUBB4A ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
TUBB4A <- tidy(TUBB4A)
WriteXLS(TUBB4A,'~/Desktop/BrainSeq/TUBB4A.xlsx')

KMT2B <- lm(KMT2B ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
KMT2B <- tidy(KMT2B)
WriteXLS(KMT2B,'~/Desktop/BrainSeq/KMT2B.xlsx')

ATP1A3 <- lm(ATP1A3 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
ATP1A3 <- tidy(ATP1A3)
WriteXLS(ATP1A3,'~/Desktop/BrainSeq/ATP1A3.xlsx')

VPS16 <- lm(VPS16 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
VPS16 <- tidy(VPS16)
WriteXLS(VPS16,'~/Desktop/BrainSeq/VPS16.xlsx')

KCTD17 <- lm(KCTD17 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
KCTD17 <- tidy(KCTD17)
WriteXLS(KCTD17,'~/Desktop/BrainSeq/KCTD17.xlsx')

TAF1 <- lm(TAF1 ~ RIN + Ethnicity + Sex + Stage, data = demographics_9.0)
TAF1 <- tidy(TAF1)
WriteXLS(TAF1,'~/Desktop/Brainseq/TAF1.xlsx')

#Heatmap
set.seed(123)
Brainseq <- data.matrix(Brainseq_Prefrontal_Heatmap_1.0[,c(2:5)])
rownames(Brainseq) = paste0(Brainseq_Prefrontal_Heatmap_1.0$Genes)
Heatmap(Brainseq, cluster_columns = FALSE,
        name = "Dorsolateral Prefrontal Cortex",
        column_names_gp = gpar(fontsize = 7),
        row_names_gp = gpar(fontsize = 7),
        row_dend_width = unit(3, "cm"), show_row_names = TRUE,
        heatmap_legend_param = list(title="T Statistic", at = seq(-10, 50, 5)))








