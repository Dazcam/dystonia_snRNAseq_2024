## Dystonia project - Kathryn Peall

## Regions of interest
  + Striatum (Merge Caudate and Putamen)
  + Cerebellum
  + Frontal Cortex

## Method

1. Download Stiletti scRNAseq data for adult brain regions
2. QC data and generate clusters using Seurat
3. Label adult cell types
4. Generate GeX violin plots for adult and fetal brain (latter from Biol Psych paper)
5. Aggregate expression for all distinct cell types (adult and fetal) and choose
a visualisation method to easily identify patterns of expression across regions 

***

## Container for Hawk

A container for all analyses was generated using the following definition file.

<details>

<summary>Definition file</summary>

```
Bootstrap: docker
From: bioconductor/bioconductor_docker:devel

%labels
    Version v0.0.1

%help
    This is a container to run Seurat 5 on Hawk

%post
    # Update the image
    apt update
    apt upgrade -y

    # for igraph
    apt install -y glpk-utils libglpk-dev

    # for sctransform
    apt install -y libicu-dev

    # for BPCells
    apt install -y libhdf5-dev

    # Install R packages
    R --no-echo -e 'remotes::install_github("bnprks/BPCells")'
    R --no-echo -e 'BiocManager::install("glmGamPoi")'
    R --no-echo -e 'install.packages("RPresto")'
    R --no-echo -e 'install.packages("Seurat")'
    R --no-echo -e 'setRepositories(ind=1:3)' # Needed to automatically install Bioconductor dependencies for Signac
    R --no-echo -e 'install.packages(c("R.utils", "Signac"))'
    R --no-echo -e 'remotes::install_github("satijalab/seurat-wrappers")'
    R --no-echo -e 'remotes::install_github("satijalab/azimuth")'
    R --no-echo -e 'BiocManager::install(c("scuttle", "scater"))'
    R --no-echo -e 'install.packages("scCustomize")'
    R --no-echo -e 'install.packages("readxl")'

    apt clean
```
</details>

***
