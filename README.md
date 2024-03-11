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
    R --no-echo -e 'install.packages(c("scCustomize", "readxl", "harmony")'

    apt clean
```
</details>

<details>

<summary>Benchmarking on Slurm</summary>

- Settings: `threads = 10, mem_mb = 40000`
- Pass: Str
- Failed (OOM): Cer, FCX

```bash
# FCX: Fail

Job ID: 56154891
Cluster: hawk
User/Group: c.c1477909/c.c1477909
State: OUT_OF_MEMORY (exit code 0)
Nodes: 1
Cores per node: 10
CPU Utilized: 00:03:50
CPU Efficiency: 9.06% of 00:42:20 core-walltime
Job Wall-clock time: 00:04:14
Memory Utilized: 0.00 MB (estimated maximum)
Memory Efficiency: 0.00% of 39.06 GB (39.06 GB/node)

# Str: Pass

Job ID: 56154892
Cluster: hawk
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 10
CPU Utilized: 01:03:37
CPU Efficiency: 9.87% of 10:44:20 core-walltime
Job Wall-clock time: 01:04:26
Memory Utilized: 28.52 GB
Memory Efficiency: 73.00% of 39.06 GB

# Cer: Fail

Job ID: 56154893
Cluster: hawk
State: OUT_OF_MEMORY (exit code 0)
Nodes: 1
Cores per node: 10
CPU Utilized: 01:40:45
CPU Efficiency: 9.85% of 17:02:20 core-walltime
Job Wall-clock time: 01:42:14
Memory Utilized: 27.76 GB
Memory Efficiency: 71.07% of 39.06 GB

```

- Settings: `threads = 20, mem_mb = 80000`
- Pass: Cer
- Failed (OOM): FCX

```bash
# FCX: Fail
Job ID: 56155003
Cluster: hawk
State: OUT_OF_MEMORY (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 06:44:21
CPU Efficiency: 4.46% of 6-07:05:00 core-walltime
Job Wall-clock time: 07:33:15
Memory Utilized: 73.40 GB
Memory Efficiency: 93.96% of 78.12 GB

# Cer: Pass
 
Cluster: hawk
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 02:30:13
CPU Efficiency: 5.15% of 2-00:38:40 core-walltime
Job Wall-clock time: 02:25:56
Memory Utilized: 34.36 GB
Memory Efficiency: 43.98% of 78.12 GB
```

- Settings: `threads = 20, mem_mb = 100000`


</details>


***
