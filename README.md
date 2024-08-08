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
- Failed (OOM): FCX
  
```bash
Job ID: 56157207
Cluster: hawk
State: OUT_OF_MEMORY (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 06:01:03
CPU Efficiency: 4.70% of 5-08:05:40 core-walltime
Job Wall-clock time: 06:24:17
Memory Utilized: 78.19 GB
Memory Efficiency: 80.06% of 97.66 GB
```

- Settings: `threads = 20, mem_mb = 200000`
- Failed (OOM): FCX

```bash
Job ID: 56157477
Cluster: hawk
State: OUT_OF_MEMORY (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 09:28:13
CPU Efficiency: 4.72% of 8-08:26:20 core-walltime
Job Wall-clock time: 10:01:19
Memory Utilized: 159.52 GB
Memory Efficiency: 81.68% of 195.31 GB
```

- Settings: `threads = 20, mem_mb = 200000`
- Pass: FCX
- using `future('multicore', workers = snakemake@threads)` within R

```bash

Job ID: 56159050
Cluster: hawk
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 10:04:51
CPU Efficiency: 4.94% of 8-11:59:00 core-walltime
Job Wall-clock time: 10:11:57
Memory Utilized: 182.53 GB
Memory Efficiency: 93.46% of 195.31 GB
```

- 2nd run on same same settings as above:
- This can take > 24hrs in queue for resources
- Seems to be a big discrepancy in the memory
  utilised: between 108 and 182 Gs.

```bash
Job ID: 56163010
Cluster: hawk
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 11:07:56
CPU Efficiency: 4.81% of 9-15:15:40 core-walltime
Job Wall-clock time: 11:33:47
Memory Utilized: 108.15 GB
Memory Efficiency: 55.37% of 195.31 GB
```
</details>

<details>

<summary>Benchmarking for hdWGCNA Slurm</summary>

Again problems with FCX in particular. Main problem is 
transposing the count matrix at the start of the script.

Tried close to max resources on the HTC cluster:

- Settings: `threads = 24, mem_mb = 240000, -p highmem`
- Fail: OOM

```bash
Job ID: 58052149
Cluster: hawk
User/Group: c.c1477909/c.c1477909
State: OUT_OF_MEMORY (exit code 0)
Nodes: 1
Cores per node: 24
CPU Utilized: 02:37:42
CPU Efficiency: 4.03% of 2-17:11:12 core-walltime
Job Wall-clock time: 02:42:58
Memory Utilized: 206.46 GB
Memory Efficiency: 88.09% of 234.38 GB
```

Tried to run it on the `highmem` partition, but even
approaching the absolute resource limits for that:

- Settings: `threads = 20, mem_mb = 300000, -p highmem`
- Fail: OOM

```bash
Job ID: 58057393
Cluster: hawk
User/Group: c.c1477909/c.c1477909
State: OUT_OF_MEMORY (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 02:32:34
CPU Efficiency: 4.87% of 2-04:13:20 core-walltime
Job Wall-clock time: 02:36:40
Memory Utilized: 236.51 GB
```

- Settings: `threads = 20, mem_mb = 350000, -p highmem`
- Fail: (exit code 1)

```bash
Job ID: 58100400
Cluster: hawk
User/Group: c.c1477909/c.c1477909
State: FAILED (exit code 1)
Nodes: 1
Cores per node: 20
CPU Utilized: 02:55:44
CPU Efficiency: 4.89% of 2-11:57:20 core-walltime
Job Wall-clock time: 02:59:52
Memory Utilized: 310.64 GB
Memory Efficiency: 90.89% of 341.80 GB
```

In this run R threw the following error:

```r
Error in mcfork() : 
  unable to fork, possible reason: Cannot allocate memory
```

This was due to the threads allocation in the R script (40) 
not matching that set by snakemake (20). The run passed the
point it usually fails at, so next run with same reosurces 
may complete.

***


