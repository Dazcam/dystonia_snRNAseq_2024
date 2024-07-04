rule download_public_data:
    input:  "../resources/sheets/Stiletti_downloads_table.xlsx"
    output: "../results/01R_objects/01seurat_{region}.rds" 
    singularity: "../resources/containers/seurat5b_latest.sif"
    params: root_dir = "../", 
            region = lambda wc: wc.get("region")
    resources: threads = 4, mem_mb = 20000 
    log:    "../results/00LOG/01download_public_data_{region}.log"
    script:
            "../scripts/dystonia_download_public_data.R"

rule prepare_data_seurat:
    input:  "../results/01R_objects/01seurat_{region}.rds"
    output: "../results/01R_objects/02seurat_{region}.rds"
    singularity: "../resources/containers/seurat5b_latest.sif"
    params: root_dir = "../",
            region = lambda wc: wc.get("region")
    resources: threads = 20, mem_mb = 100000, time="0-12:00:00"
    log:    "../results/00LOG/02prepare_data_seurat_{region}.log"
    script:
            "../scripts/dystonia_prepare_data_seurat.R"

rule annotate_clusters:
    input: "../results/01R_objects/02seurat_{region}.rds"
    output: "../results/01R_objects/03seurat_{region}.rds"
    singularity: "../resources/containers/seurat5b_latest.sif"
    params: root_dir = "../",
            region = lambda wc: wc.get("region")
    resources: threads = 5, mem_mb = 30000, time="0-12:00:00"
    log:    "../results/00LOG/03annotate_clusters_{region}.log"
    script:
            "../scripts/dystonia_annotate_clusters.R"

rule pseudobulk:
    input: "../results/01R_objects/03seurat_{region}.rds"
    output: "../results/01R_objects/aggr_exp_{region}.rds"
    singularity: "../resources/containers/seurat5b_latest.sif"
    params: root_dir = "../",
            region = lambda wc: wc.get("region")
    resources: threads = 5, mem_mb = 50000, time="0-12:00:00"
    log:    "../results/00LOG/04pseudobulk_{region}.log"
    script:
            "../scripts/dystonia_pseudobulk.R"

rule wgcna_fetal:
    input: "../resources/public_data/cameron_2023/seurat_{wgcna_region}.rds"
    output: "../results/03wgcna/dystonia_wgcna_{wgcna_region}.html"
    singularity: "../resources/containers/seurat5c_latest.sif"
    params: root_dir = "../",
            region = lambda wc: wc.get("wgcna_region")
    resources: threads = 4, mem_mb = 40000, time="0-12:00:00"
    log:    "../results/00LOG/05wgcna_{wgcna_region}.log"
    script:
            "../scripts/dystonia_wgcna.R"

#rule wgcna_adult:
#    input: "../results/01R_objects/03seurat_{region}.rds"
#    output: "../results/03wgcna/dystonia_wgcna_{wgcna_region}.html"
#    singularity: "../resources/containers/seurat5c_latest.sif"
#    params: root_dir = "../",
#            region = lambda wc: wc.get("wgcna_region")
#    resources: threads = 20, mem_mb = 40000, time="0-12:00:00"
#    log:    "../results/00LOG/05wgcna_{wgcna_region}.log"
#    script:
#            "../scripts/dystonia_pseudobulk.R"
