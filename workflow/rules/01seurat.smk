rule download_public_data:
    input:  "../resources/sheets/Stiletti_downloads_table.xlsx"
    output: "../results/01R_objects/prelim/seurat_{region}.rds" 
    singularity: "../resources/containers/seurat5b_latest.sif"
    params: root_dir = "../", 
            region = lambda wc: wc.get("region")
    resources: threads = 4, mem_mb = 20000 
    log:    "../results/00LOG/01download_public_data_{region}.log"
    script:
            "../scripts/dystonia_download_public_data.R"

rule prepare_data_seurat:
    input:  "../results/01R_objects/prelim/seurat_{region}.rds"
    output: "../results/01R_objects/seurat_{region}.rds"
    singularity: "../resources/containers/seurat5b_latest.sif"
    params: root_dir = "../",
            region = lambda wc: wc.get("region")
    resources: threads = 10, mem_mb = 40000
    log:    "../results/00LOG/02prepare_data_seurat_{region}.log"
    script:
            "../scripts/dystonia_prepare_data_seurat.R"
