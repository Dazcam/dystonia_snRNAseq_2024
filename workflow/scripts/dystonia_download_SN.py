#--------------------------------------------------------------------------------------
#
#    Download Stiletti SN h5ad file with progress bar
#
#--------------------------------------------------------------------------------------

# Substantia Nigra was requested in review. But in the intervening time, direct download
# of .rds was withdrawn from cellxgene. The R API for this dataset was corrupted. Rather
# than rejig all the previous code extract matrix and meta from h5ad and plugged these
# into Seurat. See get_dissection_data() in dystonia_functions.R

import scanpy as sc
from scipy import io
import requests
import os

# Configuration
dataset_uuid = "fd4b7829-283f-4906-8dc4-7846731894e9"
h5ad_url = f"https://datasets.cellxgene.cziscience.com/{dataset_uuid}.h5ad"
out_dir = "/Users/darren/Desktop/test/"  # Adjust to your path
out_file = os.path.join(out_dir, "SN.h5ad")

# Create output dir if it doesn't exist
os.makedirs(out_dir, exist_ok=True)

# Download with progress
print(f"Downloading Stilletti SN h5ad data from: {h5ad_url}...")
response = requests.get(h5ad_url, stream=True)
response.raise_for_status()  # Raise error on bad status

total_size = int(response.headers.get("content-length", 0))
with open(out_file, "wb") as file, requests.get(h5ad_url, stream=True) as r:  # Reuse for progress
    downloaded = 0
    for chunk in r.iter_content(chunk_size=8192):
        if chunk:
            file.write(chunk)
            downloaded += len(chunk)
            if total_size > 0:
                progress = (downloaded / total_size) * 100
                print(f"\rProgress: {progress:.1f}% ({downloaded / (1024*1024):.1f} MB)", end="", flush=True)
print(f"\nDownloaded successfully to: {out_file}")

# Load the h5ad file
print("Loading SN h5ad file into AnnData object...")
adata = sc.read(out_dir + 'SN.h5ad')

print('SN data columns:')
adata.obs.columns

print('SN counts:')
adata.X[10:20,10:20].todense()

print('SN gene meta:')
adata.var

# Export counts as mtx (note: adata.X is typically sparse; io.mmwrite handles it)
print('Saving counts matrix to SN_counts.mtx ...')
io.mmwrite('SN_counts.mtx', adata.X)

# Save metadata files (your code, adjusted for paths)
print('Saving metadata to SN_cell_meta.csv')
cell_meta = adata.obs.copy()
cell_meta['Barcode'] = cell_meta.index
cell_meta.to_csv(out_dir + 'SN_cell_meta.csv', index=None)

print('Saving metadata to SN_gene_meta.csv')
gene_meta = adata.var.copy()
gene_meta['ensembl_id'] = gene_meta.index
gene_meta.to_csv(out_dir + 'SN_gene_meta.csv', index=None)

print('All done!')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------