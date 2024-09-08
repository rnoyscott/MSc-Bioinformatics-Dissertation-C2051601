# Import required modules
import celltypist # for predicting cell types
from celltypist import models # CellTypist pre-built models
import scanpy as sc # for handling single cell RNA-seq data
import matplotlib.pyplot as plt # for plotting
import seaborn as sns # for plotting
import pandas as pd # for data manipulation
import numpy as np # for data manipulation
from scipy import stats # for statistics
import re # for regular expressions

# Ensure `importlib-metadata` is imported correctly
try:
    import importlib_metadata
except ImportError:
    import importlib.metadata as importlib_metadata

# Functions
## UMAP plotting function
def plot_umap(color_by, title, layer, legend_mod):
    fig, ax = plt.subplots(figsize=(8, 8)) # define size
    plt.subplots_adjust(wspace=0.25) # add white space
    ## Plot UMAP
    ax = sc.pl.umap(adata, color = [color_by], layer = layer, show = False, legend_fontsize = 'small', add_outline = True, s = 10, title = title)
    ## Customise the axes
    _ = ax.set_xlabel("UMAP 1", fontweight='bold')
    _ = ax.set_ylabel("UMAP 2", fontweight='bold')
    _ = ax.set_title(title, fontweight='bold', fontsize='small')
    ## Customise the legend
    if legend_mod:
        _ = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='x-small')
    plt.tight_layout()

## Predicting cell types function
def predict_cell_types(adata, model):
    sc.pp.filter_genes(adata, min_cells=3) # genes expressd in a minimum of 3 cells
    sc.pp.normalize_total(adata, target_sum=1e4) # normalise to 10,000 counts
    sc.pp.log1p(adata) # log1p transformation
    
    adata.X = adata.X.toarray() # convert to array
    predictions = celltypist.annotate(adata, model=model, majority_voting=True) # use model to transfer labels to adata object
    predictions_adata = predictions.to_adata()
    
    ## Store the results
    result_df = predictions_adata.obs[['predicted_labels', 'conf_score', 'majority_voting']]
    result_df.columns = ['label', 'score', 'voting']
    
    return result_df

# Read in the base data and check characteristics
adata_base = sc.read_h5ad('../Input/Primary/BloodCells.h5ad')
print(adata_base)
print(adata_base.var.head())

# Extract raw counts matrix from base data
raw_counts = adata_base.raw.X
print("Number of rows in raw counts:", raw_counts.shape[0]) # print number of cells
print("Number of columns in raw counts:", raw_counts.shape[1]) # print number of genes

# Extract the gene names from var
gene_names = adata_base.raw.var_names
print(len(gene_names))
    
# Extract the required metadata columns that include patient ID and condition group ID
metadata_columns = ['orig.ident', 'AllSamples', 'ID', 'cell_ID']
metadata = adata_base.obs[metadata_columns]

# Create a fresh object with raw counts and metadata
adata = sc.AnnData(X=raw_counts, obs=metadata, var=pd.DataFrame(index=gene_names))
print(adata) # check the characteristics of the new object are similar to those of the base data

adata_anno = adata.copy() # create a copy to be used in predicting cell types

# Identify all mitochondrial genes and compute percentage MT per cell
mito_genes = adata.var_names.str.startswith('MT-')

## For each cell compute fraction of counts in mito genes vs. all genes
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
## Add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

# Filter out genes and cells
sc.pp.filter_cells(adata, min_genes=500) # remove cells with < 500 genes
sc.pp.filter_genes(adata, min_cells=3) # remove genes expressed in < 3 cells

# Remove cells with MT content > 15%
adata = adata[adata.obs.percent_mito < 0.15, :]

# Store current status of counts matrix in a seperate layer
adata.layers['filtered_counts'] = adata.X.copy()

# Normalise and transform counts
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# Store current status of counts matrix in a seperate layer
adata.layers['filtered_counts_norm'] = adata.X.copy()

# Dimensionality reduction and clustering
sc.pp.highly_variable_genes(adata, n_top_genes = 3000) # identify top 3000 HVGs

adata_hvg = adata[:, adata.var.highly_variable].copy() # create a subset datasets with just the HVGs

sc.pp.regress_out(adata_hvg, ['n_counts', 'percent_mito']) # regress out the unwanted effects of MT and gene count

sc.pp.scale(adata_hvg, max_value=10) # scale data

sc.pp.pca(adata_hvg, svd_solver='randomized') # compute PCA

sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=30) # construct kNN neighbours graph

sc.tl.umap(adata_hvg) # compute UMAP

sc.tl.leiden(adata_hvg, resolution = 3, key_added = 'overcluster') # compute Leiden clustering with resolution = 3

# Map the results of dimensionality reduction and clustering to entire dataset
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']

adata.obsm['X_umap'] = adata_hvg.obsm['X_umap']

adata.obs['overcluster'] = adata_hvg.obs['overcluster']

# Predict cell types using CellTypist
## Create a custom model using Seurat Data
### PBMC3K
adata_pbmc3k = sc.read_h5ad('../Resources/CellTypist/pbmc3k.h5ad') # read in PBMC3K
annotations = pd.read_csv("../Resources/CellTypist/pbmc3k_anno.csv") # read in the associated cell annotations

annotations.set_index("Cell_ID", inplace=True) # ensure the cell IDs are the index

annotations = annotations.reindex(adata_pbmc3k.obs_names) # re-index annotations to match the order of cells present in PBMC3K

adata_pbmc3k.obs['seurat_annotations'] = annotations['seurat_annotations'].values # map the annotations to PBMC3K

adata_pbmc3k = adata_pbmc3k[~adata_pbmc3k.obs['seurat_annotations'].isna()] # remove NA values

### IFNB
adata_ifnb = sc.read_h5ad('../Resources/CellTypist/ifnb.h5ad')
annotations = pd.read_csv("../Resources/CellTypist/ifnb_anno.csv")

annotations.set_index("Cell_ID", inplace=True)

annotations = annotations.reindex(adata_ifnb.obs_names)

adata_ifnb.obs['seurat_annotations'] = annotations['seurat_annotations'].values

adata_ifnb = adata_ifnb[~adata_ifnb.obs['seurat_annotations'].isna()]

### Merge the two datasets
adata_list = [adata_pbmc3k, adata_ifnb]
adata_pbmc3k_ifnb = sc.concat(adata_list)

### Normalise to 10,000 and log1p transform
sc.pp.normalize_total(adata_pbmc3k_ifnb, target_sum = 1e4)
sc.pp.log1p(adata_pbmc3k_ifnb)

### Train the model using CellTypist
ref_model_pbmc3k = celltypist.train(adata_pbmc3k_ifnb, labels = 'seurat_annotations', use_SGD = False, feature_selection = True, top_genes = 300)

## Obtain CellTypist 'Immune_All_Low.pk1'
models.get_all_models() # retrieve all pre-built models

ref_model = models.Model.load(model = 'Immune_All_Low.pkl') # extract 'Immune_All_Low.pk1'

## Obtain predictions from both models
predictions_ref = predict_cell_types(adata_anno.copy(), ref_model)
predictions_ref_pbmc3k = predict_cell_types(adata_anno.copy(), ref_model_pbmc3k)

## Rename columns to avoid conflicts when merging
predictions_ref.columns = ['ref_label', 'ref_score', 'majority_voting']
predictions_ref_pbmc3k.columns = ['ref_label_pbmc3k', 'ref_score_pbmc3k', 'majority_voting_pbmc3k']

## Merge the two prediction datasets on their indices
merged_predictions = predictions_ref.merge(predictions_ref_pbmc3k, left_index=True, right_index=True)

## Inspect merged predictions
print(merged_predictions.head())

## Merge predictions with the single cell RNA-seq object
adata.obs = adata.obs.merge(right = merged_predictions, left_index = True, right_index = True)

# Compute markers for each predicted immune cell type, using only the highly variable genes
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.tl.rank_genes_groups(adata_hvg, 'majority_voting', method='wilcoxon', tie_correct = True)

result = adata_hvg.uns['rank_genes_groups']
groups = result['names'].dtype.names
df = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(5)

# Extract the columns that end with '_n'
n_columns = [col for col in df.columns if col.endswith('_n')]

# Extract the first gene from each '_n' column
marker_genes = [df[col].iloc[0] for col in n_columns]

sc.pl.dotplot(adata, marker_genes, layer = 'filtered_counts', groupby='majority_voting')
plt.tight_layout()

filename_dotplot = 'Predicted_CellType_Markers_Expression_Heatmap_300dpi.png'
plt.savefig(filename_dotplot, dpi=300, format='png', bbox_inches='tight')
plt.close()

# Plot UMAP of predicted cell types using 'Immune_All_Low.pk1'
plot_umap('majority_voting', 'Majority Voting (CellTypist) Immune Cell Types', 'filtered_counts', legend_mod = True)

filename_umap = f"CellTypist_MajorityVoting_CellType_UMAP_Overlay_300dpi.png"
plt.savefig(filename_umap, dpi=300, format='png')
plt.close()

# Plot UMAP of predicted cell types using the custom model
plot_umap('majority_voting_pbmc3k', 'Majority Voting (PBMC3K) Immune Cell Types', 'filtered_counts', legend_mod = True)

filename_umap = f"CellTypist_MajorityVoting(PBMC3K)_CellType_UMAP_Overlay_300dpi.png"
plt.savefig(filename_umap, dpi=300, format='png')
plt.close()

# Plot UMAP of prediction confidence scores using 'Immune_All_Low.pk1'
plot_umap('ref_score', 'Prediction Confidence Score (CellTypist High Res Immune Cells)', 'filtered_counts', legend_mod = False)

filename_umap = f"CellTypist_RefScore_UMAP_Overlay_300dpi.png"
plt.savefig(filename_umap, dpi=300, format='png')
plt.close()

# Plot UMAP of prediction confidence scores using custom model
plot_umap('ref_score_pbmc3k', 'Prediction Confidence Score (PBMC3K)', 'filtered_counts', legend_mod = False)

filename_umap = f"CellTypist_RefScore_pbmc3k_UMAP_Overlay_300dpi.png"
plt.savefig(filename_umap, dpi=300, format='png')
plt.close()

# Plot UMAP of percentage mitochondrial genes
plot_umap('percent_mito', 'Percentage Mitochondrial Genes', 'counts', legend_mod = False)

filename_umap = f"Percentage_Mito_UMAP_Overlay_300dpi.png"
plt.savefig(filename_umap, dpi=300, format='png')
plt.close()

# Plot UMAP of gene counts
plot_umap('n_genes', 'Number of Genes', 'counts', legend_mod = False)

filename_umap = f"Gene_Count_UMAP_Overlay_300dpi.png"
plt.savefig(filename_umap, dpi=300, format='png')
plt.close()

print(adata)

# Export pre-processed data
adata.write_h5ad('../Input/Secondary/BloodCells.h5ad')