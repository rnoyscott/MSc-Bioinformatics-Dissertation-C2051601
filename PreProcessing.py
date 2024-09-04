import celltypist
from celltypist import models
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
import re

# Ensure `importlib-metadata` is imported correctly
try:
    import importlib_metadata
except ImportError:
    import importlib.metadata as importlib_metadata

adata_base = sc.read_h5ad('../../input/Data Sets_py/BloodCells.h5ad')

raw_counts = adata_base.raw.X

# Extract the gene names from var
gene_names = adata_base.raw.var_names
    
# Extract the required metadata columns
metadata_columns = ['orig.ident', 'AllSamples', 'ID', 'cell_ID']
metadata = adata_base.obs[metadata_columns]

# Create a new AnnData object with raw counts and metadata
adata = sc.AnnData(X=raw_counts, obs=metadata, var=pd.DataFrame(index=gene_names))

adata_anno = adata.copy()

# Store the raw counts in another layer
adata.layers['counts'] = adata.X.copy()

sc.pp.filter_cells(adata, min_genes=500) # 500
sc.pp.filter_genes(adata, min_cells=3)

mito_genes = adata.var_names.str.startswith('MT-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1

adata = adata[adata.obs.percent_mito < 0.15, :]

sc.pp.normalize_total(adata)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['n_counts', 'percent_mito'])

sc.pp.scale(adata, max_value=10)

sc.pp.pca(adata, svd_solver='randomized')

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

sc.tl.umap(adata)

sc.tl.leiden(adata, resolution = 0.25)

models.get_all_models()

ref_model = models.Model.load(model = 'Immune_All_Low.pkl')

def predict_cell_types (adata):
    sc.pp.filter_genes(adata, min_cells = 3)
    sc.pp.normalize_total(adata, target_sum = 1e4)
    sc.pp.log1p(adata)
    
    adata.X = adata.X.toarray()
    predictions = celltypist.annotate(adata, model = ref_model, majority_voting = True)
    predictions_adata = predictions.to_adata()
    adata.obs['ref_label'] = predictions_adata.obs.loc[adata.obs.index, 'predicted_labels']
    adata.obs['ref_score'] = predictions_adata.obs.loc[adata.obs.index, 'conf_score']
    return adata.obs

predictions = predict_cell_types(adata_anno.copy())

predictions_less = predictions[['ref_label', 'ref_score', 'majority_voting', 'over_clustering']]

adata.obs = adata.obs.merge(right = predictions_less, left_index = True, right_index = True)

sc.pl.umap(adata, color = ['majority_voting'], add_outline = True, s = 10, title = 'Majority Voting Immune Cell Types', legend_fontsize = 'small')
plt.tight_layout()

filename_umap = f"CellTypist_MajorityVoting_CellType_UMAP_Overlay_300dpi.png"
plt.savefig(filename_umap, dpi=300, format='png')
plt.close()

adata.obs['predicted_celltype_major'] = adata.obs.groupby('leiden')['ref_label'].transform(lambda x: x.mode()[0])

sc.pl.umap(adata, color = ['predicted_celltype_major'], add_outline = True, s = 10, title = 'Primary Predicted Immune Cell Types')
plt.tight_layout()

filename_umap = f"CellTypist_PrimaryCellType_UMAP_Overlay_300dpi.png"
plt.savefig(filename_umap, dpi=300, format='png')
plt.close()

sc.pl.umap(adata, color = ['ref_score'], s = 10, add_outline = True, title = 'Prediction Confidence Score')
plt.tight_layout()

filename_umap = f"CellTypist_RefScore_UMAP_Overlay_300dpi.png"
plt.savefig(filename_umap, dpi=300, format='png')
plt.close()

sc.tl.rank_genes_groups(adata, 'majority_voting', method='wilcoxon', tie_correct = True)

result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
df = pd.DataFrame(
    {group + '_' + key[:1]: result[key][group]
    for group in groups for key in ['names', 'pvals']}).head(5)

# Extract the columns that end with '_n'
n_columns = [col for col in df.columns if col.endswith('_n')]

# Extract the first gene from each '_n' column
marker_genes = [df[col].iloc[0] for col in n_columns]

sc.pl.dotplot(adata, marker_genes, groupby='majority_voting')
plt.tight_layout()

filename_dotplot = 'Predicted_CellType_Markers_Expression_Heatmap_300dpi.png'
plt.savefig(filename_dotplot, dpi=300, format='png', bbox_inches='tight')
plt.close()

sc.tl.leiden(adata, resolution = 3)

adata.obs['predicted_celltype_minor'] = adata.obs.groupby('leiden')['ref_label'].transform(lambda x: x.mode()[0])

print(adata.obs['predicted_celltype_major'].unique())
print(adata.obs['predicted_celltype_minor'].unique())

adata.write_h5ad('BloodCellsnew.h5ad')

print(adata)

def create_filename(name):
    # Replace invalid characters with underscores
    return re.sub(r'[\/\\:*?"<>|]', '_', name)

# Manually check clusters
celltype_markers = pd.read_csv('cluster_predictions_markers.csv')
celltype_dict = celltype_markers.set_index('ref_label')['marker_genes'].to_dict()

for celltype, marker_genes in celltype_dict.items():
    celltype_filename = create_filename(celltype)

    # Convert marker_genes string to a list of gene names
    marker_genes_list = [gene.strip() for gene in marker_genes.split(',')]

    # Check if marker_genes are present in adata
    present_marker_genes = [gene for gene in marker_genes_list if gene in adata.var_names]
    print(f"Marker genes present in adata for {celltype}:", present_marker_genes)

    if present_marker_genes:
        # Plot UMAP
        sc.pl.umap(adata, color=present_marker_genes, s=10, add_outline=True)
        filename_umap = f'{celltype_filename}_Markers_UMAP_300dpi.png'
        plt.savefig(filename_umap, dpi=300, format='png', bbox_inches='tight')
        plt.close()


        # Plot Violin plots
        sc.pl.violin(adata, present_marker_genes, groupby='majority_voting', rotation=90)
        filename_violin = f'{celltype_filename}_Markers_Violin_300dpi.png'
        plt.savefig(filename_violin, dpi=300, format='png', bbox_inches='tight')
        plt.close()
    else:
        print(f"None of the marker genes for {celltype} are present in adata.")

labels = adata.obs[['ref_label', 'over_clustering']].groupby('over_clustering').agg(lambda x: x.mode()[0] if not x.mode().empty else pd.NA)
scores = adata.obs[['ref_score', 'over_clustering']].groupby('over_clustering').agg(lambda x: x.mean())
cluster_mapping = labels.merge(right = scores, left_index = True, right_index = True)

print(cluster_mapping)

cluster_mapping.to_csv('cluster_predictions.csv')

#sc.tl.rank_genes_groups(adata, groupby = 'over_clustering')
#marks = sc.get.rank_genes_groups_df(adata, group = None)

# Label primary cell types

# Copy the 'over_clustering' column to a new column 'over_clustering_copy'
#adata.obs['over_clustering_copy'] = adata.obs['over_clustering'].copy()

# Map the new column to the dictionary
#adata.obs['predicted_celltype_primary'] = adata.obs['over_clustering_copy'].map(dictionary_with_verified_celltypes)

# Verify sub-types and map to over_clustering
#adata.obs['predicted_celltype_secondary'] = adata.obs['over_clustering_copy'].map(dictionary_with_verified_celltypes)
