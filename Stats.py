from scipy import stats
import scikit_posthocs as sp
from scipy.stats import kruskal
import scanpy as sc
import pandas as pd
import numpy as np

# Jaccard index of CD4+ T cells
adata = sc.read_h5ad('../Input/Secondary/BloodCells.h5ad')
adata_base = sc.read_h5ad('../Input/Primary/BloodCells.h5ad')

annotations = pd.read_csv("../Resources/CellTypist/ifnb_anno.csv")
print(f"Unique labels: {annotations['seurat_annotations'].unique()}")

annotations = pd.read_csv("../Resources/CellTypist/pbmc3k_anno.csv")
print(f"Unique labels: {annotations['seurat_annotations'].unique()}")

print(f'Dimensions of base data: {adata_base.shape}')

CD4Tcells = ['Tcm/Naive helper T cells', 'Tem/Effector helper T cells']

Subset = adata[adata.obs['majority_voting'].isin(CD4Tcells)].copy()

cell_ids = Subset.obs.index.tolist()

CD4Tcells_base = ['Tcm/Naive helper T cells', 'Tem/Effector helper T cells']

CD4Tcells_base = ["CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 Proliferating"]

Subset_base = adata_base[adata_base.obs['predicted.celltype.l2'].isin(CD4Tcells_base)].copy()

cell_ids_base = Subset_base.obs.index.tolist()

set1 = set(cell_ids)
set2 = set(cell_ids_base)

# Calculate the intersection and union of the sets
intersection = set1.intersection(set2)
union = set1.union(set2)

# Calculate the Jaccard index
jaccard_index = len(intersection) / len(union)

# Display the Jaccard index
print(jaccard_index)

# Statistical difference between groups for violin plots

gene = 'KAT2B'

adata_CD4 = sc.read_h5ad('../Input/Secondary/CD4Tcells.h5ad') # read in the subset data file
gene_expression = adata_CD4[:, adata_CD4.var.index == gene].layers['filtered_counts'].toarray().flatten() # extract the gene expression data from the filtered countd layer

AllSamples = adata_CD4.obs['AllSamples']

# Combine expression data with group information into a DataFrame
data = pd.DataFrame({
    'expression': gene_expression,
    'group': AllSamples
})

# Perform Kruskal-Wallis test to determine if there are differences among groups
kruskal_result = kruskal(*[data[data['group'] == g]['expression'].values for g in data['group'].unique()])

print("Kruskal-Wallis test statistic:", kruskal_result.statistic)
print("Kruskal-Wallis p-value:", kruskal_result.pvalue)

# Perform Dunn test
dunn_result = sp.posthoc_dunn(data, val_col='expression', group_col='group', p_adjust='bonferroni')

# Display Dunn test results
print(dunn_result)