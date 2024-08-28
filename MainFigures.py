import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import gseapy as gp
import os
import glob
from gseapy.scipalette import SciPalette
from scipy import stats
import re

# Create a Figures directory to store outputs
dir = '../Output/Figures/'
os.makedirs(dir, exist_ok=True)

adata = sc.read_h5ad('../Input/Secondary/BloodCells.h5ad')
DEGfile = pd.read_csv('../Output/DGE/Scanpy/CD4Tcells/VE_DEG_CD4Tcells.csv')

celltypes = adata.obs['majority_voting'].unique()
total = len(adata.obs['majority_voting'])
print(total)

color_discrete_map= {'Viraemic': '#c38961', 'Control': 'lightgrey', 'Exposed': '#007d82'} # create a colour mapper for each patient group

# Plot clusters UMAP
fig, ax = plt.subplots(figsize=(10, 10))
sc.pl.umap(adata, color = 'overcluster', legend_loc = 'on data', legend_fontoutline=2, add_outline=True, s=10, ax=ax)
ax.set_xlabel("UMAP 1", fontweight='bold')
ax.set_ylabel("UMAP 2", fontweight='bold')
ax.set_title("Cluster Identities", fontsize=18, fontweight='bold')
# Save UMAP overlay as file
filename_umap = os.path.join(dir, f"newOverCluster_UMAP_Overlay_300dpi.png")
plt.savefig(filename_umap, dpi=300, format='png')

# Initialize an empty list to store dictionaries of proportions
proportions_list = []

# Get unique values from 'AllSamples'
unique_samples = adata.obs['AllSamples'].unique()
unique_patient = adata.obs['ID'].unique()

# Iterate over each unique sample
for sample in unique_patient:
    # Initialize a dictionary to store proportions for the current sample
    sample_proportions = {'ID': sample}
    
    # Iterate over each cell type
    for celltype in celltypes:
        # Filter adata for the current sample and celltype
        filtered_data = adata[(adata.obs['ID'] == sample) & (adata.obs['majority_voting'] == celltype)]
        
        # Count occurrences for the current sample and celltype
        count = len(filtered_data)
        print(count)
        # Calculate proportion for the current sample and celltype
        total_count_for_sample = len(adata[adata.obs['ID'] == sample])
        proportion = count / total_count_for_sample
        
        # Store proportion in the dictionary for the current celltype
        sample_proportions[celltype] = proportion
    
    # Append the dictionary of proportions for the current sample to the list
    proportions_list.append(sample_proportions)

# Convert the list of dictionaries into a DataFrame
proportions_df = pd.DataFrame(proportions_list)

# Set 'AllSamples' column as index (optional, if needed)
#proportions_df.set_index('ID', inplace=True)

# Create a new column 'Group' based on the first part of 'ID'
proportions_df['Group'] = proportions_df['ID'].str.split('_').str[0]


# Print the DataFrame
print(proportions_df)

selected_columns = [col for col in proportions_df.columns if col not in ['ID', 'Group']]

print(selected_columns)

proportions_df_long = proportions_df.melt(id_vars=['ID', 'Group'], value_vars=selected_columns, var_name='Cell Type', value_name='Proportion')

# Print the long format DataFrame
print(proportions_df_long)

kruskal_results = {}
for celltype, celltypes in proportions_df_long.groupby('Cell Type'):
    subgroup_data = celltypes.groupby('Group')['Proportion'].apply(list).values
    kruskal_result = stats.kruskal(*subgroup_data)
    kruskal_results[celltype] = kruskal_result

print(kruskal_results)

# Convert KruskalResult objects to dictionaries
kruskal_results_dict = {key: {'statistic': value.statistic, 'pvalue': value.pvalue} for key, value in kruskal_results.items()}

print(kruskal_results_dict)

sns.set(style='whitegrid')
fig, ax = plt.subplots(figsize=(18, 10))
sns.boxplot(x='Cell Type', y='Proportion', hue='Group', palette = color_discrete_map, data=proportions_df_long, )

plt.xticks(rotation=45)
ax.set_xlabel('Predicted Cell Type', fontweight='bold')
ax.set_ylabel('Proportion of Cells', fontweight='bold') 
plt.legend(title='Group', loc='center left', bbox_to_anchor=(1, 0.5))
plt.ylim(0,1)
plt.tight_layout()

# Save boxplot as file
filename_boxplot = os.path.join(dir, f"newProportion_Boxplot_300_dpi.png")
plt.savefig(filename_boxplot, dpi=300, format='png')
plt.close()
print(f"Boxplot saved to: {filename_boxplot}")

# Get unique values of 'AllSamples'
AllSamples = adata.obs['AllSamples'].unique()
print(f"Patient Groups: {AllSamples}")

n_idents = len(AllSamples)

figsize = 10
wspace = 0.25
# Adapt figure size based on number of rows and columns and added space between them
# (e.g. wspace between columns)
# Create subplots with 1 row and n_idents columns
fig, axes = plt.subplots(1, n_idents, figsize=(8 * (n_idents), 8))
plt.subplots_adjust(wspace=wspace)
# Plot UMAPs for each 'AllSamples' with the gene expression overlay
for i, AllSample in enumerate(AllSamples):
    subset = adata[adata.obs['AllSamples'] == AllSample]
    if i < n_idents - 1:
        sc.pl.umap(subset, color='majority_voting', add_outline=True, s=10, ax=axes[i], show=False, title=AllSample, legend_loc=None)
    else:
        sc.pl.umap(subset, color='majority_voting', add_outline=True, s=10, ax=axes[i], show=False, title=AllSample, legend_loc='right margin', legend_fontsize='xx-small')
        
    axes[i].set_xlabel("UMAP 1", fontweight='bold')
    axes[i].set_ylabel("UMAP 2", fontweight='bold')
    axes[i].set_title(f"{AllSample}", fontsize=18, fontweight='bold')
    axes[-1].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='x-small')

fig.suptitle('Majority Voting (CellTypist) Immune Cell Types', fontsize = 30, fontweight='bold')

# Save UMAP overlay as file
filename_umap = os.path.join(dir, f"newCellType_UMAP_Overlay_300dpi.png")
plt.savefig(filename_umap, dpi=300, format='png')

# Get unique values of 'AllSamples'
AllSamples = adata.obs['AllSamples'].unique()
print(f"Patient Groups: {AllSamples}")

n_idents = len(AllSamples)

figsize = 10
wspace = 0.25
# Adapt figure size based on number of rows and columns and added space between them
# (e.g. wspace between columns)
# Create subplots with 1 row and n_idents columns

fig, axes = plt.subplots(1, n_idents, figsize=(8 * (n_idents), 8))
plt.subplots_adjust(wspace=wspace, bottom=0.25)
# Plot UMAPs for each 'AllSamples' with the gene expression overlay
for i, AllSample in enumerate(AllSamples):
    subset = adata[adata.obs['AllSamples'] == AllSample]
    sc.pl.umap(subset, color='ID', add_outline=True, s=10, ax=axes[i], show=False, title=AllSample, legend_fontsize='x-small')
        
    axes[i].set_xlabel("UMAP 1", fontweight='bold')
    axes[i].set_ylabel("UMAP 2", fontweight='bold')
    axes[i].set_title(f"{AllSample}", fontsize=18, fontweight='bold')
    axes[i].legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), fontsize='x-small')

fig.suptitle('Immune Cell Types by Patient ID', fontsize = 30, fontweight='bold')

# Save UMAP overlay as file
filename_umap = os.path.join(dir, f"newCellType_ID_UMAP_Overlay_300dpi.png")
plt.savefig(filename_umap, dpi=300, format='png')

# Display average expression per gene against log2fc values
# Extract the gene expression matrix from adata.raw
gene_expression_matrix = adata.layers['filtered_counts'].toarray()

# Extract the gene names (assuming they are stored in adata.raw.var)
gene_names = adata.var.index

# Calculate mean expression for each gene across all cells
mean_expression = gene_expression_matrix.mean(axis=0)

# Create a DataFrame with gene names and mean expressions
mean_expression_df = pd.DataFrame({
    'Gene': gene_names,
    'MeanExpression': mean_expression
})

DEGfile_exp = pd.merge(DEGfile, mean_expression_df, on = 'Gene')

custom_params = {'axes.spines.right': False, 'axes.spines.top':False}
sns.set_theme(style='ticks', rc = custom_params)
plt.figure(figsize=(10, 6))

sns.scatterplot(x = 'MeanExpression', y = 'LogFoldChange', color = 'black', edgecolor = None, s = 30, data = DEGfile_exp)
sns.scatterplot(x = 'MeanExpression', y = 'LogFoldChange', color = '#bebada', edgecolor = None, s = 10, data = DEGfile_exp)
plt.ylabel("log2 FC")
plt.xlabel("Mean Expression")

# Save as figure
filename_exp_log2fc = os.path.join(dir, f"Expression_log2FC_300dpi.png")
plt.savefig(filename_exp_log2fc, dpi=300, format='png')

# Correlation between predicted cell type clusters
ax = sc.pl.correlation_matrix(adata, "majority_voting", cmap = 'RdBu_r', figsize=(9, 9))

filename_corrmatrix = os.path.join(dir, f"Correlation_Matrix_300dpi.png")
plt.savefig(filename_corrmatrix, dpi=300, format='png', bbox_inches='tight')

# Plot total cell number per cell type for each condition

# Initialize an empty list to store the results
results = []
celltypes = adata.obs['majority_voting'].unique()
groups = adata.obs['AllSamples'].unique()

# Iterate over each unique celltype
for celltype in celltypes:
    # Subset the data for the current celltype
    adata_celltype = adata[adata.obs['majority_voting'] == celltype]
    
    # Iterate over each unique group
    for group in groups:
        # Subset the data for the current group within the current celltype
        adata_celltype_group = adata_celltype[adata_celltype.obs['AllSamples'] == group]
        
        # Count the number of cells (rows) in this subset
        cell_count = adata_celltype_group.shape[0]
        
        # Append the result to the list
        results.append({'celltype': celltype, 'group': group, 'cell_total': cell_count})

# Convert the list of results into a DataFrame
df = pd.DataFrame(results)

# Display the resulting DataFrame
print(df)

plt.rcParams.update(plt.rcParamsDefault)
plt.figure(figsize=(12, 8))
sns.barplot(df, y = 'cell_total', x = 'celltype', palette=color_discrete_map, edgecolor = 'black', hue = 'group')
plt.xticks(rotation=80)
plt.xlabel('Predicted Cell Type', fontsize = 16, fontweight = 'bold')
plt.ylabel('Total Cell Count', fontsize = 16, fontweight = 'bold')
plt.legend(title = 'Condition', loc = 'center left', bbox_to_anchor = (1, 0.5))
plt.tight_layout()

# Save as figure
filename_exp_log2fc = os.path.join(dir, f"CellNumber_Celltype_300dpi.png")
plt.savefig(filename_exp_log2fc, dpi=300, format='png')

adata = sc.read_h5ad('../Input/Secondary/BloodCellsScores.h5ad')

# Plot naive helper T cell score
fig, ax = plt.subplots(figsize=(10, 10))
sc.pl.umap(adata, color = 'Tcm/Naive helper T cells_score', legend_loc = 'on data', legend_fontoutline=2, add_outline=True, s=10, ax=ax, show=False)
ax.set_xlabel("UMAP 1", fontweight='bold')
ax.set_ylabel("UMAP 2", fontweight='bold')
ax.set_title("Tcm/Naive Helper T Cells Score", fontsize=18, fontweight='bold')
plt.show()

# Save UMAP overlay as file
filename_umap = os.path.join(dir, f"NaiveScore_UMAP_Overlay_300dpi.png")
plt.savefig(filename_umap, dpi=300, format='png')
plt.close()

# Plot effector helper T cell score
fig, ax = plt.subplots(figsize=(10, 10))
sc.pl.umap(adata, color = 'Tem/Effector helper T cells_score', legend_loc = 'on data', legend_fontoutline=2, add_outline=True, s=10, ax=ax, show=False)
ax.set_xlabel("UMAP 1", fontweight='bold')
ax.set_ylabel("UMAP 2", fontweight='bold')
ax.set_title("Tem/Effector Helper T Cells Score", fontsize=18, fontweight='bold')
plt.show()

# Save UMAP overlay as file
filename_umap = os.path.join(dir, f"EffectorScore_UMAP_Overlay_300dpi.png")
plt.savefig(filename_umap, dpi=300, format='png')
plt.close()
