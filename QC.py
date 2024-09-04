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

# Create a QC directory to store outputs
dir = '../Output/QC/'
os.makedirs(dir, exist_ok=True)

# Set global parameters for plotting
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.style'] = 'italic'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'

# Functions
## Filename formatting function
def create_filename(name):
    ### Replace invalid characters with underscores
    return re.sub(r'[\/\\:*?"<>|]', '_', name)

## Gene name formatting function
def clean_gene_names(gene_list):
    return [gene.strip() for gene in gene_list]

# Read in the data
adata = sc.read_h5ad('../Input/Secondary/BloodCells.h5ad')

# Manually check clusters using a file containing marker genes

## Read in the file
celltype_markers = pd.read_csv('../Resources/Markers/markers.csv')

## Create dictionary of labels mapped to marker genes
celltype_dict = celltype_markers.set_index('ref_label')['marker_genes'].to_dict()

## Plot UMAP and violin plots of gene expression for each marker gene corrosponding to each cell type
for celltype, marker_genes in celltype_dict.items():
    celltype_filename = create_filename(celltype)

    ### Convert marker_genes string to a list of gene names
    marker_genes_list = [gene.strip() for gene in marker_genes.split(',')]

    ### Check if marker_genes are present in adata
    present_marker_genes = [gene for gene in marker_genes_list if gene in adata.var_names]
    print(f"Marker genes present in adata for {celltype}:", present_marker_genes)

    if present_marker_genes:
        #### Plot UMAP
        sc.pl.umap(adata, color=present_marker_genes, layer = 'filtered_counts', s=10, add_outline=True)
        filename_umap = os.path.join(dir, f'{celltype_filename}_Markers_UMAP_300dpi.png')
        plt.savefig(filename_umap, dpi=300, format='png', bbox_inches='tight')
        plt.close()


        #### Plot Violin plots
        sc.pl.violin(adata, present_marker_genes, layer = 'filtered_counts', groupby='majority_voting', rotation=90)
        filename_violin = os.path.join(dir, f'{celltype_filename}_Markers_Violin_300dpi.png')
        plt.savefig(filename_violin, dpi=300, format='png', bbox_inches='tight')
        plt.close()
    else:
        print(f"None of the marker genes for {celltype} are present in adata.")

## Set a single marker gene for each cell type to be plotted in a matrix plot
markers = ['CD19', 'FCGR3A', 'GNLY', 'CD14', 'IL1R1', 'CD34', 'ITGA2B', 'MS4A1', 'FCGR3A', 'FOXP3', 'LEF1', 'CD4', 'IL7R', 'CCL5', 'IL3RA']

sc.pl.dotplot(adata, markers, layer = 'filtered_counts', groupby='majority_voting')
plt.tight_layout()

filename_dotplot = os.path.join(dir, 'Predicted_CellType_UniqueMarkers_Expression_Heatmap_300dpi.png')
plt.savefig(filename_dotplot, dpi=300, format='png', bbox_inches='tight')
plt.close()

print(adata.obs.columns)
labels = adata.obs[['ref_label', 'ref_label_pbmc3k', 'majority_voting', 'overcluster']].groupby('overcluster').agg(lambda x: x.mode()[0] if not x.mode().empty else pd.NA)
scores = adata.obs[['ref_score', 'ref_score_pbmc3k', 'overcluster']].groupby('overcluster').agg(lambda x: x.mean())
cluster_mapping = labels.merge(right = scores, left_index = True, right_index = True)

csv_filename = os.path.join(dir, 'cluster_predictions.csv')
cluster_mapping.to_csv(csv_filename)
print(cluster_mapping)

celltype_dict = celltype_markers.set_index('ref_label')['marker_genes'].to_dict()

# Initialize an empty DataFrame for merged results
merged_results = pd.DataFrame()

# Score genes for each cell type
for celltype, marker_genes in celltype_dict.items():
    marker_genes = clean_gene_names(marker_genes.split(","))  # assuming marker_genes are comma-separated
    # Check if all marker genes are present in adata.var_names
    marker_genes = [gene for gene in marker_genes if gene in adata.var_names]
    
    if marker_genes:
        sc.tl.score_genes(adata, marker_genes, score_name=f'{celltype}_score')
        scores = adata.obs[['majority_voting', f'{celltype}_score']].groupby('majority_voting').median().reset_index()
        scores = scores.rename(columns={f'{celltype}_score': f'{celltype}_median_score'})
        
        if merged_results.empty:
            merged_results = scores
        else:
            merged_results = merged_results.merge(scores, on='majority_voting', how='outer')
    else:
        print(f"No valid marker genes found for cell type: {celltype}")

# If necessary, reset index to get a clean DataFrame
merged_results = merged_results.reset_index(drop=True)
print(f'Printing merged_results: {merged_results}')

# Save the DataFrame to CSV
csv_filename = os.path.join(dir, 'scoresv1.csv')
merged_results.to_csv(csv_filename, index=False)

adata.write_h5ad('../Input/Secondary/BloodCellsScores.h5ad')