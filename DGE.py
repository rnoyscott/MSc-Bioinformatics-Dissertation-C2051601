# Import required modules
import scanpy as sc # for handling single cell RNA-seq data
import pandas as pd # for data manipulation
import os # for directory handling

# Load in the primary data
adata = sc.read_h5ad('../Input/Secondary/BloodCells.h5ad')
print(adata) # check data characteristics
print(adata.obs['AllSamples']) # ensure all condition groups are present

# Define subsets and their filenames using a dictionary
SubsetDictionary = {
    "Tcells": ['Tcm/Naive helper T cells', 'Tem/Effector helper T cells', 'Regulatory T cells', 'Tcm/Naive cytotoxic T cells', 'Tem/Trm cytotoxic T cells'],
    "CD4Tcells": ['Tcm/Naive helper T cells', 'Tem/Effector helper T cells'],
    "CD4TCM": ['Tcm/Naive helper T cells'],
    "CD8Tcells": ['Tcm/Naive cytotoxic T cells', 'Tem/Trm cytotoxic T cells'],
    "CD8TEM": ['Tem/Trm cytotoxic T cells'],
    "pDC": ["pDC"],
    "CD14mono": ['Classical monocytes'],
    "CD16mono": ['Non-classical monocytes'],
    "NK": ['CD16+ NK cells', 'CD16- NK cells', 'ILC3']
}

# Define comparisons to be made between condition groups
ComparisonDictionary = {
    'VE': ['Viraemic', 'Exposed'],
    'VC': ['Viraemic', 'Control'],
    'EC': ['Exposed', 'Control']
}

# Set the directory for storing subsets
dir = '../../../Input/Secondary/'

# Move to the Scanpy DGE directory for storage
os.chdir('../Output/DGE/Scanpy/')
# For each comparison run DGE analysis on the primary data
for Comparison, Groups in ComparisonDictionary.items():

    ## Filter to include only 'Exposed' and 'Control' samples
    adata_Comparison = adata[adata.obs['AllSamples'].isin(Groups)].copy()


    ## Extract group labels
    Group1, Group2 = Groups

    print(adata_Comparison.obs['AllSamples'])

    ## Carry out DEG calculations
    sc.tl.rank_genes_groups(adata_Comparison, groupby='AllSamples', layer = 'filtered_counts_norm', groups = [Group1], reference = Group2, use_raw = False, method='wilcoxon', tie_correct=True, pts=True)

    ## Retrieve DEGs from the 'adata' object
    result = adata_Comparison.uns['rank_genes_groups']
    print(result)

    ## Extracting DEG names and log fold changes
    gene_names = result['names'][Group1]
    logfoldchanges = result['logfoldchanges'][Group1]
    pvals = result['pvals'][Group1]
    adj_pvals = result['pvals_adj'][Group1]
    pts_1 = result['pts'][Group1]
    pts_2 = result['pts'][Group2]

    ## Create a DataFrame for easier manipulation
    deg_df = pd.DataFrame({
        'Gene': gene_names,
        'LogFoldChange': logfoldchanges,
        'p-value': pvals,
        'adjusted p-value':adj_pvals,
        'percentage group 1': pts_1,
        'percentage group 2': pts_2
    })

    print(f"DEGs for BloodCells: {deg_df}")
    print(f"Total Number of DEGs for BloodCells: {len(deg_df)}")

    deg_df = deg_df.sort_values(by='adjusted p-value') # order by significance

    # Specify the path to save the CSV file
    output_file = f'{Comparison}_DEG_BloodCells.csv'

    # Export DEGs to CSV
    deg_df.to_csv(output_file, index=False)

    print(f"Exported BloodCells DEGs to: '{output_file}'")

# Function to subset the primary data and perform DEG analysis
def Subsetting(adata, adata_Comparison, SubsetCells, filename, Group1, Group2, Comparison):

    Subset_1 = adata[adata.obs['majority_voting'].isin(SubsetCells)].copy() # subset to the current cell type

    # Save the subset file
    Subset_1.write_h5ad(os.path.join(dir, filename + '.h5ad'))
    print(f"Saved subset to: {filename}")

    # Subset cells for the current comparison
    Subset = adata_Comparison[adata_Comparison.obs['majority_voting'].isin(SubsetCells)].copy()

    # Perform DGE analysis
    sc.tl.rank_genes_groups(Subset, groupby='AllSamples', layer='filtered_counts_norm', 
                            groups=[Group1], reference=Group2, use_raw=False, 
                            method='wilcoxon', tie_correct=True, pts=True)
    print(f"DEGs calculated for: {filename}")

    # Retrieve DEGs from the 'Subset' object
    result = Subset.uns['rank_genes_groups']

    # Extract DEG details
    gene_names = result['names'][Group1]
    logfoldchanges = result['logfoldchanges'][Group1]
    pvals = result['pvals'][Group1]
    adj_pvals = result['pvals_adj'][Group1]
    pts_1 = result['pts'][Group1]
    pts_2 = result['pts'][Group2]

    # Create a DataFrame
    deg_df = pd.DataFrame({
        'Gene': gene_names,
        'LogFoldChange': logfoldchanges,
        'p-value': pvals,
        'adjusted p-value': adj_pvals,
        'percentage group 1': pts_1,
        'percentage group 2': pts_2
    })

    print(f"DEGs for {filename}: {deg_df}")
    print(f"Total Number of DEGs for {filename}: {len(deg_df)}")

    # Sort by adjusted p-value
    deg_df = deg_df.sort_values(by='adjusted p-value')

    # Specify the path to save the CSV file
    output_file = os.path.join(filename, f'{Comparison}_DEG_{filename}.csv')

    # Export DEGs to CSV
    deg_df.to_csv(output_file, index=False)

    print(f"Exported {filename} DEGs to: {output_file}")


# Main loop for comparisons and subsets
for Comparison, Groups in ComparisonDictionary.items():
    ## Extract group labels
    Group1, Group2 = Groups
    
    ## Subset the main data for the current comparison
    adata_Comparison = adata[adata.obs['AllSamples'].isin(Groups)].copy()

    ## Iterate over each subset configuration
    for SubsetName, SubsetCells in SubsetDictionary.items():
        Subsetting(adata, adata_Comparison, SubsetCells, SubsetName, Group1, Group2, Comparison)