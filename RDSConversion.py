import anndata
import scipy
import pandas as pd

# Output directory for metadata
dir = '../Resources/Metadata/'
# Load the AnnData object
adata = anndata.read_h5ad("../Input/Secondary/BloodCells.h5ad")

# Convert index objects to lists
cell_names = adata.obs_names.tolist()
gene_names = adata.var_names.tolist()

# Save to CSV files for R to read
pd.DataFrame({
    'cell_names': cell_names
}).to_csv(os.path.join(dir, 'cell_names.csv', index=False))

pd.DataFrame({
    'gene_names': gene_names
}).to_csv(os.path.join(dir, 'gene_names.csv', index=False))