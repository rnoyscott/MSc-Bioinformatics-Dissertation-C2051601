# Functions for single cell analysis

# Function for locating groups of files, specifically h5ad and csv files
def FileLocator(input1, input2, celltype, seurat):
    # Comprehensive list of all modules used across every function
    import scanpy as sc # Scanpy is the primary module for single cell RNAseq analysis in Python
    import matplotlib.pyplot as plt # for producing customisable figures 
    import seaborn as sns # for producing customisable figures
    import pandas as pd # for data manipulation
    import numpy as np # for data manipulation
    import os # for command line arguments
    import glob # for file name pattern matching
    import re # for regular expressions
    import gseapy as gp # for ORA
    from gseapy.scipalette import SciPalette # for customisation of ORA plots
    from scipy import stats # for statistics

    adata_file = os.path.join(input1, f'{celltype}.h5ad') # concatenate the input directory with a cell type data file
    adata1_file = os.path.join(input1, 'BloodCells.h5ad') # read in the primary data file named BloodCells
    DEG_file = os.path.join(input2, f'VE_DEG_{celltype}.csv')
    #DEG_file = glob.glob(os.path.join(input2, f'*DEG*{celltype}*.csv')) # concatenate the input directory with a cell type DEG file


    adata = sc.read_h5ad(adata_file) # read in the subset data file
    adata1 = sc.read_h5ad(adata1_file) # read in the primary data file

    DEGfile = pd.read_csv(DEG_file) # read in the DEG file

    # if the DEG is the result of Seurat::FindMarkers, then reformat it to be compatible with the workflow
    if seurat:
        new_column_names = ["Gene", "p-value", "LogFoldChange", "percentage group 1", "percentage group 2", "adjusted p-value"]
        DEGfile.columns = [new_column_names[i] if i < len(new_column_names) else col for i, col in enumerate(DEGfile.columns)]

    return adata, adata1, DEGfile

# Define a function for Single Gene Analysis
def SingleGeneAnalysis(gene, cellname, width, height, dpi, h5adObject, BloodCells, output):
    # Comprehensive list of all modules used across every function
    import scanpy as sc
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    import os
    import glob
    import re
    import gseapy as gp
    from gseapy.scipalette import SciPalette
    from scipy import stats

    color_discrete_map= {'Viraemic': '#c38961', 'Control': 'lightgrey', 'Exposed': '#007d82'} # create a colour mapper for each patient group

    # Produce a UMAP overlay of BloodCells, split by condition
    AllSamples = BloodCells.obs['AllSamples'].unique() # get the patient group IDs
    print(f"Patient Groups: {AllSamples}") # print them to output

    n_idents = len(AllSamples) # set the number of plots to be created by using the number of patient groups

    fig, axes = plt.subplots(1, n_idents, figsize=(width * n_idents, height)) # define the plotting parameters

    ## Setting the maximum gene expression value depending on the input data and the highest expression value identified across condition groups
    max_gene_expression = 0 # begin by setting to 0
    for AllSample in AllSamples:
        subset = BloodCells[BloodCells.obs['AllSamples'] == AllSample] # subset to the current condition group being checked
        max_expression_in_subset = subset[:, subset.var.index == gene].layers['filtered_counts'].max() # use the filtered counts layer
        max_gene_expression = max(max_gene_expression, max_expression_in_subset) # update value depending on the maximum expression value for the current condition group

    ## Plot UMAPs for each 'AllSamples' with the gene expression overlay
    for i, AllSample in enumerate(AllSamples):
        subset = BloodCells[BloodCells.obs['AllSamples'] == AllSample] # subset the primary dataset to the condition group being iterated over
        sc.pl.umap(BloodCells, palette='lightgrey', ax=axes[i], show=False, s = 10) # plot the entire UMAP for BloodCells in grey
        sc.pl.umap(subset, color=gene, vmin=0, vmax=max_gene_expression, title=AllSample, ax=axes[i], show=False, s = 10, color_map='viridis', add_outline=True) # overlay with the UMAPs for each group
        
        ### Customise the axes
        axes[i].set_xlabel("UMAP 1", fontweight='bold')
        axes[i].set_ylabel("UMAP 2", fontweight='bold')
        axes[i].set_title(f"{AllSample}", fontsize=18, fontweight='bold')
        
    ## Save UMAP overlay as file
    filename_umap = os.path.join(output, f"{gene}_UMAP_Overlay_{dpi}dpi.png") # define a name for the output figure
    plt.savefig(filename_umap, dpi=dpi, format='png') # save the figure
    plt.close() # close the figure
    print(f"UMAP Overlay for {gene} saved to: {filename_umap}") # print a message to indicate where the figure has been saved to

    # Produce a violin plot of gene expression across groups
    ## Define figure parameters and plot
    custom_params = {'axes.spines.right': False, 'axes.spines.top':False} # define custom figure style
    sns.set_theme(style='ticks', rc = custom_params) # apply custom style to theme

    ## Produce the plot
    sc.pl.violin(h5adObject, gene, palette = color_discrete_map, groupby='AllSamples', stripplot=False, inner='box') # plot data for a given gene

    ## Customise the axes
    plt.xlabel('Condition Group', fontweight = 'bold')
    plt.ylabel(f'{gene} Expression', fontweight = 'bold')

    ## Save violin plot as file
    filename_violin = os.path.join(output, f"{gene}_Violin_Plot_{dpi}dpi.png")
    plt.savefig(filename_violin, dpi=dpi, format='png')
    plt.close()
    print(f"Violin plot for {gene} saved to: {filename_violin}")

    # Produce a boxplot displaying gene expression across patients between two groups
    ## First create pandas dataframe containing the necessary information
    gene_expression = h5adObject[:, h5adObject.var.index == gene].layers['filtered_counts'].toarray().flatten() # extract the gene expression data from the filtered countd layer
    
    ## Create the dataframe
    df = pd.DataFrame({
        'Gene Expression': gene_expression,
        'Group': h5adObject.obs['AllSamples'].values,
        'ID': h5adObject.obs['ID'].values
    })

    ## Define figure parameters and plot
    plt.figure(figsize=(width, height)) # using parameters that are user-supplied

    custom_params = {'axes.spines.right': False, 'axes.spines.top':False} # define custom figure style
    sns.set_theme(style='ticks', rc = custom_params) # apply custom style to theme

    ## Create the plot
    sns.boxplot(x='ID', y='Gene Expression', hue='Group', palette = color_discrete_map, data=df)

    ## Customise the axes
    plt.legend(title = 'Condition', loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('Patient ID', fontweight = 'bold')
    plt.ylabel(f'{gene} Expression', fontweight = 'bold')
    plt.xticks(rotation=45)    
    plt.tight_layout()

    ## Save boxplot as file
    filename_boxplot = os.path.join(output, f"{gene}_Expression_Boxplot_{dpi}dpi.png")
    plt.savefig(filename_boxplot, dpi=dpi, format='png')
    plt.close()
    print(f"Boxplot for {gene} saved to: {filename_boxplot}")

# Define a function for multiple gene analysis
def MultipleGeneAnalysis(pathwayname, cellname, genes, h5adObject, DEGfile, output, dpi=300, width=10, height=8, title=None, seurat = False):
    # Comprehensive list of all modules used across every function
    import scanpy as sc
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    import os
    import glob
    import re
    import gseapy as gp
    from gseapy.scipalette import SciPalette
    from scipy import stats

    # Set the figure size for all plots
    plt.rcParams["figure.figsize"] = (width, height)

    # For subsetting DEG file and ranking genes by absolute log2FC
    DEGfilegenes = DEGfile[DEGfile['Gene'].isin(genes)]
    DEGfilegenes['AbsLogFoldChange'] = DEGfilegenes['LogFoldChange'].abs() # convert negative values to positive for sorting
        
    DEGfilegenes = DEGfilegenes.sort_values(by='AbsLogFoldChange', ascending=False) # sort the subset data by absolute log2FC

    # If more than 20 genes then keep the top 20 ranked genes for visualisation
    if len(genes) > 20:    
        genes_top_20 = DEGfilegenes['Gene'].head(20).tolist() # get the top 10 genes based on absolute log2FC
        
        genes = genes_top_20 # update the genes list to the top 20 genes
    else:
        pass # if less than 20 genes, skip this filtering step
    try:
        # Produce heatmap grouped by patient group

        ## Set global parameters for plotting
        plt.rcParams['font.weight'] = 'bold'
        plt.rcParams['font.style'] = 'italic'
        plt.rcParams['axes.labelweight'] = 'bold'
        plt.rcParams['axes.titleweight'] = 'bold'

        ## Produce heatmap of genes with colour corrosponding to mean expression
        mp = sc.pl.matrixplot(h5adObject, genes, 'AllSamples', return_fig=True)
        mp.style(edge_color='black', cmap='PuBu', edge_lw = 0.9).show() # further customisation of heatmap
        plt.tight_layout()

        ## Save heatmap
        filename_heatmap = os.path.join(output, f"{pathwayname}_{cellname}_Group_Expression_Heatmap_{dpi}dpi.png")
        plt.savefig(filename_heatmap, dpi=dpi, format='png', bbox_inches='tight')
        plt.close()
        print(f"Expression by patient group for {pathwayname} in {cellname}, saved to: {filename_heatmap}")
    except Exception as e:
        print(f'Failed to create heatmap of genes by patient group: {e}') # provide explanation if failure to produce heatmap
    
    try:
        # Perform rank_genes_groups analysis (grouped by patient group) and plot dotplot
        ## Compute the log2FC values based on normalised counts
        sc.tl.rank_genes_groups(h5adObject, 'AllSamples', n_genes=h5adObject.layers['filtered_counts_norm'].shape[1], use_raw=False)

        ## Plot the log2FC values for each gene
        dp = sc.pl.rank_genes_groups_dotplot(
        h5adObject,
        var_names=genes,
        values_to_plot="logfoldchanges", cmap='bwr',
        vmin=-4,
        vmax=4,
        #min_logfoldchange=3,
        colorbar_title='log2 FC',
        dendrogram=False
        )

        plt.tight_layout()

        ## Save dotplot
        filename_fcheatmap = os.path.join(output, f"{pathwayname}_{cellname}_FoldChange_Heatmap_{dpi}dpi.png")
        plt.savefig(filename_fcheatmap, dpi=dpi, format='png', bbox_inches='tight')
        plt.close()
        print(f"Differential expression heatmap for {pathwayname} in {cellname}, saved to: {filename_fcheatmap}")
    except Exception as e:
        print(f'Failed to create heatmap of gene fold change by patient group: {e}')  # provide explanation if failure
    
    try:
        # Produce a volcano plot of DEGs and label specified genes
        ## Reset the global plotting parameters
        plt.rcParams.update(plt.rcParamsDefault)

        ## Define the color for each point based on the conditions
        DEGfile['color'] = 'grey'  # default color which is overwitten if thresholds are met

        ### Blue for logfoldchange < -0.25 & adjusted p value < 0.05
        DEGfile.loc[(DEGfile["LogFoldChange"] < -0.25) & (DEGfile["adjusted p-value"] < 0.05), 'color'] = 'blue'

        ### Red for logfoldchange > 0.25 & adjusted p value < 0.05
        DEGfile.loc[(DEGfile["LogFoldChange"] > 0.25) & (DEGfile["adjusted p-value"] < 0.05), 'color'] = 'red'

        ### Convert adjusted p-value using log
        DEGfile["-logQ"] = -np.log(DEGfile["adjusted p-value"].astype("float"))

        ### Plot the data as volcano plot
        fig, ax = plt.subplots()

        ax.scatter(
            x=DEGfile["LogFoldChange"],
            y=DEGfile["-logQ"],
            c=DEGfile["color"],
            s=10,
            alpha = 0.7
        )

        #### Customise the axes
        ax.set_xlabel("log2 FC", fontweight = 'bold')
        ax.set_ylabel("-log Q-value", fontweight = 'bold')

        #### Overlay with specified genes and label with gene symbols
        last_x, last_y = None, None # set the x and y coordinates

        genes_top_5 = DEGfilegenes['Gene'].head(5).tolist() # only annotate the top 5 genes according to absolute log2FC

        ##### For each row in the DEG file highlight the genes that are in the gene list
        for _, row in DEGfile.iterrows():
            if row['Gene'] in genes:
                ax.scatter(row['LogFoldChange'], row['-logQ'], c=row["color"], edgecolor='black', s=10)
                
                ###### Generate random offsets for x and y to be used in annotating
                x_offset = np.random.choice(np.concatenate((np.arange(-50, -10), np.arange(10, 50))))
                y_offset = np.random.choice(np.concatenate((np.arange(-15, -5), np.arange(5, 15))))

                ####### Ensure that current annotation does not overlap the last annotation and ensure annotations do not exceed plot border
                if last_x is not None and last_y is not None:
                    if abs(row['LogFoldChange'] - last_x) < 0.1 and abs(row['-logQ'] - last_y) < 40:
                        x_offset += 15 if x_offset > 0 else -15
                        y_offset += 15 if y_offset > 0 else -15
                    
                proposed_y_label = row['-logQ'] + y_offset
                if proposed_y_label < 5:
                    y_offset = -row['-logQ']

                ###### Annotate each of the top 5 genes
                if row['Gene'] in genes_top_5:
                    ax.annotate(
                        row['Gene'], 
                        xy=(row['LogFoldChange'], row['-logQ']), # plotting coordinates
                        xytext=(x_offset, y_offset),
                        textcoords='offset points',
                        ha='left', fontsize=8, style='italic', fontweight='bold', # customise appearance of annotation
                        arrowprops=dict(arrowstyle='-', color='black', lw=0.5, shrinkA=0, shrinkB=0)  # add line connecting annotation to point
                    )

                    last_x, last_y = row['LogFoldChange'], row['-logQ'] # update the latest annotation coordinates

        #### Add threshold lines
        ax.axhline(-np.log(0.05), c='k', lw=2, ls='--')
        ax.axvline(0.25, c='k', lw=2, ls='--')
        ax.axvline(-0.25, c='k', lw=2, ls='--')
        #### If the DEG file is from Seurat FindMarkers then limit the x axes
        if seurat:
            plt.xlim(-1, 1)

        #### If a title is not defined by user, automatically create one using the user-supplied pathway name
        if title is None:
            pathwayname_sep = pathwayname.replace("_", " ")
            pathwayname_sep_title = pathwayname_sep.title()
            title = f"{pathwayname_sep_title} in {cellname}"
        
        plt.title(title, fontweight = 'bold')
        plt.tight_layout()

        #### Save volcano plot
        filename_volcano = os.path.join(output, f"{pathwayname}_{cellname}_Volcano_Plot_{dpi}dpi.png")
        plt.savefig(filename_volcano, dpi=dpi, format='png')
        plt.close()
        print(f"Volcano plot for {pathwayname} in {cellname}, saved to: {filename_volcano}")
    except Exception as e:
        print(f'Failed to create volcano plot of genes in {cellname}: {e}') # print errors if they occur

# Define function for subsetting a list od DEGs using genes for a particular pathway
def SubsetDEGs(DEGfile, path_to_pathway_results, go_terms=None):
    # Comprehensive list of all modules used across every function
    import scanpy as sc
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    import os
    import glob
    import re
    import gseapy as gp
    from gseapy.scipalette import SciPalette
    from scipy import stats

    try:
        # Extract genes for a specific GO term and create a list
        pathways = pd.read_csv(path_to_pathway_results, sep=',')  # read in pathway file
        pathways['Term'] = pathways['Term'].apply(lambda x: re.search(r'\((.*?)\)', x).group(1)) # convert 'Term' column so that it only contains the GO accesion ID
        pathways_filt = pathways[pathways['Term'].isin(go_term)]  # filter for rows that match the GO ID
        genes = pathways_filt['Genes'].tolist() # convert the gene column to a list of genes
        genes = [gene for sublist in genes for gene in sublist.split(';')] # refine the list
        
        # Subset DEGfile based on selected genes
        subset_data = DEGfile[DEGfile['Gene'].isin(genes)]
        genes = subset_data['Gene'].tolist()
    except Exception as e:
        print(f'Failed to identify genes using {go_term}: {e}') # print errors if they occur
    
    # Deal with excessively long gene lists as these are an issue when plotting figures
    if len(genes) > 20:
        subset_data['AbsLogFoldChange'] = subset_data['LogFoldChange'].abs() # convert negative log2FC values to positive for sorting
        
        subset_data = subset_data.sort_values(by='AbsLogFoldChange', ascending=False) # sort the subset data by absolute log2FC
        
        genes_top_20 = subset_data['Gene'].head(20).tolist() # get the top 20 genes based on absolute log2FC
        
        genes = genes_top_20 # update the genes list to the top 20 genes
    return genes

# Define function for ORA
def Pathways (cellname, width, height, dpi, comparison, DEGfile, path, h5adObject, output):
    # Comprehensive list of all modules used across every function
    import scanpy as sc
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    import os
    import glob
    import re
    import gseapy as gp
    from gseapy.scipalette import SciPalette
    from scipy import stats
    
    # Create empty vectors for storing file names
    filename_enr_up = None
    filename_enr_dw = None

    # Create background gene set
    
    ## Genes from experiment
    genes_exp = h5adObject.var_names.tolist() # convert gene symbols to a list
    print(f"Number of experimental genes: {len(genes_exp)}") # check how many there are
    
    ## Genes from GSEA C7 Immunological datasets
    genes_C7 = [] # create an empty vector which will store the genes

    files = glob.glob(os.path.join(path, 'GSE22229*')) # read in all files from the background directory that are prefixed by the unique ID for C7 datasets
    print(f"Number of GSE files found: {len(files)}") # check how many there are

    ### For each of the located files, extract the gene symbols
    for file in files:
        df = pd.read_csv(file, sep='\t', skiprows=17, header=None)
        genes = df.iloc[0, 1].split(',')
        print(f'Number of genes in {file}: {len(genes)}')
        genes_C7.extend(genes) # append to the vector created earlier

    print(f"Number of genes from GSE files: {len(genes_C7)}")

    genes_C7_unique = list(set(genes_C7)) # remove duplicate genes

    ## Combine the sets to get unique genes from both sources
    combined_genes_unique = list(set(genes_exp).union(genes_C7_unique))
    print(f"Number of unique combined genes: {len(combined_genes_unique)}")

    # Subset up or down regulated genes
    degs_sig = DEGfile[DEGfile['adjusted p-value'] < 0.05] # filter for significant DEGs
    degs_up = degs_sig[degs_sig.LogFoldChange > 0] # define upregulated DEGs
    degs_dw = degs_sig[degs_sig.LogFoldChange < 0] # define downregulated DEGs

    print(f"Number of significant DEGs: {len(degs_sig)}")
    print(f"Number of upregulated DEGs: {len(degs_up)}")
    print(f"Number of downregulated DEGs: {len(degs_dw)}")

    # Enrichr API for performing ORA

    ## Only perform if there were significantly upregulated DEGs found
    if len(degs_up) > 0:
        enr_up = gp.enrichr(gene_list=degs_up.Gene.tolist(), # genes from the significantly upregulated genes
                            gene_sets=['GO_Biological_Process_2023'], # the most recent library of GO BPs
                            outdir=None,
                            background=combined_genes_unique, # the background gene set defined earlier
                            cutoff = 0.5) # a threshold for pathway analysis

        enr_up_df = enr_up.results # store the results in a dataframe

        print(f"Number of upregulated pathways: {len(enr_up_df)}") # check how many pathways were identified

        ### Save as a CSV file
        filename_enr_up = os.path.join(output, f"{cellname}_{comparison}_ORA_upregulated.csv")
        enr_up_df.to_csv(filename_enr_up)

        ## Dotplot of upregulated pathways if the pathways dataframe is not empty
        if not enr_up_df.empty:
            
            ### Define plotting parameters
            fig, ax = plt.subplots(figsize=(width, height))

            ### Plot a lenient threshold of 0.5 and save to the specified output directory
            gp.dotplot(enr_up_df, title=f"Upregulated Pathways {comparison}", cmap=plt.cm.autumn_r, cutoff=0.9, ofname=os.path.join(output, f"{cellname}_{comparison}ORA_Plot_Upregulated_{dpi}dpi.png"))
            plt.close()
    else:
        print("No upregulated genes, enrichment analysis not run.")
    
    # Identical process but for downregulated DEGs
    if len(degs_dw) > 0:
        enr_dw = gp.enrichr(gene_list=degs_dw.Gene.tolist(),
                            gene_sets=['GO_Biological_Process_2023'],
                            outdir=None,
                            background=combined_genes_unique,
                            cutoff = 0.5)
        enr_dw_df = enr_dw.results
        print(f"Number of downregulated pathways: {len(enr_dw_df)}")
        filename_enr_dw = os.path.join(output, f"{cellname}_{comparison}_ORA_downregulated.csv")
        enr_dw_df.to_csv(filename_enr_dw)
        # Dotplot of downregulated pathways
        if not enr_dw_df.empty:
            fig, ax = plt.subplots(figsize=(width, height))
            gp.dotplot(enr_dw_df, title=f"Downregulated Pathways {comparison}", cmap=plt.cm.winter_r, cutoff=0.9, ofname=os.path.join(output, f"{cellname}_{comparison}_ORA_Plot_Downregulated_{dpi}dpi.png"))
    else:
        print("No downregulated genes, enrichment analysis not run.")

    # Print saved file paths
    print(f"Pathways for {cellname} saved to:")
    if filename_enr_up:
        print(f"  Upregulated: {filename_enr_up}")
    if filename_enr_dw:
        print(f"  Downregulated: {filename_enr_dw}")