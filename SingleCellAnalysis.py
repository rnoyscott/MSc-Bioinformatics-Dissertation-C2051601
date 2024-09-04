# Import modules
import argparse # for creating command line arguments
import os # for input/output directory handling
import glob # for reading in files
import scanpy as sc # for reading in the h5ad data files containing single cell data
import pandas as pd # for dataframe manipulation

# Import custom functions stored in Functions.py
from Functions import FileLocator, SingleGeneAnalysis, MultipleGeneAnalysis, Pathways, SubsetDEGs

# Call the custom functions and define them inside 'main'
def main(args):
    # Uses the output argument to check whether an output directory already exists, if not then create one
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Procedure if pathways is selected from function calls
    if args.command == 'pathways':
        isFile = os.path.isfile(args.input1) # check whether the input is a directory or a file name
        
        ## Initialise vectors for storing data
        dataset = {}
        DEGfile = {}

        ## If the input is a file, do the following:
        if isFile:
            data_name = os.path.basename(args.input1).replace('.h5ad', '') # create a data ID name
            data = sc.read_h5ad(args.input1) # read in the data file
            DEGfile = pd.read_csv(args.input2) # read in the DEG file

            ### Use the Pathways function
            print(f"Calculating pathways for {data_name}...")
            Pathways(data_name, args.width, args.height, args.dpi, DEGfile, args.background, data, args.output)
            print(f"{data_name} completed.")
        
        ## If the input is a directory, do the following:
        else:
            ### Create a list of files in the input directory that are suffixed by .h5ad (the RDS equivalent for python)
            files = glob.glob(os.path.join(args.input1, '*.h5ad'))
            print(f"Number of datasets found: {len(files)}")

            ### Read each file in
            datasets = {}
            for file in files:
                var_name = os.path.basename(file).replace('.h5ad', '') # extract the base name and remove the '.h5ad' extension
                
                data = sc.read_h5ad(file) # read in the data file
                
                datasets[var_name] = data # pass the base name to a vector for storage

            ## The same procedure but for the corrosponding DEG files
            files = glob.glob(os.path.join(args.input2, '*DEG*.csv'))
            print(f"Number of DEG files found: {len(files)}")

            DEGfiles = {}

            for file in files:
                basename = os.path.basename(file)
                baseparts = basename.split('_')
                var_name = baseparts[2].split('.')[0]
                # Comparison suffix to use in naming output files
                comparison = file.split(os.path.sep)[-1].split('_')[0]
                
                DEGfile = pd.read_csv(file)
                if args.seurat:
                    new_column_names = ["Gene", "p-value", "LogFoldChange", "percentage group 1", "percentage group 2", "adjusted p-value"]
                    DEGfile.columns = [new_column_names[i] if i < len(new_column_names) else col for i, col in enumerate(DEGfile.columns)]
                
                DEGfiles[var_name] = DEGfile

            ## For each pair of dataset and DEG file, perform the pathways function
            for dataset_name in datasets:
                ### Only if a pair is identified
                if dataset_name in DEGfiles:
                    print(f"Calculating pathways for {dataset_name}...")
                    Pathways(dataset_name, comparison, DEGfiles[dataset_name], args.background, datasets[dataset_name], args.output)
                    print(f"{dataset_name} completed.")
                ### If there is a mismatch, do not run pathways for the current data file and print an error message
                else:
                    print(f"DEG file for {dataset_name} not found.")
    
    # Procedure if subset is selected from function calls
    if args.command == 'subset':
        adata, adata1, DEGfile = FileLocator(args.input1, args.input2, args.comparison, args.celltype, args.seurat) # read in the data files that are specific to the given cell type

        genes = SubsetDEGs(DEGfile, args.input3, go_term=args.goterm) # get a subset of genes from the DEG file that belong to the GO term
        print(genes)

        MultipleGeneAnalysis(args.pathwaytitle, args.celltype, genes, adata, DEGfile, output=args.output, dpi=args.dpi, width=args.width, height=args.height, seurat=args.seurat) # perform multiple gene analysis on the cell type and genes identified using the GO term

    # Procedure if a manual list of genes is passed
    if args.command == 'multiple':
        if args.genes:
            adata, adata1, DEGfile = FileLocator(args.input1, args.input2, args.comparison, args.celltype, args.seurat) # read in the data files that are specific to the given cell type
            MultipleGeneAnalysis(args.pathwaytitle, args.celltype, args.genes, adata, DEGfile, output=args.output, dpi=args.dpi, width=args.width, height=args.height, seurat=args.seurat) # perform multiple gene analysis on the cell type and genes that were manually passed

    # Procedure if single is selected from function calls
    if args.command == 'single':
        if args.gene:
            adata, adata1, DEGfile = FileLocator(args.input1, args.input2, args.comparison, args.celltype, args.seurat) # read in the data files that are specific to the given cell type

            SingleGeneAnalysis(args.gene, args.celltype, width=args.width, height=args.height, dpi=args.dpi, h5adObject=adata, BloodCells=adata1, output=args.output) # perform single gene analysis on the cell type and gene of interest

# Define the command line argument structure
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='A comprehensive single cell RNA-seq analysis pipeline.')
    subparsers = parser.add_subparsers(dest='command', help='Sub-commands')

    # Pathways subparser
    pathways_parser = subparsers.add_parser('pathways', help='Compute up/downregulated pathways for different immune cell types.')
    pathways_parser.add_argument('-b', '--background', type=str, required=True, help='Path to directory containing background gene sets.')
    pathways_parser.add_argument('-i1', '--input1', type=str, required=True, help='Path to input directory containing datasets.')
    pathways_parser.add_argument('-i2', '--input2', type=str, required=True, help='Path to input directory containing DEG files.')
    pathways_parser.add_argument('-o', '--output', type=str, required=True, help='Path to output directory.')
    pathways_parser.add_argument('-srt', '--seurat', action='store_true', help='Whether DEGs were calculated with Seurat.')

    # Subset subparser
    subset_parser = subparsers.add_parser('subset', help='Subset analysis using pathway enrichment results.')
    subset_parser.add_argument('-i1', '--input1', type=str, required=True, help='Path to input directory containing datasets.')
    subset_parser.add_argument('-i2', '--input2', type=str, required=True, help='Path to input directory containing DEG files.')
    subset_parser.add_argument('-i3', '--input3', type=str, help='Path to input directory containing pathway results.')
    subset_parser.add_argument('-o', '--output', type=str, required=True, help='Path to output directory.')
    subset_parser.add_argument('-wd', '--width', type=int, default=10, help='Width of output figures.')
    subset_parser.add_argument('-ht', '--height', type=int, default=8, help='Height of output figures.')
    subset_parser.add_argument('-d', '--dpi', type=int, default=300, help='DPI of output figures.')
    subset_parser.add_argument('-pt', '--pathwaytitle', type=str, required=True, help='Title of pathway, used to name associated outputs.')
    subset_parser.add_argument('-ct', '--celltype', type=str, required=True, help='Title of cell type.')
    subset_parser.add_argument('-gt', '--goterm', nargs='+', help='GO accession ID.')
    subset_parser.add_argument('-srt', '--seurat', action='store_true', help='Whether the DEGs were calculated with Seurat.')
    subset_parser.add_argument('-cmp', '--comparison', type=str, required=True, help='Groupwise comparison suffix (either VE, VC, or EC).')
    
    # Multiple subparser
    multiple_parser = subparsers.add_parser('multiple', help='Multiple gene analysis.')
    multiple_parser.add_argument('-i1', '--input1', type=str, required=True, help='Path to input directory containing datasets.')
    multiple_parser.add_argument('-i2', '--input2', type=str, required=True, help='Path to input directory containing DEG files.')
    multiple_parser.add_argument('-o', '--output', type=str, required=True, help='Path to output directory.')
    multiple_parser.add_argument('-wd', '--width', type=int, default=10, help='Width of output figures.')
    multiple_parser.add_argument('-ht', '--height', type=int, default=8, help='Height of output figures.')
    multiple_parser.add_argument('-d', '--dpi', type=int, default=300, help='DPI of output figures.')
    multiple_parser.add_argument('-pt', '--pathwaytitle', type=str, required=True, help='Title of pathway, used to name associated outputs.')
    multiple_parser.add_argument('-ct', '--celltype', type=str, required=True, help='Title of cell type.')
    multiple_parser.add_argument('-gs', '--genes', nargs='+', help='List of genes for analysis.')
    multiple_parser.add_argument('-srt', '--seurat', action='store_true', help='Whether the DEGs were calculated with Seurat.')
    multiple_parser.add_argument('-cmp', '--comparison', type=str, required=True, help='Groupwise comparison suffix (either VE, VC, or EC).')

    # Single subparser
    single_parser = subparsers.add_parser('single', help='Single gene analysis.')
    single_parser.add_argument('-i1', '--input1', type=str, required=True, help='Path to input directory containing datasets.')
    single_parser.add_argument('-i2', '--input2', type=str, required=True, help='Path to input directory containing DEG files.')
    single_parser.add_argument('-o', '--output', type=str, required=True, help='Path to output directory.')
    single_parser.add_argument('-wd', '--width', type=int, default=10, help='Width of output figures.')
    single_parser.add_argument('-ht', '--height', type=int, default=8, help='Height of output figures.')
    single_parser.add_argument('-d', '--dpi', type=int, default=300, help='DPI of output figures.')
    single_parser.add_argument('-ct', '--celltype', type=str, required=True, help='Title of cell type.')
    single_parser.add_argument('-g', '--gene', type=str, help='Gene name for analysis.')
    single_parser.add_argument('-srt', '--seurat', action='store_true', help='Whether the DEGs were calculated with Seurat.')
    single_parser.add_argument('-cmp', '--comparison', type=str, required=True, help='Groupwise comparison suffix (either VE, VC, or EC).')

    args = parser.parse_args()
    main(args)

