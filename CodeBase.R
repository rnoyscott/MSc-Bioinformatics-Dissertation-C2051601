#Install the packages required for our analysis. Their usage is given where they are loaded.
#install.packages("Seurat")
#install.packages("qdapTools")
#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("cowplot")
#install.packages("tidyverse")
#install.packages("dplyr")
#install.packages("RColorBrewer")
#install.packages("SCpubr")
#BiocManager::install('biomaRt')
#BiocManager::install('clusterProfiler')
#BiocManager::install('org.Hs.eg.db')
#install.packages("GOSemSim")
#BiocManager::install('reticulate')
#remotes::install_github('satijalab/seurat-data')

#Load required packages and their dependencies.
#For loading and holding the data
library(Seurat)

#For adding patient IDs based on cell IDs
library(qdapTools)

#For plotting all graphs
library(ggplot2)
library(ggpubr)
library(cowplot)

#For functions
library(tidyverse)
library(dplyr)

#To give me pretty colours
library(RColorBrewer)

#For heatmaps
library(SCpubr)

#Enrichment Analysis
library(biomaRt)
library(clusterProfiler)

#Human genome reference
library(org.Hs.eg.db)
library(GOSemSim)

#For converting anndata back into RDS
library(reticulate)

#For curating data to produce the custom CellTypist model
library(SeuratData)

#Function for defining the top 10 GO pathways using a given list of DEGs
#Will also output all pathways, with corresponding description, gene ratio, p-value, adjusted p-value, gene IDs etc. as a CSV file
Pathways <- function(cellname, width, height, dpi, DEGfolder, ORAfolder, ont, path, seuratobject) {
  dir <- paste0(ORAfolder, cellname, "/")
  
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  
  # First, we define our background gene set, comprised of:
  # All protein-coding genes from our Seurat dataset subset (e.g. CD8 T cell subset)
  # Genes gathered from https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=C7 (GSEA immunological datasets)
  
  # List all files in the directory (this should be where the GSEA C7 files are stored)
  file_names <- list.files(path = path, pattern = ".tsv$")
  
  # Initialise a vector to store all genes from the GSEA C7 datasets
  C7_genes <- c()
  
  # Loop through the files and read the genes into a vector called C7_genes
  for (file_name in file_names) {
    # Read the data
    data <- readr::read_tsv(file.path(path, file_name), skip = 17, col_names = FALSE)
    
    # Extract genes
    genes <- strsplit(as.character(data[1, 2]), ",")[[1]]
    
    # Append to C7_genes
    C7_genes <- c(C7_genes, genes)
  }
  
  # Remove duplicate genes
  C7_genes <- unique(C7_genes)
  
  # Extract gene names from provided Seurat object and remove the duplicates
  exp_genes <- rownames(seuratobject@assays$RNA@counts)
  exp_genes <- sub("\\..*", "", exp_genes) # for handling multiple transcripts (e.g. GENE1.1, GENE1.2)
  exp_genes <- unique(exp_genes)
  
  # Connect to the Ensembl database
  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Query for gene names and biotypes based on unique identifiers
  gene_annotations <- biomaRt::getBM(
    attributes = c("entrezgene_id", "hgnc_symbol", "gene_biotype"),
    filters = "hgnc_symbol",
    values = exp_genes, 
    mart = ensembl)
  
  # Filter for protein-coding genes in the query results
  protein_coding_genes <- gene_annotations[gene_annotations$gene_biotype == "protein_coding", "hgnc_symbol"]
  
  # Subset original gene list extracted from the provided Seurat object (this will remove all non-protein-coding genes)
  exp_genes_protein_coding <- exp_genes[exp_genes %in% protein_coding_genes]
  
  # Merge the two lists (GSEA C7 & Seurat-associated genes) together and create an identifier column to state where they came from (either GSEA C7, or the Seurat object)
  BG_genes <- tibble(
    gene = c(C7_genes, exp_genes_protein_coding),
    source = c(rep("C7_genes", length(C7_genes)), rep("exp_genes", length(exp_genes_protein_coding)))
  )
  
  # Finally, remove all duplicated genes
  BG_genes <- BG_genes[!duplicated(BG_genes$gene), ]
  
  colnames(BG_genes)[1] <- "hgnc_symbol"
  
  # Identify corresponding ENTREZ ID
  gene_annotations <- biomaRt::getBM(
    attributes = c('entrezgene_id', 'hgnc_symbol'),
    filters = 'hgnc_symbol',
    values = BG_genes$hgnc_symbol,
    mart = ensembl)
  
  # Merge ENTREZ ID with gene name and output background to global environment, in case user wants to assess it
  BG_genes <- BG_genes |>
    dplyr::left_join(gene_annotations, by = 'hgnc_symbol')
  
  # Ensure that IDs are treated as character strings, important for the enrichment analysis
  BG_genes <- as.character(BG_genes$entrezgene_id)
  
  # The ensuing code will perform GSEA for a given list of DEGs, using the background gene set just defined
  
  # First letter in DEG file is upregulated when red
  # e.g. for VE, red means upregulated in Viraemic + downregulated in Exposed
  # e.g. for VE, blue means downregulated in Viraemic + upregulated in Exposed
  
  #Read in the DEG files  
  filelist <- list.files(path = DEGfolder, full.names = TRUE)
  for (DEGfile in filelist) {
    degs <- read.csv(DEGfile)
    
    # Standardise the header for genes
    colnames(degs)[1] <- "hgnc_symbol"
    
    # Identify corresponding ENTREZ ID
    gene_annotations <- biomaRt::getBM(
      attributes = c('entrezgene_id', 'hgnc_symbol'),
      filters = 'hgnc_symbol',
      values = degs$hgnc_symbol,
      mart = ensembl)
    
    # Merge ENTREZ ID with gene name
    degs <- degs |>
      dplyr::left_join(gene_annotations, by = 'hgnc_symbol')
    
    #Remove na from data set
    degs <- degs |>
      drop_na(entrezgene_id)
    degs <- degs |>
      drop_na(hgnc_symbol)
    
    # Filter for downregulated genes
    degs_downreg <- degs |>
      dplyr::filter(avg_log2FC < 0)
    
    # Filter significant upregulated genes based on adjusted p-value
    sigOE <- dplyr::filter(degs_downreg, p_val_adj < 0.05)
    sigOE_genes <- as.character(sigOE$entrezgene_id)
    
    # Check if there are upregulated genes
    if (!is.null(sigOE_genes)) {
      
      # Perform GO enrichment analysis on the significant genes
      GOen_down <- clusterProfiler::enrichGO(
        gene = sigOE_genes,        # Significant genes
        universe = BG_genes,       # Background genes
        keyType = 'ENTREZID',      # Annotation database
        OrgDb = org.Hs.eg.db,      # Reference species
        ont = ont,                 # Ontology (e.g., BP, MF, CC)
        pAdjustMethod = 'BH',      # P-value adjustment method
        qvalueCutoff = 0.05,       # Q-value cutoff
        readable = TRUE            # Convert gene IDs to gene names
      )
      
      # Convert enrichment results to a dataframe
      GOen_down_test <- data.frame(GOen_down)
      
      # Check if there are any enriched pathways
      if (nrow(GOen_down_test) > 0) {
        
        # Simplify the GO enrichment results
        GOen_down_simplified <- simplify(
          GOen_down, 
          cutoff = 0.7, 
          measure = "Wang", 
          semData = godata("org.Hs.eg.db", ont = ont)
        )
        
        # Convert the simplified results to a dataframe
        cluster_summary <- data.frame(GOen_down_simplified)
        
        # Assign condition comparison key based on the value of 'i'
        if (i == 1) {
          condition_suffix = '_EC'
        } else if (i == 2) {
          condition_suffix = '_VC'
        } else if (i == 3) {
          condition_suffix = '_VE'
        } else {
          stop('Error: Value of i is too large.')
        }
        
        condition_1 <- substr(condition_suffix, start=2, stop=2)
        condition_2 <- substr(condition_suffix, start=3, stop=3)
        
        # Write the simplified enrichment results to a CSV file
        write.csv(cluster_summary, paste0(dir, cellname, condition_suffix, "_ORA_downregulated_summary.csv"))
        
        # Create a dot plot of the top 10 pathways associated with upregulation
        downreg <- dotplot(GOen_down_simplified) +
          labs(title = paste0('Top 10 Pathways Associated with Downregulation in', condition_1, 'Compared to', condition_2)) +
          geom_text_repel(
            aes(label = Count),
            box.padding = unit(0.7, 'lines'),
            hjust = 0.30, 
            fontface = 'bold', 
            size = 4, 
            color = 'black'
          )
        
        upsetdown <- enrichplot::upsetplot(GOen_down_simplified, 10)
        upsetdown <- upsetdown +
          geom_bar(fill = 'grey', col = 'black') +
          theme_classic()
        
      } else {
        downreg <- NULL
        upsetdown <- NULL
        message('Unable to simplify pathways.')
      }
    } else {
      downreg <- NULL
      upsetdown <- NULL
      message('No upregulated genes to analyze.')
    }
    
    # Repeat for upregulated pathways
    degs_upreg <- degs |>
      dplyr::filter(avg_log2FC > 0)
    
    # Filter significant upregulated genes based on adjusted p-value
    sigOE <- dplyr::filter(degs_upreg, p_val_adj < 0.05)
    sigOE_genes <- as.character(sigOE$entrezgene_id)
    
    # Check if there are upregulated genes
    if (!is.null(sigOE_genes)) {
      
      # Perform GO enrichment analysis on the significant genes
      GOen_up <- clusterProfiler::enrichGO(
        gene = sigOE_genes,        # Significant genes
        universe = BG_genes,       # Background genes
        keyType = 'ENTREZID',      # Annotation database
        OrgDb = org.Hs.eg.db,      # Reference species
        ont = ont,                 # Ontology (e.g., BP, MF, CC)
        pAdjustMethod = 'BH',      # P-value adjustment method
        qvalueCutoff = 0.05,       # Q-value cutoff
        readable = TRUE            # Convert gene IDs to gene names
      )
      
      # Convert enrichment results to a dataframe
      GOen_up_test <- data.frame(GOen_up)
      
      # Check if there are any enriched pathways
      if (nrow(GOen_up_test) > 0) {
        
        # Simplify the GO enrichment results
        GOen_up_simplified <- simplify(
          GOen_up, 
          cutoff = 0.7, 
          measure = "Wang", 
          semData = godata("org.Hs.eg.db", ont = ont)
        )
        
        # Convert the simplified results to a dataframe
        cluster_summary <- data.frame(GOen_up_simplified)
        
        # Assign condition comparison key based on the value of 'i'
        if (i == 1) {
          condition_suffix = '_EC'
        } else if (i == 2) {
          condition_suffix = '_VC'
        } else if (i == 3) {
          condition_suffix = '_VE'
        } else {
          stop('Error: Value of i is too large.')
        }
        
        condition_1 <- substr(condition_suffix, start=2, stop=2)
        condition_2 <- substr(condition_suffix, start=3, stop=3)
        
        # Write the simplified enrichment results to a CSV file
        write.csv(cluster_summary, paste0(dir, cellname, condition_suffix, "_ORA_upregulated_summary.csv"))
        
        # Create a dot plot of the top 10 pathways associated with upregulation
        upreg <- dotplot(GOen_up_simplified) +
          labs(title = paste0('Top 10 Pathways Associated with Upregulation in', condition_1, 'Compared to', condition_2)) +
          geom_text_repel(
            aes(label = Count),
            box.padding = unit(0.7, 'lines'),
            hjust = 0.30, 
            fontface = 'bold', 
            size = 4, 
            color = 'black'
          )
        upsetup <- enrichplot::upsetplot(GOen_up_simplified, 10)
        upsetup <- upsetup +
          geom_bar(fill = 'grey', col = 'black') +
          theme_classic()
        
      } else {
        upreg <- NULL
        upsetup <- NULL
        message('Unable to simplify pathways.')
      }
    } else {
      upreg <- NULL
      upsetup <- NULL
      message('No upregulated genes to analyze.')
    }
    
    if (!is.null(upreg) && !is.null(downreg)) {
      # Plot upreg/downreg pathway results together in same window
      plots <- cowplot::plot_grid(upreg, downreg)
      
      # Save this as a file
      name <- paste(dir, cellname, condition_suffix, '_pathways_', dpi, 'dpi','.bmp', sep = '')
      bmp(filename = name,
          width = width, height = height, units = 'in', pointsize = 12, res = dpi)
      print(plots)
      dev.off()
    if (!is.null(upsetup) && !is.null(upsetdown)) {
      plots <- cowplot::plot_grid(NULL, upsetup, NULL, upsetdown, nrow = 2, ncol = 2)
      
      # Save this as a file
      name <- paste(dir, cellname, condition_suffix, '_pathways_upset_', dpi, 'dpi','.bmp', sep = '')
      bmp(filename = name,
          width = width, height = height, units = 'in', pointsize = 12, res = dpi)
      print(plots)
      dev.off()
    }
    }
    i = i + 1
  }
}

adata2RDS <- function(PrimaryData, MetaData) {
  library(reticulate)
  
  # Use reticulate to load the AnnData object
  anndata <- import("anndata")
  adata <- anndata$read_h5ad(PrimaryData)
  
  filtered_counts <- as.matrix(adata$layers[['filtered_counts_norm']])
  
  cell_names_df <- read_csv(file.path(MetaData, 'cell_names.csv'))
  gene_names_df <- read_csv(file.path(MetaData, 'gene_names.csv'))
  genes <- gene_names_df$gene_names
  cells <- cell_names_df$cell_names
  
  # Assign row and column names
  colnames(filtered_counts) <- genes
  rownames(filtered_counts) <- cells
  
  filtered_counts <- t(filtered_counts)
  
  # Convert AnnData to a Seurat object
  BloodCells <- Seurat::CreateSeuratObject(counts = filtered_counts, project = "P1")
  
  metadata_df <- as.data.frame(adata$obs)
  row.names(metadata_df) <- metadata_df$cell_ID
  
  # Add metadata and other information
  BloodCells <- AddMetaData(BloodCells, metadata = metadata_df)
  
  # Extract UMAP embeddings
  umap_embeddings <- adata$obsm[['X_umap']]
  umap_embeddings_df <- as.data.frame(umap_embeddings)
  rownames(umap_embeddings_df) <- cells
  
  
  # Add UMAP embeddings to Seurat object
  BloodCells[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = as.matrix(umap_embeddings_df),
    key = "UMAP_",
    assay = DefaultAssay(BloodCells)
  )
  
  # Export to global environment
  assign("BloodCells", BloodCells, envir = .GlobalEnv)
  
  rm(filtered_counts, cell_names_df, gene_names_df, cells, genes, metadata_df, umap_embeddings, umap_embeddings_df)
}

Subsetting <- function(SubsetCells, filename, groups, group1, group2) {
  # Subset the cells
  Subset <- subset(BloodCells, subset = majority_voting %in% SubsetCells)
  
  # Define base directory and create path
  base_dir <- '../Output/DGE/Seurat/'
  dir_path <- file.path(base_dir, filename)
  
  # Create the directory if it does not already exist
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Carry out DEG calculations
  DEG <- FindMarkers(Subset, logfc.threshold = 0, min.pct = 0, ident.1 = group1, ident.2 = group2, group.by = "AllSamples")
  
  # Format DEG results
  DEG_df <- data.frame(genename = row.names(DEG), DEG)
  
  # Create DEG filename
  file_path <- file.path(dir_path, paste0(groups, "_DEG_", filename, ".csv"))
  
  # Save as csv
  write.csv(DEG_df, file_path, row.names = FALSE)
}