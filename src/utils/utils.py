# Ana Carsi 2024
import os
import pandas as pd
import gseapy as gp
import mygene
import csv
import numpy as np
import requests
import time


"""
Dataset overview:
    56 breast cancer cell lines were profiled
    The data represents gene expression levels in these cell lines
    Each cell line is in an unperturbed, baseline state
File structure:
    Each folder corresponds to a different breast cancer cell line
    Inside each folder is a .txt file containing gene expression data for that specific cell line
"""

def load_data(path_input: str) -> pd.DataFrame:
    """
    Function to combine the cell line files into a single DataFrame.
    
    Arguments:
    - path_input: The root directory containing folders with .txt files.
    - path_output: The directory where the combined data will be saved.
    
    Returns:
    - combined_data: The combined DataFrame.
    """
    if os.path.exists('combined_data.txt'):
        print("Combined data already exists. Loading from file.")
        return pd.read_csv('combined_data.txt', sep='\t', index_col=0)
    else:
        combined_data = pd.DataFrame()
        parent_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
        output_file = os.path.join(parent_dir, 'data', 'combined_data.txt')
        print(f"Saving combined data to {output_file}")

        # Iterate through all folders in the root directory
        for folder_name in os.listdir(path_input):
            folder_path = os.path.join(path_input, folder_name)

            if os.path.isdir(folder_path):
                # Look for a .txt file in the folder
                txt_files = [f for f in os.listdir(folder_path) if f.endswith('.txt')]
                
                if txt_files:
                    # Load the data from the .txt file (taking the first .txt file found)
                    file_path = os.path.join(folder_path, txt_files[0])
                    df = pd.read_csv(file_path, sep='\t', header=0)

                    # If the combined_data is empty, initialize it using the first column as the index
                    if combined_data.empty:
                        combined_data = pd.DataFrame(index=df.iloc[:, 0])
                        print(f"Index set to {df.columns[0]}")
                    
                    # Strip the .txt extension from the folder name for the column name
                    column_name = txt_files[0].replace('.txt', '')
                    
                    # Add the second column from the current dataframe to combined_data
                    if len(df) == len(combined_data):
                        combined_data[column_name] = df.iloc[:, 1].values
                    else:
                        print(f"Row mismatch in {folder_name}. Skipping this file.")

        # Save the combined data to a file
        combined_data.to_csv(output_file, sep='\t')

        print(f"Combined data has been saved to {output_file}")
        return combined_data

def classify_genes_in_pathways():
    """ 
    Find which genes are in which KEGG pathways.
    """
    data = pd.read_csv('combined_data.txt', sep='\t', index_col=0)
    genes = data.index.tolist()

    # Convert Ensembl IDs to gene symbols
    mg = mygene.MyGeneInfo()
    gene_info = mg.querymany(genes, scopes='ensembl.gene', fields='symbol', species='human')
    ensembl_to_symbol = {gene['query']: gene.get('symbol', '') for gene in gene_info if 'symbol' in gene} # create a dictionary of Ensembl IDs to gene symbols

    # Load KEGG pathway gene sets
    kegg_gene_sets = gp.get_library('KEGG_2021_Human')
    gene_pathway_map = {}

    for pathway, pathway_genes in kegg_gene_sets.items():
        # Find genes that are both in data and in the pathway
        common_genes = set(ensembl_to_symbol.values()).intersection(pathway_genes)
        
        # If there are common genes, add them to the gene_pathway_map
        for gene_symbol in common_genes:
            ensembl_ids = [ensembl for ensembl, symbol in ensembl_to_symbol.items() if symbol == gene_symbol]
            for ensembl_id in ensembl_ids:
                if ensembl_id not in gene_pathway_map:
                    gene_pathway_map[ensembl_id] = []
                gene_pathway_map[ensembl_id].append(pathway)

    # Create a DataFrame from the gene_pathway_map
    result_df = pd.DataFrame.from_dict(gene_pathway_map, orient='index')
    result_df.columns = [f'Pathway_{i+1}' for i in range(result_df.shape[1])]

    result_df.to_csv('gene_pathway_mappings.csv')
    
    total_genes = len(genes)
    mapped_genes = len(gene_pathway_map)
    print(f"Total genes in your data: {total_genes}")
    print(f"Genes mapped to KEGG pathways: {mapped_genes}")
    print(f"Percentage of genes mapped: {mapped_genes/total_genes*100:.2f}%")

    # Print the first few rows of the result
    print("\nFirst few gene-pathway mappings:")
    print(result_df.head())

def log_stability_data(filename, iteration, A, B, G_A, G_B):
    """
    Logs the L1 and Linf norms of matrices A, B, G_A, and G_B for each iteration.
    """
    with open(filename, mode='a', newline='') as file:
        writer = csv.writer(file)
        
        # Compute L1 and L∞ norms for each matrix
        A_l1 = np.linalg.norm(A, ord=1)
        A_inf = np.linalg.norm(A, ord=np.inf)
        B_l1 = np.linalg.norm(B, ord=1)
        B_inf = np.linalg.norm(B, ord=np.inf)
        G_A_l1 = np.linalg.norm(G_A, ord=1)
        G_A_inf = np.linalg.norm(G_A, ord=np.inf)
        G_B_l1 = np.linalg.norm(G_B, ord=1)
        G_B_inf = np.linalg.norm(G_B, ord=np.inf)
        
        # Write iteration and norms to the CSV
        writer.writerow([iteration] + 
                        [A_l1, A_inf] + 
                        [B_l1, B_inf] + 
                        [G_A_l1, G_A_inf] + 
                        [G_B_l1, G_B_inf])

def init_stability_log(filename: str, k: int):
    """
    Initializes the stability log file with headers.
    """
    headers = ['Iteration'] + \
              [f"A_l1", f"A_inf"] + \
              [f"B_l1", f"B_inf"] + \
              [f"G_A_l1", f"G_A_inf"] + \
              [f"G_B_l1", f"G_B_inf"]
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(headers)

def classify_cancerous_celllines() -> pd.DataFrame:
    """
    Create a new dataset with the classification between immortalized, cancerous or unclassified in breast cancer tumours.

    Returns:
    - classification_df: A DataFrame with the classification of each cell line.
    """
    # The first column contains the cell line names (i.e., GSM1172844_184A1)
    data = pd.read_csv("combined_data.txt", sep="\t", index_col=0)
    cell_lines = data.columns.tolist()
    cell_line_names = []

    # Loop through the cell lines and classify each one
    for cell_line in cell_lines:
        cell_line_name = cell_line.split("_")[1]  # Extract the cell line name from the column name
        print(f"Classifying {cell_line_name}...")
        cell_line_names.append(cell_line_name)
    
    store_cell_lines(pd.DataFrame(cell_line_names))
        
    print("Classification complete.")
    return pd.DataFrame(cell_line_names)

def store_cell_lines(cell_line_names: pd.DataFrame) -> None:
    """
    Store the cell lines in a CSV file.

    Arguments:
    - cell_line_names: The DataFrame containing the cell line names.
    """
    file_path = os.path.join(os.getcwd(), "..", "data", "cell_lines")
    if not os.path.exists(file_path):
        os.makedirs(file_path)
    try: 
        with open(os.path.join(file_path, "cell_line_names.csv"), "w") as f:
            cell_line_names.to_csv(f, index=False)
    except Exception as e:
        print(f"Error storing cell lines: {e}")

    print(f"Cell lines stored in {file_path}")
    