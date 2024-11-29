import argparse
import os
import pandas as pd
from utils.utils import load_data, classify_cancerous_celllines, classify_genes_in_pathways
from src.ccaman.preprocessing import log2_transform, preprocessing
from cca_functions import run_experiment, standard_cca
from visualize import plot_pca, plot_results #TODO: plot_cca_scores, plot_canonical_correlations, plot_explained_variance

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Breast Cancer Subtype Analysis Pipeline")

    parser.add_argument(
        '--file_path', 
        type=str, 
        required=True, 
        help='Path to the dataset'
    )
    
    parser.add_argument(
        '--log_transform', 
        action='store_true', 
        help='Apply log2 transformation to gene expression data'
    )

    parser.add_argument(
        '--run_cca', 
        action='store_true', 
        help='Run Canonical Correlation Analysis (CCA)'
    )

    parser.add_argument(
        '--plot_results', 
        action='store_true', 
        help='Plot CCA and other results'
    )

    return parser.parse_args()

def main():
    args = parse_args()

    # Load data
    data = load_data(args.file_path)
    print("Data loaded successfully.")
    print(data.head())

    # Log2 Transformation if specified
    if args.log_transform:
        log2_data = log2_transform(data)
        print("Log2 transformation applied.")
        print(log2_data.head())

    # Classify cancerous cell lines
    data = classify_cancerous_celllines()
    print("Cell lines classified.")
    print(data.head())

    # Preprocessing
    file_path_combined = os.path.join(os.getcwd(), "..", "data", "combined_data.txt")
    filepath_pathway = os.path.join(os.getcwd(), "..", "data", "gene_pathway_mappings.csv")
    gene_expression = pd.read_csv(file_path_combined, sep="\t", index_col=0)
    gene_pathway_mappings = pd.read_csv(filepath_pathway, index_col=0)

    X, Y = preprocessing(gene_expression, gene_pathway_mappings)
    n, p1 = X.shape
    _, p2 = Y.shape
    k = min(p1, p2, 10)

    # Run CCA if specified
    if args.run_cca:
        k_values = [5, 10, 15, 20]
        retraction_methods = ["cholesky", "polar"]
        results, scores = run_experiment(X, Y, k_values, retraction_methods)
        print("CCA experiment completed.")

    # Plot results if specified
    if args.plot_results:
        plot_results(results)
        print("Results plotted.")

if __name__ == "__main__":
    main()
