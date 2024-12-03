# Main script for running the CCA analysis pipeline
# Ana Carsi 2024
import argparse
import os
import pandas as pd
import argparse
import os
import pandas as pd
from preprocessing import log2_transform, preprocessing
from ccaman import CCAMan 

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Breast Cancer Subtype Analysis Pipeline"
    )

    parser.add_argument(
        "--file_path", type=str, required=False, help="Path to the dataset"
    )

    parser.add_argument(
        "--log_transform",
        action="store_true",
        help="Apply log2 transformation to gene expression data",
    )

    parser.add_argument(
        "--run_cca",
        action="store_true",
        help="Run Canonical Correlation Analysis (CCA)",
    )

    return parser.parse_args()

def main():
    args = parse_args()

    # Initialize the CCAMan class and its pipeline
    ccaman = CCAMan(log_file="ccaman.log")
    ccaman.analyze()
    """
    # Log2 Transformation if specified
    if args.log_transform:
        log2_data = log2_transform(data)
        print("Log2 transformation applied.")
        print(log2_data.head())

    # Preprocessing for CCA
    file_path_combined = os.path.join(os.getcwd(), "..", "data", "combined_data.txt")
    filepath_pathway = os.path.join(
        os.getcwd(), "..", "data", "gene_pathway_mappings.csv"
    )
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
    """

if __name__ == "__main__":
    main()
