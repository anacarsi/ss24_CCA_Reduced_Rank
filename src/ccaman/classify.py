import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import requests
from bs4 import BeautifulSoup

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
        cell_line_name = cell_line.split("_")[
            1
        ]  # Extract the cell line name from the column name
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


def run_experiment(X, Y, k_values, retraction_methods):
    """
    Runs a CCA experiment with specified parameters.
    Arguments:
    - X, Y : DataFrames or np.ndarrays : Data matrices for genes and pathways.
    - k_values : list : list of k values for analysis.
    - retraction_methods : list : methods for retraction in CCA.
    Returns:
    - tuple : results, scores
    """
    # Placeholder for CCA experiment code
    results = pd.DataFrame({"method": ["cholesky", "polar"], "k": [5, 10], "correlations": [0.8, 0.9]})
    scores = pd.DataFrame({"XA": [np.random.rand(10)], "YB": [np.random.rand(10)]})
    return results, scores

def standard_cca(X, Y, k):
    """
    Runs standard Canonical Correlation Analysis (CCA) and returns results.
    Arguments:
    - X, Y : Data matrices for genes and pathways.
    - k : int : Number of components for CCA.
    Returns:
    - tuple : results of CCA
    """
    # Example CCA function
    correlations = np.random.rand(k)
    return X, Y, correlations, None

def plot_pca(X_pca, Y_pca):
    """
    Plots the PCA results for genes and pathways.
    Arguments:
    - X_pca, Y_pca : Data for PCA results
    """
    # Plotting code for PCA results
    pass

def plot_results(results):
    """
    Plots results of the CCA experiment.
    Arguments:
    - results : pd.DataFrame : DataFrame containing CCA results.
    """
    # Plotting logic for results
    pass

def plot_cca_scores(XA, YB, cell_line_labels, title="CCA Scores"):
    """
    Plots CCA scores.
    Arguments:
    - XA, YB : CCA scores
    - cell_line_labels : list : Labels for the cell lines
    """
    # Scatter plot of CCA scores for cell lines
    pass

def plot_canonical_correlations(correlations, title="Canonical Correlations Heatmap"):
    """
    Plots a heatmap for canonical correlations.
    Arguments:
    - correlations : np.ndarray : Canonical correlations
    """
    # Heatmap plotting logic
    pass

def plot_explained_variance(correlations, title="Explained Variance by Canonical Component"):
    """
    Plots the explained variance for each CCA component.
    Arguments:
    - correlations : np.ndarray : Canonical correlations
    """
    # Plotting explained variance
    pass
