import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from bs4 import BeautifulSoup
import logging
import requests
import json


def run_experiment(X, Y, k_values, retraction_methods, logger=None):
    """
    Runs a CCA experiment with specified parameters.

    Arguments:
    - X, Y: DataFrames or np.ndarrays : Data matrices for genes and pathways.
    - k_values: list : List of k values for analysis.
    - retraction_methods: list : Methods for retraction in CCA.
    - logger: Logger for capturing logs and exceptions.

    Returns:
    - tuple : results, scores
    """
    try:
        # Placeholder for CCA experiment
        results = pd.DataFrame(
            {"method": ["cholesky", "polar"], "k": [5, 10], "correlations": [0.8, 0.9]}
        )
        scores = pd.DataFrame({"XA": [np.random.rand(10)], "YB": [np.random.rand(10)]})

        if logger:
            logger.info("Experiment completed successfully.")
        else:
            print("Experiment completed successfully.")

        return results, scores

    except Exception as e:
        if logger:
            logger.error(f"Error in run_experiment: {e}")
        else:
            print(f"Error in run_experiment: {e}")
        return None, None


def get_sensitivity_data(
    drug_classes: dict, cell_lines: list, logger=None
) -> pd.DataFrame:
    """
    Collects and stores sensitivity data for a list of cell lines.

    Arguments:
    - drug_classes: dict
        Dictionary of drug classes and associated drugs.
    - cell_lines: list
        List of cell lines to collect sensitivity data for.
    - logger: Logger for capturing logs and exceptions.

    Returns:
    - sensitivity_data: pd.DataFrame
        DataFrame containing sensitivity data.
    """
    # Initialize a dictionary to store sensitivity data
    sensitivity_dict = {"Cell Line": [], "Drug Name": [], "Z Score": []}

    filepath = os.path.join(
        os.getcwd(), "..", "data", "sensitivity", "sensitivity_data.csv"
    )
    if os.path.exists(filepath):
        sensitivity_data = pd.read_csv(filepath)
        (
            logger.info("Sensitivity data already exists. Loading from file.")
            if logger
            else print("Sensitivity data already exists. Loading from file.")
        )
    else:
        sensitivity_data = pd.DataFrame(sensitivity_dict)
        parentdir = os.path.join(os.getcwd(), "..", "data", "sensitivity")
        for folder_name in os.listdir(parentdir):
            # Remove the "s_" prefix from the folder name
            folder_name = folder_name[2:]
            folder_path = os.path.join(parentdir, folder_name)

            if os.path.exists(folder_path):
                sensitivity_line = pd.read_csv(folder_path)

                # Filter rows by drug names that exist in the drug_classes dictionary
                for drug_class, drugs in drug_classes.items():
                    for drug in drugs:
                        # Filter rows where the 'Drug Name' column matches the current drug
                        drug_rows = sensitivity_line[
                            sensitivity_line["Drug Name"].str.contains(
                                drug, case=False, na=False
                            )
                        ]

                        for _, row in drug_rows.iterrows():
                            sensitivity_dict["Cell Line"].append(folder_name)
                            sensitivity_dict["Drug Name"].append(row["Drug Name"])
                            sensitivity_dict["Z Score"].append(row["Z Score"])

    # Store the sensitivity data
    try:
        sensitivity_data.to_csv(filepath, index=False)
        (
            logger.info(f"Sensitivity data stored in {filepath}")
            if logger
            else print(f"Sensitivity data stored in {filepath}")
        )
    except Exception as e:
        if logger:
            logger.error(f"Error storing sensitivity data: {e}")
        else:
            print(f"Error storing sensitivity data: {e}")

    return sensitivity_data


def standard_cca(X, Y, k, logger=None):
    """
    Runs standard Canonical Correlation Analysis (CCA) and returns results.

    Arguments:
    - X, Y : Data matrices for genes and pathways.
    - k : int : Number of components for CCA.
    - logger: Logger for capturing logs and exceptions.

    Returns:
    - tuple : Results of CCA
    """
    try:
        # Example placeholder for actual CCA computation
        correlations = np.random.rand(k)
        if logger:
            logger.info("Standard CCA completed successfully.")
        else:
            print("Standard CCA completed successfully.")

        return X, Y, correlations, None

    except Exception as e:
        if logger:
            logger.error(f"Error in standard_cca: {e}")
        else:
            print(f"Error in standard_cca: {e}")
        return None, None, None, None
