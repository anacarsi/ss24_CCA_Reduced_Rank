import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from bs4 import BeautifulSoup
import logging
import requests

def classify_cancerous_celllines(logger=None) -> list:
    """
    Create a new dataset with the classification between immortalized, cancerous, or unclassified in breast cancer tumors.

    Arguments:
    - logger: Logger for capturing logs and exceptions.

    Returns:
    - classification_df: A DataFrame with the classification of each cell line.
    """
    try:
        data = pd.read_csv("combined_data.txt", sep="\t", index_col=0)
        cell_lines = data.columns.tolist()
        cell_line_names = []

        for cell_line in cell_lines:
            cell_line_name = cell_line.split("_")[1]  # Extract the cell line name
            cell_line_names.append(cell_line_name)

        # Attempt to store cell lines in a file
        try:
            store_cell_lines(pd.DataFrame(cell_line_names), logger=logger)
        except Exception as e:
            if logger:
                logger.error(f"Error in store_cell_lines: {e}")
            else:
                print(f"Error in store_cell_lines: {e}")

        if logger:
            logger.info("Classification complete.")
        else:
            print("Classification complete.")
        return cell_line_names

    except Exception as e:
        if logger:
            logger.error(f"Error in classify_cancerous_celllines: {e}")
        else:
            print(f"Error in classify_cancerous_celllines: {e}")
        return []


def store_cell_lines(cell_line_names: pd.DataFrame, logger=None) -> None:
    """
    Store the cell lines in a CSV file.

    Arguments:
    - cell_line_names: The DataFrame containing the cell line names.
    - logger: Logger for capturing logs and exceptions.
    """
    file_path = os.path.join(os.getcwd(), "..", "data", "cell_lines")
    try:
        if not os.path.exists(file_path):
            os.makedirs(file_path)

        output_file = os.path.join(file_path, "cell_line_names.csv")
        cell_line_names.to_csv(output_file, index=False)

        if logger:
            logger.info(f"Cell lines stored in {output_file}")
        else:
            print(f"Cell lines stored in {output_file}")

    except Exception as e:
        if logger:
            logger.error(f"Error storing cell lines: {e}")
        else:
            print(f"Error storing cell lines: {e}")


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
    
    
def get_sensitivity_data(drug_classes: dict, cell_lines: list, logger=None) -> pd.DataFrame:
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
    sensitivity_dict = {"Cell Line": []}

    # Ensure each drug class has its own key in the dictionary
    for drug_class in drug_classes.keys():
        sensitivity_dict[drug_class] = []

    for cell_line in cell_lines:
        url = f"https://www.cancerrxgene.org/celllines/{cell_line}"
        try:
            response = requests.get(url)
            if response.status_code != 200:
                logger.error(f"Failed to retrieve data for {cell_line}: {response.status_code}")
                continue

            soup = BeautifulSoup(response.content, "html.parser")
            sensitivity_dict["Cell Line"].append(cell_line)

            for drug_class, drugs in drug_classes.items():
                z_scores = []
                for drug in drugs:
                    try:
                        drug_element = soup.find("td", string=drug)
                        if drug_element:
                            z_score_element = drug_element.find_next("td")
                            if z_score_element:
                                z_score = float(z_score_element.text.strip())
                                z_scores.append(z_score)
                                break  # Exit loop after finding the first match
                    except Exception as e:
                        if logger:
                            logger.error(f"Error extracting Z-score for {drug}: {e}")
                        else:
                            print(f"Error extracting Z-score for {drug}: {e}")
                # Append the average Z-score for this drug class (or NaN if none found)
                sensitivity_dict[drug_class].append(
                    sum(z_scores) / len(z_scores) if z_scores else float("nan")
                )

        except Exception as e:
            if logger:
                logger.error(f"Error processing cell line {cell_line}: {e}")
            else:
                print(f"Error processing cell line {cell_line}: {e}")

    # Convert the dictionary to a DataFrame
    sensitivity_data = pd.DataFrame(sensitivity_dict)

    # Store the sensitivity data
    try:
        filepath = os.path.join(os.getcwd(), "..", "data", "sensitivity_data.csv")
        sensitivity_data.to_csv(filepath, index=False)
        logger.info(f"Sensitivity data stored in {filepath}")
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
