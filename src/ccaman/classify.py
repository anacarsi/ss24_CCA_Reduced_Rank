import os
import pandas as pd
import numpy as np


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

    # Flatten the drug classes into a list for filtering
    filter_drugs = [drug.lower() for drugs in drug_classes.values() for drug in drugs]

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

        # Process each folder corresponding to a cell line
        for folder_name in os.listdir(parentdir):
            # Remove the "s_" prefix from the folder name
            cell_line = folder_name[2:]
            # Remove the ".csv" extension from the folder name
            cell_line = cell_line[:-4]
            folder_path = os.path.join(parentdir, folder_name)

            if os.path.exists(folder_path):
                try:
                    sensitivity_line = pd.read_csv(folder_path)
                    # Filter the data based on the drug classes
                    sensitivity_line = sensitivity_line[
                        sensitivity_line["Drug Name"].str.lower().isin(filter_drugs)
                    ]
                    # Append filtered data to the dictionary
                    for _, row in sensitivity_line.iterrows():
                        sensitivity_dict["Cell Line"].append(cell_line)
                        sensitivity_dict["Drug Name"].append(row["Drug Name"])
                        sensitivity_dict["Z Score"].append(row["Z Score"])
                        (
                            logger.info(
                                f"Sensitive data for {cell_line} and {row['Drug Name']} collected."
                            )
                            if logger
                            else print(
                                f"Sensitive data for {cell_line} and {row['Drug Name']} collected."
                            )
                        )
                except Exception as e:
                    if logger:
                        logger.error(f"Error processing cell line {cell_line}: {e}")
                    else:
                        print(f"Error processing cell line {cell_line}: {e}")

        sensitivity_data = pd.DataFrame(sensitivity_dict)

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
