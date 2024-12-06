# Ana Carsi 2024
import os
import pandas as pd
import csv
import numpy as np


"""
Dataset overview:
    56 breast cancer cell lines were profiled
    The data represents gene expression levels in these cell lines
    Each cell line is in an unperturbed, baseline state
File structure:
    Each folder corresponds to a different breast cancer cell line
    Inside each folder is a .txt file containing gene expression data for that specific cell line
"""


# External function modified to accept a logger
def combine_data(logger=None) -> pd.DataFrame:
    """
    Function to combine the cell line files into a single DataFrame with exception handling.

    Arguments:
    - logger: Optional logger instance for logging exceptions.

    Returns:
    - genes_to_cellline: The combined DataFrame, or an empty DataFrame if an error occurs.
    """
    try:
        filepath = os.path.join(
            os.getcwd(), "..", "data", "GSE48213", "genes_to_cellline.csv"
        )
        if os.path.exists(filepath):
            (
                logger.info("Combined data already exists. Loading from file.")
                if logger
                else print("Combined data already exists. Loading from file.")
            )
            return pd.read_csv(filepath, sep="\t", index_col=0)
        else:
            print("Combining data 2")
            genes_to_cellline = pd.DataFrame()
            (
                logger.info(f"Saving combined data to {filepath}")
                if logger
                else print(f"Saving combined data to {filepath}")
            )
            parentdir = os.path.join(filepath, "..")
            for folder_name in os.listdir(parentdir):
                folder_path = os.path.join(parentdir, folder_name)
                if os.path.isdir(folder_path):
                    txt_files = [
                        f for f in os.listdir(folder_path) if f.endswith(".txt")
                    ]
                    if txt_files:
                        file_path = os.path.join(folder_path, txt_files[0])
                        try:
                            df = pd.read_csv(file_path, sep="\t", header=0)

                            if genes_to_cellline.empty:
                                genes_to_cellline = pd.DataFrame(index=df.iloc[:, 0])
                                (
                                    logger.info(f"Index set to {df.columns[0]}")
                                    if logger
                                    else print(f"Index set to {df.columns[0]}")
                                )

                            column_name = txt_files[0].replace(".txt", "")
                            column_name = column_name.split("_")[1]
                            if len(df) == len(genes_to_cellline):
                                genes_to_cellline[column_name] = df.iloc[:, 1].values
                            else:
                                (
                                    logger.warning(
                                        f"Row mismatch in {folder_name}. Skipping this file."
                                    )
                                    if logger
                                    else print(
                                        f"Row mismatch in {folder_name}. Skipping this file."
                                    )
                                )
                        except Exception as e:
                            (
                                logger.error(f"Error processing file {file_path}: {e}")
                                if logger
                                else print(f"Error processing file {file_path}: {e}")
                            )

            genes_to_cellline.to_csv(filepath, sep="\t")
            (
                logger.info(f"Combined data has been saved to {filepath}")
                if logger
                else print(f"Combined data has been saved to {filepath}")
            )

            print("Wrote to file")
            return genes_to_cellline
    except Exception as e:
        if logger:
            logger.error(f"Error in load_data: {e}")
        else:
            print(f"Error in load_data: {e}")
        return pd.DataFrame()  # Return an empty DataFrame in case of an error


def log_stability_data(
    filename: str,
    iteration: int,
    A: np.ndarray,
    B: np.ndarray,
    G_A: np.ndarray,
    G_B: np.ndarray,
    logger=None,
):
    """
    Logs the L1 and Linf norms of matrices A, B, G_A, and G_B for each iteration, with exception handling.
    """
    try:
        with open(filename, mode="a", newline="") as file:
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
            writer.writerow(
                [iteration]
                + [A_l1, A_inf]
                + [B_l1, B_inf]
                + [G_A_l1, G_A_inf]
                + [G_B_l1, G_B_inf]
            )
            (
                logger.info(f"Stability data logged at iteration {iteration}")
                if logger
                else print(f"Stability data logged at iteration {iteration}")
            )
    except Exception as e:
        if logger:
            logger.error(f"Error logging stability data at iteration {iteration}: {e}")
        else:
            print(f"Error logging stability data at iteration {iteration}: {e}")


def init_stability_log(filename: str, k: int, logger=None):
    """
    Initializes the stability log file with headers, with exception handling.
    """
    try:
        headers = (
            ["Iteration"]
            + [f"A_l1", f"A_inf"]
            + [f"B_l1", f"B_inf"]
            + [f"G_A_l1", f"G_A_inf"]
            + [f"G_B_l1", f"G_B_inf"]
        )
        with open(filename, mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(headers)
            (
                logger.info(f"Stability log file {filename} initialized with headers")
                if logger
                else print(f"Stability log file {filename} initialized with headers")
            )
    except Exception as e:
        if logger:
            logger.error(f"Error initializing stability log file {filename}: {e}")
        else:
            print(f"Error initializing stability log file {filename}: {e}")


def get_unique_drugs(sensitivity_data: pd.DataFrame) -> list:
    """
    Extracts unique drug names from the sensitivity data DataFrame.

    Arguments:
    - sensitivity_data: pd.DataFrame : DataFrame containing sensitivity data.

    Returns:
    - list : List of unique drug names.
    """
    return sensitivity_data["Drug Name"].unique().tolist()
