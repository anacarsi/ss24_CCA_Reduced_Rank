# Ana Carsi 2024
import os
import pandas as pd
import gseapy as gp
import mygene
import csv
import numpy as np
import requests
import time
import logging


"""
Dataset overview:
    56 breast cancer cell lines were profiled
    The data represents gene expression levels in these cell lines
    Each cell line is in an unperturbed, baseline state
File structure:
    Each folder corresponds to a different breast cancer cell line
    Inside each folder is a .txt file containing gene expression data for that specific cell line
"""

import logging
import os
import pandas as pd

class CCAMan:
    def __init__(self, log_file="ccaman.log"):
        # Set up the logger for the class
        self.logger = logging.getLogger("CCAMan")
        self.logger.setLevel(logging.DEBUG)
        handler = logging.FileHandler(log_file)
        handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
        self.logger.addHandler(handler)

    def load_data(self, path_input: str) -> pd.DataFrame:
        """
        Calls the global load_data function and ensures it logs through this instance's logger.
        """
        try:
            return load_data(path_input, logger=self.logger)
        except Exception as e:
            self.logger.error(f"Error in load_data: {e}")
            raise  # Optionally, re-raise the exception if it needs to be handled higher up

# External function modified to accept a logger
def load_data(path_input: str, logger=None) -> pd.DataFrame:
    """
    Function to combine the cell line files into a single DataFrame with exception handling.

    Arguments:
    - path_input: The root directory containing folders with .txt files.
    - logger: Optional logger instance for logging exceptions.

    Returns:
    - combined_data: The combined DataFrame, or an empty DataFrame if an error occurs.
    """
    try:
        if os.path.exists("combined_data.txt"):
            logger.info("Combined data already exists. Loading from file.") if logger else print("Combined data already exists. Loading from file.")
            return pd.read_csv("combined_data.txt", sep="\t", index_col=0)
        else:
            combined_data = pd.DataFrame()
            parent_dir = os.path.abspath(os.path.join(os.getcwd(), ".."))
            output_file = os.path.join(parent_dir, "data", "combined_data.txt")
            logger.info(f"Saving combined data to {output_file}") if logger else print(f"Saving combined data to {output_file}")

            for folder_name in os.listdir(path_input):
                folder_path = os.path.join(path_input, folder_name)

                if os.path.isdir(folder_path):
                    txt_files = [f for f in os.listdir(folder_path) if f.endswith(".txt")]

                    if txt_files:
                        file_path = os.path.join(folder_path, txt_files[0])
                        try:
                            df = pd.read_csv(file_path, sep="\t", header=0)

                            if combined_data.empty:
                                combined_data = pd.DataFrame(index=df.iloc[:, 0])
                                logger.info(f"Index set to {df.columns[0]}") if logger else print(f"Index set to {df.columns[0]}")

                            column_name = txt_files[0].replace(".txt", "")

                            if len(df) == len(combined_data):
                                combined_data[column_name] = df.iloc[:, 1].values
                            else:
                                logger.warning(f"Row mismatch in {folder_name}. Skipping this file.") if logger else print(f"Row mismatch in {folder_name}. Skipping this file.")
                        except Exception as e:
                            logger.error(f"Error processing file {file_path}: {e}") if logger else print(f"Error processing file {file_path}: {e}")

            combined_data.to_csv(output_file, sep="\t")
            logger.info(f"Combined data has been saved to {output_file}") if logger else print(f"Combined data has been saved to {output_file}")
            return combined_data
    except Exception as e:
        if logger:
            logger.error(f"Error in load_data: {e}")
        else:
            print(f"Error in load_data: {e}")
        return pd.DataFrame()  # Return an empty DataFrame in case of an error


def log_stability_data(filename: str, iteration: int, A: np.ndarray, B: np.ndarray, G_A: np.ndarray, G_B: np.ndarray, logger=None):
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
            logger.info(f"Stability data logged at iteration {iteration}") if logger else print(f"Stability data logged at iteration {iteration}")
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
            logger.info(f"Stability log file {filename} initialized with headers") if logger else print(f"Stability log file {filename} initialized with headers")
    except Exception as e:
        if logger:
            logger.error(f"Error initializing stability log file {filename}: {e}")
        else:
            print(f"Error initializing stability log file {filename}: {e}")
