"""
(GSE48213): Identifying gene expression patterns associated with 
different breast cancer subtypes
In file: 
1. Column 1 (EnsEMBL_Gene_ID): 
    unique identifier for each gene from the Ensembl database.
2. Column 2 (e.g., MDAMB453): 
    expression value for each gene in the specific cell line.

These are normalized read counts or FPKM/TPM values (Fragments/Transcripts Per Kilobase Million).
Higher values indicate higher expression of the gene in that cell line, zero values indicate that the gene is not expressed (or expression is below detection threshold)
"""

import logging
from .utils.utils import combine_data
from .classify import get_sensitivity_data
import json
import pandas as pd
import os
from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import StandardScaler
import numpy as np
import matplotlib.pyplot as plt


class CCAMan:
    def __init__(self, log_file="ccaman.log"):
        if os.path.exists(log_file):
            os.remove(log_file)
        self.logger = logging.getLogger("CCAMan")
        self.logger.setLevel(logging.DEBUG)
        handler = logging.FileHandler(log_file)
        handler.setFormatter(
            logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        )
        self.logger.addHandler(handler)
        self.logger.info("CCAMan initialized")
        self.cell_line_names = []
        self.genes_to_cellline = self.load_data()
        self.logger.info("Data loaded successfully.")

    def load_data(self) -> pd.DataFrame:
        """
        Calls the global load_data function and ensures it logs through this instance's logger.
        """
        try:
            return combine_data(logger=self.logger)
        except Exception as e:
            self.logger.error(f"Error in load_data: {e}")
            raise e

    def analyze(self):
        # Load the cell lines to extract sensitivity data
        filepath = os.path.join(
            os.getcwd(), "..", "data", "cell_lines", "cell_line_names_filtered.json"
        )
        try:
            with open(filepath, "w") as json_file:
                self.cell_line_names = json.load(json_file)
        except Exception as e:
            if self.logger:
                self.logger.error(f"Error loading filtered cell line names: {e}")
            else:
                print(f"Error loading filtered cell line names: {e}")

        # Extract the sensitivity data
        drug_classes = {
            "HER2 inhibitors": ["Lapatinib"],
            "Hormone therapy": ["Tamoxifen"],
            "PARP inhibitors": ["Olaparib"],
            "CDK4/6 inhibitors": ["Palbociclib"],
            "PI3K inhibitors": ["Alpelisib"],
        }
        self.sensitivity_data = get_sensitivity_data(
            drug_classes, self.cell_line_names, self.logger
        )
