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
import sys
from ..utils.utils import load_data
from classify import classify_cancerous_celllines, get_sensitivity_data
import json
import pandas as pd
import os


class CCAMan:
    def __init__(self, log_file="ccaman.log"):
        self.logger = logging.getLogger("CCAMan")
        self.logger.setLevel(logging.DEBUG)
        handler = logging.FileHandler(log_file)
        handler.setFormatter(
            logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        )
        self.logger.addHandler(handler)
        self.logger.info("CCAMan initialized")
        try:
            with open("config.json", "r") as f:
                self.config = json.load(f)
        except Exception as e:
            self.logger.error(f"Error loading config file: {e}")
            sys.exit(1)
        self.cell_line_names = []

    def load_data(self) -> pd.DataFrame:
        """
        Calls the global load_data function and ensures it logs through this instance's logger.
        """
        try:
            return load_data(self.config["file_path"], logger=self.logger)
        except Exception as e:
            self.logger.error(f"Error in load_data: {e}")
            raise e

    def analyze(self):
        # Load data
        self.data = self.load_data()
        self.logger.info("Data loaded successfully.")

        # Load the cell lines to extract sensitivity data
        filepath = os.join(
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
            "HER2 inhibitors": ["lapatinib"],
            "Hormone therapy": ["tamoxifen", "anastrozole", "letrozole", "exemestane"],
            "PARP inhibitors": ["olaparib"],
            "CDK4/6 inhibitors": ["palbociclib"],
            "PI3K inhibitors": ["alpelisib"],
        }
        self.sensitivity_data = get_sensitivity_data(
            drug_classes, self.cell_lines_names, self.logger
        )
        print(self.sensitivity_data.head())
