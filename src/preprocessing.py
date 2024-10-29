"""
Research Goal:
Use CCA to identify sets of genes that best differentiate between NSCLC patients and healthy controls, showing genetic factors in NSCLC development.
For this analysis, we use the TCGA-LUAD (Lung Adenocarcinoma) dataset, which is a type of NSCLC. 
This data is publicly available through the GDC Data Portal and UCSC Xena Browser.

"""

import os
import numpy as np
import pandas as pd
import logging

def log2_transform(df: pd.DataFrame) -> pd.DataFrame:
    """
    Applies log2 transformation to gene expression data.
    
    Arguments:
    - df: DataFrame where the first column is gene IDs and the rest are expression values
    
    Returns:
    - Transformed DataFrame with log2 applied to expression values.
    """
    # Exclude the first column
    gene_ids = df.iloc[:, 0]
    expression_data = df.iloc[:, 1:]
    
    log2_expression_data = np.log2(expression_data + 1) # avoid log(0)
    
    # Combine gene IDs and transformed data
    transformed_df = pd.concat([gene_ids, log2_expression_data], axis=1)
    return transformed_df

def preprocessing(gene_expression: pd.DataFrame, gene_pathway_mappings: pd.DataFrame) -> tuple:
    """  
    Preprocess gene expression data and calculate pathway activity scores.

    Parameters:
    ----------------
    gene_expression: pd.DataFrame
        Gene expression data.
    gene_pathway_mappings: pd.DataFrame
        Gene-pathway mappings.

    Returns:
    ----------------
    X_scaled: np.ndarray
        Preprocessed gene expression data.
    Y_scaled: np.ndarray
        Preprocessed pathway activity data.
    """
    # Check for missing values
    print("Missing values in gene expression data:")
    print(gene_expression.isnull().sum().sum())
    gene_expression.replace([np.inf, -np.inf], np.nan, inplace=True)
    gene_expression = gene_expression.dropna()

    # Apply log2 transformation
    gene_expression = log2_transform(gene_expression)

    # Calculate pathway activity scores
    def calculate_pathway_score(pathway_genes, expression_data):
        return expression_data.loc[pathway_genes].mean()

    pathway_scores = {}
    for pathway in gene_pathway_mappings.columns:
        pathway_genes = gene_pathway_mappings[gene_pathway_mappings[pathway].notna()].index
        pathway_scores[pathway] = gene_expression.apply(lambda col: calculate_pathway_score(pathway_genes, col))

    pathway_activity = pd.DataFrame(pathway_scores)

    # Ensure gene expression and pathway activity have the same samples
    common_samples = gene_expression.columns.intersection(pathway_activity.index)
    X = gene_expression[common_samples].T
    Y = pathway_activity.loc[common_samples]

    scaler = StandardScaler()
    X_scaled = pd.DataFrame(scaler.fit_transform(X), index=X.index, columns=X.columns)
    Y_scaled = pd.DataFrame(scaler.fit_transform(Y), index=Y.index, columns=Y.columns)

    print("Shape of X (gene expression):", X_scaled.shape)
    print("Shape of Y (pathway activity):", Y_scaled.shape)

    return X_scaled, Y_scaled

def classify(df: pd.DataFrame) -> pd.DataFrame:
    """
    Classify cell lines into subtypes.
    """
    classification = {
        'Luminal A': ['ZR751', 'ZR7530', 'ZR75B'],
        'Luminal B': ['BT474', 'CAMA1', 'HCC1428'],
        'HER2-enriched': ['AU565', 'UACC812', 'UACC893', 'MDAMB453'],
        'Basal-like': ['BT549', 'HCC1143', 'HCC1395'],
        'Normal-like': ['184A1', '184B5'],
        'Unclassified': ['21MT1', '21MT2', '21NT', '21PT', '600MPE', 'BT483', 'EFM192A', 'EFM192B', 'EFM192C', 'HCC1419']
    }
    def classify_cell_line(cell_line):
        for subtype, cell_lines in classification.items():
            if any(cl in cell_line for cl in cell_lines):
                return subtype
            return 'Unclassified'

    subtypes = {col: classify_cell_line(col.split('_')[1]) for col in df.columns}
    subtype_df = pd.DataFrame.from_dict(subtypes, orient='index', columns=['Subtype']) # new dataframe with cell line names as index

    subtype_df = ignore_subtypes(subtype_df)

    print(subtype_df)

    # Save the classifications
    subtype_df.to_csv('cell_line_subtypes.csv')

def ignore_subtypes(df_classified: pd.DataFrame) -> pd.DataFrame:
    """
    Ignore cell lines that have Unclassified subtype.
    """
    df_classified = df_classified.drop(columns=['Unclassified'])
    return df_classified

class Data_Processor:
    """
    Class to process the gene expression data.
    """
    def __init__(self):
        self.logger = self.setup_logger()
        self.logger.info("Data Processor initialized.")
        # self.X1, self.X2 = self.load_data()

    def setup_logger(self):
        """
        Set up the logger for the script.
        """
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.INFO)
        return logger

    def load_data(self) -> tuple:
        """
        Load the gene expression data.
        1. IlluminaHiSeq* (n=576) TCGA Hub, log2(x+1) transformed and RSEM normalized, not normalized across other cancer types, which is appropriate for distinguishing LUAD from normal samples within this cohort.
        2. Phenotypes (n=706) TCGA Hub -> which samples are tumor and which are normal + other clinical information.
    
        Returns:
        --------
        tuple
            Data matrices X1 and X2
        """
        first_filename = os.path.join("..", 'data', 'illumina', "HiSeqV2")
        second_filename = os.path.join("..", 'data', 'phenotypes', "TCGA.LUAD.sampleMap_LUAD_clinicalMatrix")

        try: 
            with open(first_filename, "r") as f:
                X1 = pd.read_csv(f, sep="\t", index_col=0)

            self.logger.info(f"Shape of Illumina data: {X1.shape}")
        except FileNotFoundError:
            raise FileNotFoundError(f"File {first_filename} not found.")
        
        try:
            with open(second_filename, "r") as f:
                X2 = pd.read_csv(f, sep="\t", index_col=0)

            self.logger.info(f"Shape of Phenotype data: {X2.shape}")         
        except FileNotFoundError:
            raise FileNotFoundError(f"File {second_filename} not found.")
        
        self.show_data(self)
        self.logger.info("Data loaded successfully.")

        return X1, X2
        
    def data_processing(self, X1: np.ndarray, X2: np.ndarray) -> tuple:
        """
        Process the gene expression data.
        
        Parameters:
        -----------
        X1 : np.ndarray
            Data matrix of size n x p
        X2 : np.ndarray
            Data matrix of size n x q
            
        Returns:
        --------
        tuple
            Processed data matrices X1 and X2
        """
        # Normalize gene expression data (log2 transformation) TODO: normalize all or inside sets?
        X1 = np.log2(X1 + 1)
        X2 = np.log2(X2 + 1)

        # Handle missing values
        X1 = np.nan_to_num(X1)
        X2 = np.nan_to_num(X2)
        
        return X1, X2
    
    def show_data(self):
        """
        Display the gene expression data.
        """
        print(f"Shape of Illumina data: {self.X1.shape}")
        print(f"Shape of Phenotype data: {self.X2.shape}")
        print("\n----------------------------------------------")
        print(f"First rows of gene expression profiles measured using RNA-Seq;")
        print(self.X1.head())
        print("\n----------------------------------------------")
        print(f"First rows of clinical data on phenotypes;")
        print(self.X2.head())

    