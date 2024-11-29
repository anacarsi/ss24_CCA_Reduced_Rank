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
from sklearn.preprocessing import StandardScaler


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

    log2_expression_data = np.log2(expression_data + 1)  # avoid log(0)

    # Combine gene IDs and transformed data
    transformed_df = pd.concat([gene_ids, log2_expression_data], axis=1)
    return transformed_df


def preprocessing(
    gene_expression: pd.DataFrame, gene_pathway_mappings: pd.DataFrame
) -> tuple:
    """
    Preprocess gene expression data and calculate pathway activity scores.

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
        pathway_genes = gene_pathway_mappings[
            gene_pathway_mappings[pathway].notna()
        ].index
        pathway_scores[pathway] = gene_expression.apply(
            lambda col: calculate_pathway_score(pathway_genes, col)
        )

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
        "Luminal A": ["ZR751", "ZR7530", "ZR75B"],
        "Luminal B": ["BT474", "CAMA1", "HCC1428"],
        "HER2-enriched": ["AU565", "UACC812", "UACC893", "MDAMB453"],
        "Basal-like": ["BT549", "HCC1143", "HCC1395"],
        "Normal-like": ["184A1", "184B5"],
        "Unclassified": [
            "21MT1",
            "21MT2",
            "21NT",
            "21PT",
            "600MPE",
            "BT483",
            "EFM192A",
            "EFM192B",
            "EFM192C",
            "HCC1419",
        ],
    }

    def classify_cell_line(cell_line):
        for subtype, cell_lines in classification.items():
            if any(cl in cell_line for cl in cell_lines):
                return subtype
            return "Unclassified"

    subtypes = {col: classify_cell_line(col.split("_")[1]) for col in df.columns}
    subtype_df = pd.DataFrame.from_dict(
        subtypes, orient="index", columns=["Subtype"]
    )  # new dataframe with cell line names as index

    subtype_df = ignore_subtypes(subtype_df)

    print(subtype_df)

    # Save the classifications
    subtype_df.to_csv("cell_line_subtypes.csv")


def ignore_subtypes(df_classified: pd.DataFrame) -> pd.DataFrame:
    """
    Ignore cell lines that have Unclassified subtype.
    """
    df_classified = df_classified.drop(columns=["Unclassified"])
    return df_classified

