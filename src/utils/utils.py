import os
import pandas as pd

"""
Dataset overview:
    56 breast cancer cell lines were profiled
    The data represents gene expression levels in these cell lines
    Each cell line is in an unperturbed, baseline state
File structure:
    Each folder corresponds to a different breast cancer cell line
    Inside each folder is a .txt file containing gene expression data for that specific cell line
"""


def load_data(path_input: str, path_output: str) -> pd.DataFrame:
    """
    Function to combine the cell lines files into a single DataFrame.
    """
    combined_data = pd.DataFrame()
    output_file = os.path.join(path_output, 'combined_data.txt')

    # Iterate through all folders in the root directory
    for folder_name in os.listdir(path_input):
        folder_path = os.path.join(path_input, folder_name)

        if os.path.isdir(folder_path):
            # Look for a .txt file in the folder
            txt_files = [f for f in os.listdir(folder_path) if f.endswith('.txt')]
            
            if txt_files:
                file_path = os.path.join(folder_path, txt_files[0])
                df = pd.read_csv(file_path, sep='\t', header=0)
                
                # If this is the first file, use its first column as the index
                if combined_data.empty:
                    combined_data = pd.DataFrame(index=df.iloc[:, 0])
                
                # Add the second column to the combined data
                combined_data[folder_name] = df.iloc[:, 1].values

    # Save the combined data to a file
    combined_data.to_csv(output_file, sep='\t')

    print(f"Combined data has been saved to {output_file}")
