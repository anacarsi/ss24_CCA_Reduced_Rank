# Ana Carsi 2024
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
    Function to combine the cell line files into a single DataFrame.
    
    Arguments:
    - path_input: The root directory containing folders with .txt files.
    - path_output: The directory where the combined data will be saved.
    
    Returns:
    - combined_data: The combined DataFrame.
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
                # Load the data from the .txt file (taking the first .txt file found)
                file_path = os.path.join(folder_path, txt_files[0])
                df = pd.read_csv(file_path, sep='\t', header=0)

                # If the combined_data is empty, initialize it using the first column as the index
                if combined_data.empty:
                    combined_data = pd.DataFrame(index=df.iloc[:, 0])
                    print(f"Index set to {df.columns[0]}")
                
                # Strip the .txt extension from the folder name for the column name
                column_name = txt_files[0].replace('.txt', '')
                
                # Add the second column from the current dataframe to combined_data
                if len(df) == len(combined_data):
                    combined_data[column_name] = df.iloc[:, 1].values
                else:
                    print(f"Row mismatch in {folder_name}. Skipping this file.")

    # Save the combined data to a file
    combined_data.to_csv(output_file, sep='\t')

    print(f"Combined data has been saved to {output_file}")
    return combined_data
