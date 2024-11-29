import logging
import numpy as np

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
        first_filename = os.path.join("..", "data", "illumina", "HiSeqV2")
        second_filename = os.path.join(
            "..", "data", "phenotypes", "TCGA.LUAD.sampleMap_LUAD_clinicalMatrix"
        )

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