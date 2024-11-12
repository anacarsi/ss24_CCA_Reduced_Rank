# Apply PCA to X and Y
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt

def apply_pca(X: np.ndarray, Y: np.ndarray, n_components: int) -> tuple:
    """
    Apply PCA to X and Y.

    Parameters:
    ----------
    X : np.ndarray
        The first dataset (genes with expressions).
    Y : np.ndarray
        The second dataset (pathways of the genes).
    n_components : int
        The number of components to compute.

    Returns:
    -------
    tuple : (X_pca, Y_pca)
        PCA-transformed X and Y.
    """
    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(X)
    Y_pca = pca.fit_transform(Y)
    return X_pca, Y_pca

def plot_pca(X: np.ndarray, Y: np.ndarray):
    """
    Plot PCA-transformed X and Y.

    Parameters:
    ----------
    X : np.ndarray
        PCA-transformed X.
    Y : np.ndarray
        PCA-transformed Y.
    title : str
        The title of the plot.
    """
    fig, ax = plt.subplots()
    ax.scatter(X[:, 0], X[:, 1], label='X', color='pink')
    ax.scatter(Y[:, 0], Y[:, 1], label='Y', color='lightpink')
    ax.set_title("PCA for gene classification on cell lines", fontsize=16)
    ax.legend()
    plt.show()