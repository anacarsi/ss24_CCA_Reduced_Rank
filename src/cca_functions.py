# Ana Carsi 2024
# Project: Review on Adaptive Canonical Correlation Analysis
# Description: This file contains the functions to perform Canonical Correlation Analysis (CCA) on the Stiefel Manifold.


import numpy as np
import matplotlib.pyplot as plt
from pymanopt import Problem
from pymanopt.manifolds import Stiefel, Product
from pymanopt.optimizers import TrustRegions
from pymanopt.function import Callable

def cca_on_stiefel(X: np.ndarray, Y: np.ndarray, k: int):
    """
    Perform CCA on the Stiefel Manifold.

    Parameters:
    -----------
    X : np.ndarray
        Data matrix of size n x p
    Y : np.ndarray
        Data matrix of size n x q
    k : int
        Number of canonical components to extract

    Returns:
    --------
    tuple
        Optimal points on the Stiefel manifold (U, V)
    """
    n, p = X.shape
    q = Y.shape[1]
    
    M1 = Stiefel(p, k)
    M2 = Stiefel(q, k)
    manifold = Product([M1, M2])
    
    def cost(point):
        U, V = point
        cost_value = -np.trace((U.T @ X.T @ Y @ V) @ (U.T @ X.T @ Y @ V).T)
        print(f"The cost function value is {cost_value}")
        return cost_value
    
    # Create a problem instance
    problem = Problem(manifold=manifold, cost=cost)
    optimizer = TrustRegions()
    
    Xopt = optimizer.run(problem).point
    
    return Xopt

def cca_with_tracking(X: np.ndarray, Y: np.ndarray, k: int, max_iter: int = 100):
    """
    Perform CCA on the Stiefel Manifold with tracking of objective values.

    Parameters:
    -----------
    X : np.ndarray
        Data matrix of size n x p
    Y : np.ndarray
        Data matrix of size n x q
    k : int
        Number of canonical components to extract
    max_iter : int, optional
        Maximum number of iterations (default is 100)

    Returns:
    --------
    list
        List of objective values over iterations
    """
    obj_values = []
    for _ in range(max_iter):
        Wx, Wy = cca_on_stiefel(X, Y, k)
        obj = -np.trace(np.dot(np.dot(Wx.T, X.T @ Y), Wy) @ np.dot(Wx.T, X.T @ Y), Wy).T
        obj_values.append(obj)
    
    return obj_values

def canonical_correlation_heatmap(X: np.ndarray, Y: np.ndarray, k: int):
    """
    Plot the canonical correlation heatmap.

    Parameters:
    -----------
    X : np.ndarray
        Data matrix of size n x p
    Y : np.ndarray
        Data matrix of size n x q
    k : int
        Number of canonical components to extract
    """
    Wx, Wy = cca_on_stiefel(X, Y, k)
    
    Ux = X @ Wx
    Uy = Y @ Wy
    
    corr_matrix = np.corrcoef(Ux.T, Uy.T)
    
    plt.figure(figsize=(10, 8))
    plt.imshow(corr_matrix, cmap='viridis')
    plt.colorbar()
    plt.xlabel
