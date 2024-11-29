# Ana Carsi 2024
# Project: Review on Adaptive Canonical Correlation Analysis
# Description: This file contains the functions to perform Canonical Correlation Analysis (CCA) on the Stiefel Manifold.


import numpy as np
import matplotlib.pyplot as plt
from pymanopt import Problem
from pymanopt.manifolds import Stiefel, Product
from pymanopt.optimizers import TrustRegions
import pandas as pd
from sklearn.cross_decomposition import CCA
import time
from utils.utils import init_stability_log, log_stability_data


def standard_cca(X: np.ndarray, Y: np.ndarray, n_components: int) -> tuple:
    """
    Perform Standard CCA.

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
    tuple : (XA, YB, correlations, runtime)
        Canonical scores for X and Y, canonical correlations, and runtime.
    """
    start_time = time.time()
    cca = CCA(n_components=n_components)
    X_c, Y_c = cca.fit_transform(X, Y)
    end_time = time.time()

    # Compute canonical correlations for comparison
    correlations = [
        np.corrcoef(X_c[:, i], Y_c[:, i])[0, 1] for i in range(n_components)
    ]

    return X_c, Y_c, correlations, end_time - start_time


from scipy.linalg import cholesky, polar


def cca_objective(X: np.ndarray, Y: np.ndarray, A: np.ndarray, B: np.ndarray) -> float:
    """
    Compute the objective function for CCA.

    Parameters:
    ----------
    X : np.ndarray
        The first dataset (genes with expressions).
    Y : np.ndarray
        The second dataset (pathways of the genes).
    A : np.ndarray
        The first projection matrix.
    B : np.ndarray
        The second projection matrix.
    """
    XA = X @ A
    YB = Y @ B
    corr = np.sum(XA * YB) / (n - 1)
    return -corr  # negative because we're maximizing


def cca_gradient(X: np.ndarray, Y: np.ndarray, A: np.ndarray, B: np.ndarray) -> tuple:
    """
    Compute the gradient of the objective function for Canonical Correlation Analysis (CCA).
    """
    n_samples = X.shape[0]

    # Check matrix dimensions
    if X.shape[1] != A.shape[0] or Y.shape[1] != B.shape[0] or A.shape[1] != B.shape[1]:
        raise ValueError(
            "Dimension mismatch: Please ensure X, Y, A, and B are compatible with shape requirements."
        )

    # Compute projections
    XA = X @ A  # Shape: (n_samples, k)
    YB = Y @ B  # Shape: (n_samples, k)

    # Compute gradients
    grad_A = -X.T @ YB / (n_samples - 1)
    grad_B = -Y.T @ XA / (n_samples - 1)

    return grad_A, grad_B


def cholesky_qr_retraction(
    X: np.ndarray, G: np.ndarray, xi: np.ndarray, epsilon=1e-6
) -> np.ndarray:
    """
    Retract a point on the Stiefel manifold using the Cholesky QR retraction. Add regularization to avoid singular matrices.
    """
    Z = (X + xi).T @ G @ (X + xi)
    Z += epsilon * np.eye(Z.shape[0])  # Add regularization
    try:
        L = cholesky(Z, lower=True)
    except np.linalg.LinAlgError:
        # If Cholesky fails, fall back to Polar decomposition
        return polar_retraction(X, G, xi)
    retracted_point = (X + xi) @ np.linalg.inv(L.T)
    return retracted_point


def polar_retraction(X: np.ndarray, G: np.ndarray, xi: np.ndarray) -> np.ndarray:
    """
    Retract a point on the Stiefel manifold using the polar retraction.
    """
    Z = (X + xi).T @ G @ (X + xi)
    U, _ = polar(Z)
    retracted_point = (X + xi) @ np.linalg.inv(U.T)
    return retracted_point


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


def riemannian_gradient_descent(
    X,
    Y,
    k,
    retraction_method,
    max_iter=10,
    lr=0.8,
    tol=1e-3,
    stability_log="stability.csv",
):
    """
    Perform Riemannian gradient descent to solve CCA.

    Parameters:
    ----------
    X : np.ndarray
        The first dataset (genes with expressions).
    Y : np.ndarray
        The second dataset (pathways of the genes).
    retraction_method : str
        The retraction method to use. Either 'cholesky' or 'polar'.
    max_iter : int
        The maximum number of iterations.
    lr : float
        The learning rate.
    tol : float
        The tolerance for convergence.
    """
    n, p1 = X.shape
    _, p2 = Y.shape

    # Initialize A and B
    A = np.random.randn(p1, k)
    B = np.random.randn(p2, k)
    A, _ = np.linalg.qr(A)
    B, _ = np.linalg.qr(B)

    G_A = np.eye(p1)  # Metric for A
    G_B = np.eye(p2)  # Metric for B

    start_time = time.time()

    # Initialize the stability log file
    init_stability_log(stability_log, k)

    for i in range(max_iter):
        # Log stability data for each iteration
        log_stability_data(stability_log, i, A, B, G_A, G_B)

        old_obj = cca_objective(X, Y, A, B)
        grad_A, grad_B = cca_gradient(X, Y, A, B)

        # Clip gradients to avoid too large updates
        max_grad_norm = 1e2
        grad_A = np.clip(grad_A, -max_grad_norm, max_grad_norm)
        grad_B = np.clip(grad_B, -max_grad_norm, max_grad_norm)

        # Retraction with checks
        if retraction_method == "cholesky":
            A_new = cholesky_qr_retraction(A, G_A, -lr * grad_A)
            B_new = cholesky_qr_retraction(B, G_B, -lr * grad_B)
        elif retraction_method == "polar":
            A_new = polar_retraction(A, G_A, -lr * grad_A)
            B_new = polar_retraction(B, G_B, -lr * grad_B)

        # Check for NaNs or Infs in retracted matrices
        if (
            np.isnan(A_new).any()
            or np.isnan(B_new).any()
            or np.isinf(A_new).any()
            or np.isinf(B_new).any()
        ):
            print(
                "Warning: Retraction resulted in NaN or Inf values. Reducing learning rate."
            )
            lr *= 0.5
            continue  # Skip this iteration with a reduced learning rate

        A, B = A_new, B_new  # Update if values are valid

        # Objective convergence check
        new_obj = cca_objective(X, Y, A, B)
        if abs(new_obj - old_obj) < tol:
            break

    end_time = time.time()

    # Compute canonical correlations
    XA = X @ A
    YB = Y @ B
    correlations = np.diag(XA.T @ YB) / (n - 1)

    return A, B, correlations, end_time - start_time, XA, YB


def run_experiment(
    X: np.ndarray, Y: np.ndarray, k_values: list, retraction_methods: list
):
    results = []
    scores = []
    for k in k_values:
        for method in retraction_methods:
            A, B, correlations, runtime, XA, YB = riemannian_gradient_descent(
                X, Y, k, method
            )
            results.append(
                {
                    "k": k,
                    "method": method,
                    "correlations": correlations,
                    "runtime": runtime,
                }
            )
            scores.append({"k": k, "XA": XA, "YB": YB})
    return pd.DataFrame(results), pd.DataFrame(scores)
