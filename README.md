# **Multivariate Reduced-Rank Models**  
Welcome to the **Multivariate Reduced-Rank Models** repository! This project explores **Canonical Correlation Analysis (CCA)** in machine learning, focusing on numerical methods and optimization on matrix manifolds. It's a continuation of a prior project (2023) on Principal Component Analysis (PCA), expanding into CCA for high-dimensional data and data visualization.

---

## **🌟 Goals of the Project**  
1. **CCA for High-Dimensional Data**: Understand and apply CCA to explore relationships between datasets with many features.
2. **Optimization on Matrix Manifolds**: Use **Cholesky QR-based retraction** to improve the efficiency of manifold-based optimization methods. I use the Stiefel manifold for orthogonalization reasons.
3. **CCA for Cancer Cell Pathways Classification**: Apply CCA for meaningful visualizations of multi-view datasets and demonstrate code applications.

The dataset contains gene expression data from breast cancer cell lines and corresponding KEGG pathway information. The data is derived from the GSE48213 dataset, which includes gene expression profiles of various breast cancer cell lines.
- Gene Expression Data:Genes (Ensembl IDs) x Cell lines (e.g., GSM1172844_184A1, GSM1172845_184B5, etc.)
- [KEGG Pathway Information](https://www.genome.jp/kegg/pathway.html): Genes (Ensembl IDs) x KEGG pathways

---

## **🔧 Foundations**  
This project builds on the following key concepts:
- **Canonical Correlation Analysis (CCA)**: A method to identify and quantify relationships between two sets of high-dimensional variables.
- **Matrix Manifolds**: Focus on the **Stiefel manifold**, a space of orthonormal matrices.
- **Riemannian Optimization**: Optimization methods on curved spaces using retractions to maintain feasibility.
- **Cholesky Retraction**: A more efficient way of projecting points back onto a manifold compared to traditional methods like polar decomposition.

---

## **📚 Mathematical Background**  
This project revolves around the **generalized Stiefel manifold** \(\text{St}(n, p, B) = \{ X \in \mathbb{R}^{n \times p} : X^T B X = I_p\}\), where **\(B\)** is a positive definite matrix.  
### Key Equations:
- **CCA**: Maximizes the correlation between two sets of canonical variables:
  \[
  \text{max}_{A, B} \, \text{corr}(X A, Y B) \quad \text{such that} \quad A^T A = I, \, B^T B = I
  \]
- **Cholesky Retraction**: Retracts points from the tangent space to the manifold by:
  \[
  R_X(\xi) = X L^{-T}
  \]
  where **\(L\)** is the Cholesky factor obtained from \(Y = I + \xi^T B \xi\).

---

## **💻 Requirements**  
To run the code and experiments, make sure you have the following installed:
- **Python 3.11**
- **NumPy**: For numerical operations.
- **SciPy**: For matrix decompositions and linear algebra routines.
- **Matplotlib**: For plotting and visualization.
- **Pandas**: For handling datasets.
- **PyManopt**: For manifold optimization.

You can install the dependencies with:
```bash
pip install numpy scipy matplotlib pandas pymanopt

