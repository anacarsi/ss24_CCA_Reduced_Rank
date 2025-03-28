# **Manifold-Based CCA for Gene-Drug Sensitivity**  
Welcome to this repository:) This project presents **Canonical Correlation Analysis (CCA)** from the prespective of optimization on matrix manifolds. It's a continuation of a prior project (2023) on Principal Component Analysis (PCA), expanding into CCA for high-dimensional data and data visualization. A study comparing the final optimization methods used is presented.
 <p align="center"><img src="images/results_retraction_comparison.png" alt="Retraction comparison for GSE48213 Dataset" />
</p>

---

## **🌟 Goals of the Project**  
1. **CCA for High-Dimensional Data**: Apply CCA to find relationships between high dimensional datasets.
2. **Optimization on Matrix Manifolds**: Use **Cholesky QR-based retraction** to improve the efficiency of manifold-based optimization methods. I use the Stiefel manifold for orthogonalization reasons.
3. **CCA for Cancer Cell Pathways Classification**: Apply CCA to determine which genes suppose a predisposition to be more or less sensitive breast cancer drugs.
   
---

## **🔧 Disclaimer and Refs**  
Inspiration for this project was obtained from the work of Florian Yger et al. on [Adaptive Canonical Correlation Analysis](https://arxiv.org/pdf/1206.6453). The application of the algorithm came through research favoured by Wenxing Hu, Dongdong Lin et al. on [PubMed](https://pubmed.ncbi.nlm.nih.gov/29364120/). Researchers from the DFKZ at Heidelberg University have also contributed to the perspectives achieved in this work.  

The dataset used contains gene expression data from breast cancer cell lines and corresponding KEGG pathway information. The data is derived from the GSE48213 dataset, which includes gene expression profiles of various breast cancer cell lines.
- [Gene Expression Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48213):Genes (Ensembl IDs) x Cell lines (e.g., GSM1172844_184A1, GSM1172845_184B5, etc.)].
- [KEGG Pathway Information](https://www.genome.jp/kegg/pathway.html): Genes (Ensembl IDs) x KEGG pathways

The list of drugs against breast cancer analyzed are currently:
1. HER2 Inhibitors:
    - Example: Lapatinib.
    - Importance:
        Crucial for HER2-positive breast cancer, which is an aggressive subtype but highly treatable with targeted therapies like HER2 inhibitors.
        HER2-positive cancers make up about 15-20% of breast cancer cases.
    - Potential ranking: Very important in the context of HER2-positive patients.

2. Hormone Therapy:
    - Examples: Tamoxifen, anastrozole, letrozole, exemestane.
    - Importance:
        Dominant for hormone receptor-positive (HR+) breast cancer, which constitutes about 70-80% of breast cancer cases.
        Hormone therapies are a cornerstone in treatment for postmenopausal and premenopausal HR+ cases.
    - Potential ranking: Likely the most important due to the high prevalence of HR+ breast cancer.

3. PARP Inhibitors:
    - Example: Olaparib.
   - Importance:
        Used primarily in patients with BRCA1/BRCA2 mutations, which are less common but significant for targeted treatment.
        Effective in treating triple-negative breast cancer (TNBC) with BRCA mutations, a challenging subtype to manage.
    - Potential ranking: Important but more niche compared to HER2 inhibitors and hormone therapy.

4. CDK4/6 Inhibitors:
    - Example: Palbociclib.
    - Importance:
        Used in HR+/HER2- advanced or metastatic breast cancer in combination with hormone therapy.
        Rapidly gaining prominence in first-line treatment for advanced HR+ cases.
    - Potential ranking: Very important, especially in advanced cases of HR+ breast cancer.

5. PI3K Inhibitors:
    - Example: Alpelisib.
    - Importance:
        Targets PIK3CA-mutated HR+/HER2- advanced breast cancer, a subset of hormone receptor-positive cancers.
        Approved for patients who progress on endocrine therapy.
    - Potential ranking: Important but more targeted for specific mutations.

Due to the sparse nature of gene data, sparse CCA on the Stiefel manifold is applied to distinguish different patterns of sensitivity across cell lines.
 <p align="center"><img src="images/results_stiefel_lapatinib_alpelisib.png" alt="Projection on main CCA components for HER2 and PI3K Inhibitors" />
</p>

---


## **📚 Mathematical Background**  
This project revolves around the **generalized Stiefel manifold** \(\text{St}(n, p, B) = \{ X \in \mathbb{R}^{n \times p} : X^T B X = I_p\}\), where **\(B\)** is a positive definite matrix.  

- **CCA**: Maximizes the correlation between two sets of canonical variables.

- **Cholesky Retraction**: Retracts points from the tangent space to the manifold by:

  $$
  R_X(\xi) = X L^{-T}
  $$

  where **\( L \)** is the Cholesky factor obtained from:

  $$
  Y = I + \xi^T B \xi
  $$
  
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
pip install -e .
