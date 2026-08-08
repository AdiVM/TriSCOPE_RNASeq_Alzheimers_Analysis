# TriSCOPE_RNASeq_Alzheimers_Analysis
We developed a machine learning pipeline for single-nucleus RNA-seq data to predict Alzheimer’s disease status across diverse brain cell types. Key approaches include predictive modeling, differential expression analysis, and trajectory inference.


Contents:
1. PM_Gex.py – Predictive modeling using gene expression features only.
2. PM_Multimodal.py – Predictive modeling using both gene expression and demographic metadata.
3. Differential_Expression.R – DEG analysis (full methods described in paper).
4. Trajectory_Inference.py – Diffusion pseudotime–based trajectory inference (run with 7 random seeds).
5. figure{N}.ipynb – Jupyter notebooks for generating figures from the manuscript. Each notebook corresponds to the numbered figure in the main or supplemental text.
6. main_tables.ipynb - Notebook for generating donor and cohort summary.
8. README.md – Overview of repository.

Expected Outputs:
1. Differential Expression (Random-Effects GLMM)
For each major cell type, the differential expression pipeline outputs a CSV file (poisson_DE_results_PFC_<celltype>.csv) containing gene-level effect sizes, standard errors, p-values, FDR-adjusted p-values, log2 fold changes, and DEG calls. These results are derived from Poisson generalized linear mixed models with a donor-level random intercept and RUV-based latent covariates, enabling donor-aware inference while controlling for technical confounding.

2. Predictive Modeling (Multimodal and Gene-Only)
Predictive modeling scripts output serialized model objects (joblib files) containing trained classifiers, feature importance rankings, and metadata required for downstream analysis. These objects also store train/test predicted probabilities across multiple random splits, enabling reproducible ROC, calibration, and feature-importance analyses without retraining.

3. Trajectory Inference (Diffusion Pseudotime)
Trajectory inference outputs a set of summary CSV files with one entry per subcluster, each containing diffusion pseudotime statistics and the top 1% of genes most strongly associated with trajectory progression. These summaries enable direct comparison of trajectory-associated genes across subclusters and random seeds without requiring access to full cell-level embeddings.

4. GWAS–eQTL Colocalization (coloc)
The colocalization pipeline outputs a single CSV file (coloc_results_AD_GWAS_GTEx_v10_*.csv) summarizing posterior probabilities (PP0–PP4), SNP counts, genomic coordinates, tissue context, and gene–cell-type associations for each tested locus. This file provides a compact, gene-centric view of shared genetic architecture between Alzheimer’s disease GWAS signals and GTEx brain eQTLs, suitable for downstream filtering and enrichment analyses.


Experimental Design
1. Predictive Modeling (PM): Each PM script was run five times with random seeds 1–5. All reported metrics reflect results aggregated across these five splits.
2. Differential Expression: Conducted using Poisson GLMMs, controlling for technical and donor covariates. Exact details of model specification and thresholds are provided in the Methods section of the manuscript.
3. Trajectory Inference: Trajectories were computed with seven random seeds, as specified in Trajectory_Inference.py.

## Installation

Requires Python 3.9.21 and R 4.4.2 (R 4.5.1 for the colocalization step).
Pinned package versions for all three environments are listed in requirements.txt.

    git clone https://github.com/AdiVM/TriSCOPE_RNASeq_Alzheimers_Analysis.git
    cd TriSCOPE_RNASeq_Alzheimers_Analysis
    pip install -r requirements.txt

R packages are installed from CRAN and Bioconductor at the versions listed in
requirements.txt.

Input and output paths are defined at the top of each script and must be
updated to point at local copies of the data.

## Running the code

Predictive modeling, per cell type and per split (split index 1 to 5):

    python PM_Gex.py --exp_type maximal --cell_type <cell_type> --split_index 1
    python PM_Multimodal.py --exp_type maximal --cell_type <cell_type> --split_index 1

Differential expression:

    Rscript Differential_Expression.R

Trajectory inference:

    python Trajectory_Inference.py

GWAS-eQTL colocalization:

    Rscript revision/code/Coloc.R

Figures and tables are produced by the Jupyter notebooks in the repository
root, which read the outputs written by the scripts above.

## Test data

triscope_synthetic_demo.zip contains a synthetic microglial expression matrix
and donor metadata under Data/, together with the expected model outputs under
Expected_Outputs/ and Analysis_Results/, so the predictive modeling pipeline
can be run end to end without access to the controlled-access cohorts.

 
