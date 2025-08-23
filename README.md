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


Experimental Design
1. Predictive Modeling (PM): Each PM script was run five times with random seeds 1–5. All reported metrics reflect results aggregated across these five splits.
2. Differential Expression: Conducted using Poisson GLMMs, controlling for technical and donor covariates. Exact details of model specification and thresholds are provided in the Methods section of the manuscript.
3. Trajectory Inference: Trajectories were computed with seven random seeds, as specified in Trajectory_Inference.py.

 
