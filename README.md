# TriSCOPE_RNASeq_Alzheimers_Analysis
We developed a machine learning pipeline for single-nucleus RNA-seq data to predict Alzheimer’s disease status across diverse brain cell types. Key approaches include predictive modeling, differential expression analysis, and trajectory inference.


Contents:
	•	PM_Gex.py – Predictive modeling using gene expression features only.
	•	PM_Multimodal.py – Predictive modeling using both gene expression and demographic metadata.
	•	Differential_Expression.R – DEG analysis (full methods described in paper).
	•	Trajectory_Inference.py – Diffusion pseudotime–based trajectory inference (run with 7 random seeds).
	•	figure{N}.ipynb – Jupyter notebooks for generating figures from the manuscript. Each notebook corresponds to the numbered figure in the main or supplemental text.
  • main_tables.ipynb - Notebook for generating donor and cohort summary.
	•	data_tables.ipynb / Validation_data_tables.ipynb – Scripts for generating and formatting supplemental tables.
	•	README.md – Overview of repository.


Experimental Design
	•	Predictive Modeling (PM):
	  1. Each PM script was run five times with random seeds 1–5.
	  2. All reported metrics reflect results aggregated across these five splits.
	•	Differential Expression:
	  1. Conducted using Poisson GLMMs, controlling for technical and donor covariates.
	  2. Exact details of model specification and thresholds are provided in the Methods section of the manuscript.
	•	Trajectory Inference:
	  1. Trajectories were computed with seven random seeds, as specified in Trajectory_Inference.py.

 
