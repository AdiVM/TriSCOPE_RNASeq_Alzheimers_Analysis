#!/bin/bash
#SBATCH --job-name=rosmap427_norm
#SBATCH --account=patel
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=/home/adm808/Revision_nat_comms/fresh_later/rosmap_427/rosmap427_norm_%j.out
#SBATCH --error=/home/adm808/Revision_nat_comms/fresh_later/rosmap_427/rosmap427_norm_%j.err

# Same R environment the PFC427 DEG ran in (reads Seurat RDS, writes parquet via arrow)
module load gcc/14.2.0
module load R/4.4.2
export R_LIBS_USER=$HOME/R/library

Rscript /home/adm808/Revision_nat_comms/fresh_later/rosmap_427/rosmap_normalization_matrix.R
