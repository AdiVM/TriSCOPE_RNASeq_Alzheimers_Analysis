#!/bin/bash
#SBATCH --job-name=rosmap427_inspect
#SBATCH --account=patel
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=160G
#SBATCH --time=02:00:00
#SBATCH --output=/home/adm808/Revision_nat_comms/fresh_later/rosmap_427/rosmap427_inspect_%j.out
#SBATCH --error=/home/adm808/Revision_nat_comms/fresh_later/rosmap_427/rosmap427_inspect_%j.err

# Same R environment the PFC427 DEG ran in
module load gcc/14.2.0
module load R/4.4.2
export R_LIBS_USER=$HOME/R/library

Rscript /home/adm808/Revision_nat_comms/fresh_later/rosmap_427/inspect_rosmap427.R
