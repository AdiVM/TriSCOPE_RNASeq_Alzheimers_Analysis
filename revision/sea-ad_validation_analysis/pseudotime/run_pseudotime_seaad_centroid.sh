#!/bin/bash
# Per-subcluster pseudotime runner. Resource flags (--partition, --mem, --time,
# --job-name, --output, --error) are supplied by submit_pseudotime_seaad_centroid.sh
# and override these fallback defaults.

#SBATCH --account=sabatini
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=pseudotime_%x_%j.out
#SBATCH --error=pseudotime_%x_%j.err

set -euo pipefail

SUBCLUSTER="${1:?usage: run_pseudotime_seaad_centroid.sh <subcluster>}"

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-NA} subcluster=${SUBCLUSTER}"
echo "[$(date)] SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE:-NA} partition=${SLURM_JOB_PARTITION:-NA}"

# Pin BLAS / OpenMP thread pools to the SLURM allocation.
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

# Reproducibility: fix Python hash randomization.
export PYTHONHASHSEED=0

source /n/groups/patel/adithya/scenv/bin/activate

python /home/adm808/Revision_nat_comms/fresh_later/pseudotime/pseudotime_seaad_centroid.py \
  --subcluster "$SUBCLUSTER"

echo "[$(date)] done subcluster=${SUBCLUSTER}"
