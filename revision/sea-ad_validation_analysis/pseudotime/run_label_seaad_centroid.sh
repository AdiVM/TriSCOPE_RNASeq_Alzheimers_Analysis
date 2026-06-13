#!/bin/bash
# Stage A runner (centroid variant): score one SEA-AD broad cell type by Pearson
# correlation to Mathys per-subpop centroids. Sibling of run_label_seaad_subpops.sh.
# Resource flags (--partition, --mem, --time, --job-name, --output, --error)
# are supplied per-broad by submit_label_seaad_centroid.sh and override these
# fallback defaults.

#SBATCH --account=sabatini
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=seaad_centroid_%x_%j.out
#SBATCH --error=seaad_centroid_%x_%j.err

set -euo pipefail

BROAD="${1:?usage: run_label_seaad_centroid.sh <broad_cell_type>}"

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-NA} broad=${BROAD}"
echo "[$(date)] SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE:-NA} partition=${SLURM_JOB_PARTITION:-NA}"

# Pin BLAS / OpenMP thread pools to the SLURM allocation so numpy doesn't
# oversubscribe (default is host core count, not cpus-per-task).
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

# Reproducibility: fix Python hash randomization so dict / set iteration
# orders are bit-identical across re-runs (matters for any future code path
# that iterates a set; current code uses sets only for `in` membership tests).
export PYTHONHASHSEED=0

source /n/groups/patel/adithya/scenv/bin/activate

python /home/adm808/Revision_nat_comms/fresh_later/pseudotime/label_seaad_subpops_centroid.py \
  --broad_cell_type "$BROAD"

echo "[$(date)] done broad=${BROAD}"
