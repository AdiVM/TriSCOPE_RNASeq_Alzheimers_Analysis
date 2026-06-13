#!/bin/bash
# Stage A runner: score one SEA-AD broad cell type against Mathys subpops.
# Resource flags (--partition, --mem, --time, --job-name, --output, --error)
# are supplied per-broad by submit_label_seaad_subpops.sh and override these
# fallback defaults.

#SBATCH --account=sabatini
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=seaad_label_%x_%j.out
#SBATCH --error=seaad_label_%x_%j.err

set -euo pipefail

BROAD="${1:?usage: run_label_seaad_subpops.sh <broad_cell_type>}"

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-NA} broad=${BROAD}"
echo "[$(date)] SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE:-NA} partition=${SLURM_JOB_PARTITION:-NA}"

# Pin BLAS / OpenMP thread pools to the SLURM allocation so numpy doesn't
# oversubscribe (default is host core count, not cpus-per-task).
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

source /n/groups/patel/adithya/scenv/bin/activate

python /home/adm808/Revision_nat_comms/fresh_later/pseudotime/label_seaad_subpops.py \
  --broad_cell_type "$BROAD"

echo "[$(date)] done broad=${BROAD}"
