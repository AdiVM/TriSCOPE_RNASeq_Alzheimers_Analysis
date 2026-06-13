#!/bin/bash
#SBATCH --job-name=SEAAD_DEG
#SBATCH --account=patel
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=249G
#SBATCH --time=5-00:00:00
#SBATCH --output=SEAAD_DEG_%x_%j.out
#SBATCH --error=SEAAD_DEG_%x_%j.err

# ----------------------------
# Arguments
# ----------------------------
CELLTYPE=$1
CHUNK_ID=$2
N_CHUNKS=$3
SET_ID=$4

echo "Running DEG for cell type: $CELLTYPE"
echo "Chunk $CHUNK_ID of $N_CHUNKS"

# ----------------------------
# Modules
# ----------------------------
module load gcc/14.2.0
module load R/4.4.2

export R_LIBS_USER=$HOME/R/library

# ----------------------------
# Run
# ----------------------------
if [ -z "$SET_ID" ]; then
  Rscript /home/adm808/Revision_nat_comms/fresh_later/DEG_SEAAD_chunking.R \
    "$CELLTYPE" "$CHUNK_ID" "$N_CHUNKS"
else
  Rscript /home/adm808/Revision_nat_comms/fresh_later/DEG_SEAAD_chunking.R \
    "$CELLTYPE" "$CHUNK_ID" "$N_CHUNKS" "$SET_ID"
fi
