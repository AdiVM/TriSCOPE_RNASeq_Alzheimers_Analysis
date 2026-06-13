#!/bin/bash
# First-pass submitter: every (celltype, chunk) on the SHORT partition (12h).
# Anything that times out or fails can be retried later via submit_SEAAD_DEG_retry_medium.sh.

set -u

N_CHUNKS=100
CELLTYPES=("Opc" "Ast" "Mic" "Oli" "Inh")
# NOTE: "Ex" excluded — single-MTX is 4.4B nnz, exceeds R sparse 32-bit index limit.
# Handle Ex separately after re-extracting into per-set MTXs.

RUN_SCRIPT=/home/adm808/Revision_nat_comms/fresh_later/run_SEAAD_DEG_single_short.sh

n_submitted=0
for ct in "${CELLTYPES[@]}"; do
  for chunk_id in $(seq 1 $N_CHUNKS); do
    sbatch "$RUN_SCRIPT" "$ct" "$chunk_id" "$N_CHUNKS"
    n_submitted=$((n_submitted + 1))
  done
done

echo "Submitted $n_submitted jobs (N_CHUNKS=$N_CHUNKS, celltypes=${CELLTYPES[*]})."
