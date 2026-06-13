#!/bin/bash
# Submit SEA-AD DEG chunked jobs: 100 chunks per cell type.

N_CHUNKS=100
CELLTYPES=("Opc" "Ast" "Mic" "Oli" "Inh" "Ex")

RUN_SCRIPT=/home/adm808/Revision_nat_comms/fresh_later/run_SEAAD_DEG_single.sh

for ct in "${CELLTYPES[@]}"; do
  for chunk_id in $(seq 1 $N_CHUNKS); do
    sbatch "$RUN_SCRIPT" "$ct" "$chunk_id" "$N_CHUNKS"
  done
done
