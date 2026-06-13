#!/bin/bash
# Submit all Ex chunks across 3 donor-split sets on SHORT (12h).
# Mirrors 427 Ex_set1/set2/set3 architecture.
#
# Run AFTER extraction (submit_SEAAD_extract_Ex_split.sh) completes successfully.

set -u

N_CHUNKS=100
N_SETS=3
RUN_SCRIPT=/home/adm808/Revision_nat_comms/fresh_later/run_SEAAD_DEG_single_short.sh

# Pre-flight: confirm each set's MTX exists and is non-empty.
for s in $(seq 1 $N_SETS); do
  mtx="/n/scratch/users/a/adm808/SEAAD_DEG_input/Ex_set${s}_counts.mtx"
  if [ ! -s "$mtx" ]; then
    echo "ERROR: $mtx missing or empty. Run SEAAD_extract_Ex_split.py first."
    exit 1
  fi
done

n_submitted=0
for s in $(seq 1 $N_SETS); do
  for chunk_id in $(seq 1 $N_CHUNKS); do
    sbatch "$RUN_SCRIPT" "Ex" "$chunk_id" "$N_CHUNKS" "$s"
    n_submitted=$((n_submitted + 1))
  done
done

echo "Submitted $n_submitted Ex jobs (N_SETS=$N_SETS, N_CHUNKS=$N_CHUNKS)."
