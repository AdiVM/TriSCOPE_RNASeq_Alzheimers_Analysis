#!/bin/bash
# Submit specific Ex chunks for one donor-split set on SHORT.
# Carbon copy of submit_SEAAD_DEG_Ex_short.sh, but takes <set> <chunk_id> [<chunk_id> ...]
# instead of looping over all sets and chunks.
# Uses the fairshare run script (mem=124G, time=4:00:00).
#
# Usage:
#   ./submit_SEAAD_DEG_Ex_short_chunks.sh 2 47
#   ./submit_SEAAD_DEG_Ex_short_chunks.sh 2 47 48 51 73
#   ./submit_SEAAD_DEG_Ex_short_chunks.sh 2 $(seq 47 100)

set -u

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <set> <chunk_id> [<chunk_id> ...]"
  exit 1
fi

SET_ID="$1"
shift
CHUNKS=("$@")

N_CHUNKS=100
RUN_SCRIPT=/home/adm808/Revision_nat_comms/fresh_later/run_SEAAD_DEG_single_short_fairshare.sh

# Pre-flight: confirm this set's MTX exists and is non-empty.
mtx="/n/scratch/users/a/adm808/SEAAD_DEG_input/Ex_set${SET_ID}_counts.mtx"
if [ ! -s "$mtx" ]; then
  echo "ERROR: $mtx missing or empty. Run SEAAD_extract_Ex_split.py first."
  exit 1
fi

n_submitted=0
for chunk_id in "${CHUNKS[@]}"; do
  sbatch "$RUN_SCRIPT" "Ex" "$chunk_id" "$N_CHUNKS" "$SET_ID"
  n_submitted=$((n_submitted + 1))
done

echo "Submitted $n_submitted Ex jobs (SET=$SET_ID, N_CHUNKS=$N_CHUNKS, chunks=${CHUNKS[*]})."
