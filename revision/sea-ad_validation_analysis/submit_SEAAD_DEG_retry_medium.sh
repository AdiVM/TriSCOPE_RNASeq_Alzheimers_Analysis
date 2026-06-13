#!/bin/bash
# Retry submitter: scans the output dir for any missing (celltype, set, chunk)
# CSV and resubmits ONLY those jobs on the MEDIUM partition (5d).
#
# Handles:
#   - 5 non-Ex cell types (set_id always 1)
#   - Ex with 3 donor-split sets (set_id 1/2/3)
#
# Usage:
#   bash submit_SEAAD_DEG_retry_medium.sh            # actually resubmit
#   bash submit_SEAAD_DEG_retry_medium.sh --dry-run  # just report missing

set -u

N_CHUNKS=100
NON_EX_CELLTYPES=("Opc" "Ast" "Mic" "Oli" "Inh")
EX_SETS=(1 2 3)

OUT_DIR="/n/groups/patel/adithya/Alz_Outputs/Final_Outputs/seaad_DEG_results"
RUN_SCRIPT="/home/adm808/Revision_nat_comms/fresh_later/run_SEAAD_DEG_single_medium.sh"
MANIFEST="${OUT_DIR}/retry_manifest_$(date +%Y%m%d_%H%M%S).txt"

mkdir -p "$OUT_DIR"

DRY_RUN=0
if [ "${1:-}" = "--dry-run" ]; then
  DRY_RUN=1
  echo "DRY-RUN: nothing will be submitted."
fi

echo "Scanning $OUT_DIR for missing chunks (N_CHUNKS=$N_CHUNKS)..."
{
  echo "# Retry manifest generated $(date)"
  echo "# N_CHUNKS=$N_CHUNKS"
  echo "# non-Ex celltypes: ${NON_EX_CELLTYPES[*]}"
  echo "# Ex sets: ${EX_SETS[*]}"
  echo "# Columns: celltype set_id chunk_id status"
} > "$MANIFEST"

total_missing=0
total_present=0

# ---- non-Ex (set_id=1) ----
for ct in "${NON_EX_CELLTYPES[@]}"; do
  ct_missing=0
  ct_present=0
  for chunk_id in $(seq 1 $N_CHUNKS); do
    expected="${OUT_DIR}/poisson_DE_results_SEAAD_${ct}_set1_chunk${chunk_id}_of_${N_CHUNKS}.csv"
    if [ -s "$expected" ]; then
      ct_present=$((ct_present + 1))
      echo "$ct 1 $chunk_id PRESENT" >> "$MANIFEST"
    else
      ct_missing=$((ct_missing + 1))
      echo "$ct 1 $chunk_id MISSING" >> "$MANIFEST"
      if [ $DRY_RUN -eq 0 ]; then
        sbatch "$RUN_SCRIPT" "$ct" "$chunk_id" "$N_CHUNKS"
      fi
    fi
  done
  echo "  $ct: present=$ct_present missing=$ct_missing"
  total_present=$((total_present + ct_present))
  total_missing=$((total_missing + ct_missing))
done

# ---- Ex (set_id=1/2/3) ----
for s in "${EX_SETS[@]}"; do
  ex_missing=0
  ex_present=0
  for chunk_id in $(seq 1 $N_CHUNKS); do
    expected="${OUT_DIR}/poisson_DE_results_SEAAD_Ex_set${s}_chunk${chunk_id}_of_${N_CHUNKS}.csv"
    if [ -s "$expected" ]; then
      ex_present=$((ex_present + 1))
      echo "Ex $s $chunk_id PRESENT" >> "$MANIFEST"
    else
      ex_missing=$((ex_missing + 1))
      echo "Ex $s $chunk_id MISSING" >> "$MANIFEST"
      if [ $DRY_RUN -eq 0 ]; then
        sbatch "$RUN_SCRIPT" "Ex" "$chunk_id" "$N_CHUNKS" "$s"
      fi
    fi
  done
  echo "  Ex set$s: present=$ex_present missing=$ex_missing"
  total_present=$((total_present + ex_present))
  total_missing=$((total_missing + ex_missing))
done

echo ""
echo "Summary: present=$total_present  missing=$total_missing  total=$((total_present + total_missing))"
echo "Manifest written to: $MANIFEST"

if [ $DRY_RUN -eq 0 ] && [ $total_missing -gt 0 ]; then
  echo "Resubmitted $total_missing jobs on partition=medium."
fi
