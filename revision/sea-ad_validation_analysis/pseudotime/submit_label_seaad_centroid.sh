#!/bin/bash
# Stage A submitter (centroid variant): one SLURM job per SEA-AD broad cell type.
# Sibling of submit_label_seaad_subpops.sh; same per-broad mem budgets observed
# on the 2026-05-20 score_genes run plus headroom for the intermediate X_shared
# allocation in the centroid pipeline.
#
# Usage:
#   ./submit_label_seaad_centroid.sh             # submits all six
#   ./submit_label_seaad_centroid.sh Ex Oli      # submits only listed broads

set -euo pipefail

RUN=/home/adm808/Revision_nat_comms/fresh_later/pseudotime/run_label_seaad_centroid.sh

# broad  partition  mem    time
declare -A PART=( [Ex]=short [Inh]=short [Oli]=short [Ast]=short [Mic]=short [Opc]=short )
declare -A MEM=(  [Ex]=249G  [Inh]=160G  [Oli]=160G   [Ast]=140G   [Mic]=124G   [Opc]=116G  )
declare -A TIME=( [Ex]=12:00:00 [Inh]=12:00:00 [Oli]=12:00:00 [Ast]=12:00:00 [Mic]=12:00:00 [Opc]=12:00:00 )

if [ "$#" -gt 0 ]; then
  BROADS=("$@")
else
  BROADS=(Ex Inh Oli Ast Mic Opc)
fi

for broad in "${BROADS[@]}"; do
  if [ -z "${PART[$broad]:-}" ]; then
    echo "Unknown broad: $broad (valid: Ex Inh Oli Ast Mic Opc)" >&2
    exit 2
  fi
  echo "Submitting $broad -> partition=${PART[$broad]} mem=${MEM[$broad]} time=${TIME[$broad]}"
  sbatch \
    --job-name="seaad_centroid_${broad}" \
    --partition="${PART[$broad]}" \
    --mem="${MEM[$broad]}" \
    --time="${TIME[$broad]}" \
    --output="seaad_centroid_${broad}_%j.out" \
    --error="seaad_centroid_${broad}_%j.err" \
    "$RUN" "$broad"
done
