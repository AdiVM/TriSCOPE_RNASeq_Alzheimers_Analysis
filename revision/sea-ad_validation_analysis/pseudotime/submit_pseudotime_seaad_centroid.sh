#!/bin/bash
# Submit one SLURM job per Mathys subcluster (41 jobs total when called with no args).
# Per-broad memory matches the labeling-pipeline budget because the dominant memory
# cost is loading the SEA-AD per-broad expression matrix.
#
# Usage:
#   ./submit_pseudotime_seaad_centroid.sh             # submit all 41 subclusters
#   ./submit_pseudotime_seaad_centroid.sh Ex5 In3     # submit only the listed ones

set -euo pipefail

RUN=/home/adm808/Revision_nat_comms/fresh_later/pseudotime/run_pseudotime_seaad_centroid.sh

# Broad-level memory budgets (same as submit_label_seaad_centroid.sh).
declare -A BROAD_MEM=( [Ex]=249G [Inh]=160G [Oli]=160G [Ast]=140G [Mic]=124G [Opc]=116G )

# All 41 Mathys subpops (centroid pipeline labels SEA-AD cells with these names).
ALL_SUBCLUSTERS=(
  Ex0 Ex1 Ex2 Ex3 Ex4 Ex5 Ex6 Ex7 Ex8 Ex9 Ex11 Ex12 Ex14
  In0 In1 In2 In3 In4 In5 In6 In7 In8 In9 In10 In11
  Ast0 Ast1 Ast2 Ast3
  Oli0 Oli1 Oli3 Oli4 Oli5
  Mic0 Mic1 Mic2 Mic3
  Opc0 Opc1 Opc2
)

declare -A SUB_TO_BROAD=(
  [Ex0]=Ex [Ex1]=Ex [Ex2]=Ex [Ex3]=Ex [Ex4]=Ex [Ex5]=Ex [Ex6]=Ex
  [Ex7]=Ex [Ex8]=Ex [Ex9]=Ex [Ex11]=Ex [Ex12]=Ex [Ex14]=Ex
  [In0]=Inh [In1]=Inh [In2]=Inh [In3]=Inh [In4]=Inh [In5]=Inh
  [In6]=Inh [In7]=Inh [In8]=Inh [In9]=Inh [In10]=Inh [In11]=Inh
  [Ast0]=Ast [Ast1]=Ast [Ast2]=Ast [Ast3]=Ast
  [Oli0]=Oli [Oli1]=Oli [Oli3]=Oli [Oli4]=Oli [Oli5]=Oli
  [Mic0]=Mic [Mic1]=Mic [Mic2]=Mic [Mic3]=Mic
  [Opc0]=Opc [Opc1]=Opc [Opc2]=Opc
)

if [ "$#" -gt 0 ]; then
  TARGETS=("$@")
else
  TARGETS=("${ALL_SUBCLUSTERS[@]}")
fi

for sub in "${TARGETS[@]}"; do
  broad="${SUB_TO_BROAD[$sub]:-}"
  if [ -z "$broad" ]; then
    echo "Unknown subcluster: $sub (valid: ${ALL_SUBCLUSTERS[*]})" >&2
    exit 2
  fi
  mem="${BROAD_MEM[$broad]}"
  echo "Submitting $sub (broad=$broad, mem=$mem)"
  sbatch \
    --job-name="pseudotime_${sub}" \
    --partition=short \
    --mem="$mem" \
    --time=4:00:00 \
    --output="pseudotime_${sub}_%j.out" \
    --error="pseudotime_${sub}_%j.err" \
    "$RUN" "$sub"
done
