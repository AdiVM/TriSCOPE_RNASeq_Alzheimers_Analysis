#!/bin/bash
#SBATCH --job-name=cognitive_tests
#SBATCH --output=/n/scratch/users/a/adm808/slurm_logs/job_%x_%A_%a.out
#SBATCH --error=/n/scratch/users/a/adm808/slurm_logs/job_%x_%A_%a.err
#SBATCH --nodes=1
#SBATCH --time=08:00:00
#SBATCH --partition=short
#SBATCH --mem=64G

# Activate your environment
source /n/groups/patel/adithya/scenv/bin/activate
export PYTHONUNBUFFERED=1

# Define the experiment type
experiment="maximal"

# Define the list of cell types and number of splits
cell_types=("Inh" "Ex")
#  "Oli" "Ast" "Mic" "Opc" removing these cuz they were run add them back if you rerun this
num_splits=5

#Loop over each cell type and each split
for ct in "${cell_types[@]}"; do
    for split in $(seq 1 $num_splits); do
        echo "Submitting job for $ct, split $split"
        sbatch --account=patel ein_Multirun_cell_on_cell_genes_seaad.sh "$experiment" "$ct" "$split"
        sleep 1
    done
done


echo "All jobs submitted."
