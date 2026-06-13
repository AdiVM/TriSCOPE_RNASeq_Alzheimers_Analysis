#!/bin/bash
#SBATCH --job-name=all_rosmap
#SBATCH --output=/n/scratch/users/a/adm808/slurm_logs/all_%J.out
#SBATCH --error=/n/scratch/users/a/adm808/slurm_logs/all_%J.err
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --partition=medium
#SBATCH --mem=200G

# Load modules (modify if necessary)
#module load gcc/9.2.0
#module load graphviz/3.0.0

# conda init bash
#source ~/.bashrc
#conda activate scenv
#python -m venv scenv
source /n/groups/patel/adithya/scenv/bin/activate

# Set PYTHONUNBUFFERED to ensure immediate flushing of print statements
export PYTHONUNBUFFERED=1


# running a single cell_type
experiment=$1
# sample_number=$2
ct=$2
split_number=$3
echo "Running leave one out experiment with experiment $experiment and sample: $ct and $split_number"
#echo "Running Maximal experiment"
#python3 runAutoML_celltype.py --cell_type "$cell_type"
# python3 runAutoML_celltype_Syn16.py --cell_type "$cell_type"
python3 /home/adm808/Revision_nat_comms/fresh_later/rosmap_427/Multirun_cell_on_cell_genes_all_rosmap.py --exp_type "$experiment" --cell_type "$ct" --split_index "$split_number"
#--cell_type "$ct"


