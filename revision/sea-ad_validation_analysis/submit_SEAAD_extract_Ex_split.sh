#!/bin/bash
#SBATCH --job-name=seaad_ex_split
#SBATCH --account=patel
#SBATCH --output=seaad_ex_split_%J.out
#SBATCH --error=seaad_ex_split_%J.err
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --partition=short
#SBATCH --mem=248G

source /n/groups/patel/adithya/scenv/bin/activate

python /home/adm808/Revision_nat_comms/fresh_later/SEAAD_extract_Ex_split.py
