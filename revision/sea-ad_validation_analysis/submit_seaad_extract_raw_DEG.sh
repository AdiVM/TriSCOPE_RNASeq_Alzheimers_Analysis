#!/bin/bash
#SBATCH --job-name=seaad_deg_extract
#SBATCH --output=seaad_deg_extract_%J.out
#SBATCH --error=seaad_deg_extract_%J.err
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --partition=short
#SBATCH --mem=248G

source /n/groups/patel/adithya/scenv/bin/activate

python /home/adm808/Revision_nat_comms/fresh_later/SEAAD_extract_raw_for_DEG.py