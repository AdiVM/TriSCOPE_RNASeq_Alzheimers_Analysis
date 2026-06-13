#!/bin/bash
#SBATCH --job-name=seaad_extract
#SBATCH --output=seaad_extract_%J.out
#SBATCH --error=seaad_extract_%J.err
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --partition=short
#SBATCH --mem=200G

source /n/groups/patel/adithya/scenv/bin/activate

python /home/adm808/Revision_nat_comms/fresh_later/sea_ad_normalization.py