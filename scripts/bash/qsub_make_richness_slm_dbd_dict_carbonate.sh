#!/bin/bash
#SBATCH --mail-user=williamrshoemaker@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=dbd-richness
#SBATCH --mem=20G


module unload python
module load anaconda
conda activate macroeco_phylo


python /N/u/wrshoema/Carbonate/macroeco_phylo/scripts/make_richness_slm_dbd_dict_carbonate.py




# sbatch qsub_make_phylo_dict_carbonate.sh

# squeue -u wrshoema
