#!/bin/bash
#SBATCH --mail-user=williamrshoemaker@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=test
#SBATCH --mem=20G


module unload python
module load anaconda
conda activate macroeco_phylo


python /N/u/wrshoema/Carbonate/macroeco_phylo/scripts/make_phylo_dict_carbonate.py 0.047729609235051915



# sbatch qsub_make_phylo_dict_carbonate.sh

# squeue -u wrshoema
