#!/bin/bash
#SBATCH --mail-user=williamrshoemaker@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=96:00:00
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=sim-slm-test
#SBATCH --mem=80G


module unload python
module load anaconda
conda activate macroeco_phylo



python /N/u/wrshoema/Carbonate/macroeco_phylo/scripts/run_svd_and_make_diversity_slm_dbd_simulation_phylo_dict_carbonate.py 'human_gut_metagenome'
