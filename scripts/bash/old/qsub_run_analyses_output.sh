#!/bin/bash
#$ -N run_analyses
#$ -e /u/project/ngarud/wrshoema/macroeco_phylo/scripts/run_analyses_error
#$ -o /u/project/ngarud/wrshoema/macroeco_phylo/scripts/run_analyses_output
#$ -l h_data=8G
#$ -l time=48:00:00
#$ -l highp
#$ -m bea

. /u/local/Modules/default/init/modules.sh

module unload python
module load python/3.6.1

module load anaconda

source activate macroeco_phylo

python3 /u/home/w/wrshoema/project-ngarud/macroeco_phylo/scripts/run_analyses.py
