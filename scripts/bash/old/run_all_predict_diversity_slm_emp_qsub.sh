#!/bin/bash



declare -a environments=("freshwater_metagenome" "marine_sediment_metagenome" "marine_metagenome" "soil_metagenome" "human_oral_metagenome" "human_skin_metagenome" "microbial_mat_metagenome" "freshwater_sediment_metagenome" "human_gut_metagenome")


for environment in "${environments[@]}"
do

  bash_out="/N/u/wrshoema/Carbonate/macroeco_phylo/scripts/bash/qsub_predict_diversity_slm_emp_carbonate_${environment}.sh"
  bash_error="/N/u/wrshoema/Carbonate/macroeco_phylo/scripts/bash/qsub_predict_diversity_slm_emp_carbonate_${environment}.error"

  if [ -f $bash_out ]; then
      rm $bash_out
  fi

  echo '#!/bin/bash' >> $bash_out
  echo '#SBATCH --mail-user=williamrshoemaker@gmail.com' >> $bash_out
  echo '#SBATCH --nodes=1' >> $bash_out
  echo '#SBATCH --ntasks-per-node=1' >> $bash_out
  echo '#SBATCH --cpus-per-task=1' >> $bash_out
  echo '#SBATCH --time=48:00:00' >> $bash_out
  echo '#SBATCH --mail-type=FAIL,BEGIN,END' >> $bash_out
  echo "#SBATCH --job-name=div-${environment}" >> $bash_out
  echo '#SBATCH --mem=32G' >> $bash_out
  echo "#SBATCH --error=${bash_error}" >> $bash_out

  echo 'module unload python' >> $bash_out
  echo 'module load anaconda' >> $bash_out
  echo 'conda activate macroeco_phylo' >> $bash_out

  echo "python /N/u/wrshoema/Carbonate/macroeco_phylo/scripts/predict_diversity_slm_emp_carbonate.py ${environment}" >> $bash_out

  #sbatch ${bash_out}

done


#python /N/u/wrshoema/Carbonate/macroeco_phylo/scripts/make_taxon_dict_carbonate.py "human gut metagenome" "richness"
