#!/bin/bash



declare -a environments=("freshwater_metagenome" "marine_sediment_metagenome" "marine_metagenome" "soil_metagenome" "human_oral_metagenome" "human_skin_metagenome" "microbial_mat_metagenome" "freshwater_sediment_metagenome" "human_gut_metagenome")
#declare -a environments=("freshwater metagenome" "marine sediment metagenome" "marine metagenome" "soil metagenome" "human oral metagenome" "human skin metagenome" "microbial mat metagenome" "freshwater sediment metagenome" "human gut lmetagenome")

#declare -a measures=("diversity")
#declare -a distances=("soil_metagenome" "human_skin_metagenome")


counter=0


for environment in "${environments[@]}"
do

    bash_out="/N/u/wrshoema/Carbonate/macroeco_phylo/scripts/bash/qsub_make_cov_phylo_dict_carbonate_${counter}.sh"
    if [ -f $bash_out ]; then
        rm $bash_out
    fi

    echo '#!/bin/bash' >> $bash_out
    echo '#SBATCH --mail-user=williamrshoemaker@gmail.com' >> $bash_out
    echo '#SBATCH --nodes=1' >> $bash_out
    echo '#SBATCH --ntasks-per-node=1' >> $bash_out
    echo '#SBATCH --cpus-per-task=8' >> $bash_out
    echo '#SBATCH --time=96:00:00' >> $bash_out
    echo '#SBATCH --mail-type=FAIL,BEGIN,END' >> $bash_out
    echo "#SBATCH --job-name=int-${environment}" >> $bash_out
    echo '#SBATCH --mem=20G' >> $bash_out

    echo 'module unload python' >> $bash_out
    echo 'module load anaconda' >> $bash_out
    echo 'conda activate macroeco_phylo' >> $bash_out

    echo "python /N/u/wrshoema/Carbonate/macroeco_phylo/scripts/run_make_cov_phylo_dict_carbonate.py ${environment}" >> $bash_out

    sbatch ${bash_out}

    let counter++
done



#python /N/u/wrshoema/Carbonate/macroeco_phylo/scripts/make_taxon_dict_carbonate.py "human gut metagenome" "richness"
