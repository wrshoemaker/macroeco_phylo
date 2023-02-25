#!/bin/bash


# 0.01,  0.010500551788704102

declare -a distances=("freshwater_metagenome" "marine_sediment_metagenome" "marine_metagenome" "soil_metagenome" "human_oral_metagenome" "human_skin_metagenome" "microbial_mat_metagenome" "freshwater_sediment_metagenome" "human_gut_metagenome")



#"0.09006280202112786" "0.09457091168586576" "0.09930467558623954" "0.10427538888537682" "0.10949491212781594" "0.11497569953977363" "0.12073082874598748" "0.12677403197404075" "0.13311972882062445" "0.13978306065792132" "0.14677992676220697" "0.15412702225087496" "0.1618418779184062" "0.16994290206633514" "0.17844942442702214" "0.1873817422860384" "0.19676116891321516" "0.20661008441791714" "0.21695198914988656" "0.22781155977307543" "0.23921470814626386" "0.25118864315095807" #"0.26376193561409494" "0.27696458648046407" "0.2908280983975128" "0.30538555088334157" "0.320671679257246" "0.33672295752114223" "0.3535776853896366" "0.3712760796764005" "0.3898603702549074" "0.4093749008225011" "0.4298662347082279" "0.45138326597689776")

#declare -a distances=("0.010500551788704102")

counter=0

for distance in "${distances[@]}"
do

  bash_out="/N/u/wrshoema/Carbonate/macroeco_phylo/scripts/bash/qsub_make_taxon_dict_carbonate_${counter}.sh"
  if [ -f $bash_out ]; then
      rm $bash_out
  fi

  echo '#!/bin/bash' >> $bash_out
  echo '#SBATCH --mail-user=williamrshoemaker@gmail.com' >> $bash_out
  echo '#SBATCH --nodes=1' >> $bash_out
  echo '#SBATCH --ntasks-per-node=1' >> $bash_out
  echo '#SBATCH --cpus-per-task=8' >> $bash_out
  echo '#SBATCH --time=72:00:00' >> $bash_out
  echo '#SBATCH --mail-type=FAIL,BEGIN,END' >> $bash_out
  echo "#SBATCH --job-name=rank-${distance}" >> $bash_out
  echo '#SBATCH --mem=15G' >> $bash_out

  echo 'module unload python' >> $bash_out
  echo 'module load anaconda' >> $bash_out
  echo 'conda activate macroeco_phylo' >> $bash_out

  echo "python /N/u/wrshoema/Carbonate/macroeco_phylo/scripts/make_taxon_dict_carbonate.py ${distance}" >> $bash_out

  sbatch ${bash_out}

  let counter++

done




#python /N/u/wrshoema/Carbonate/macroeco_phylo/scripts/make_taxon_dict_carbonate.py "human gut metagenome"
