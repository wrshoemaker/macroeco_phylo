#!/bin/bash


# 0.01,  0.010500551788704102

declare -a distances=("freshwater_metagenome" "marine_sediment_metagenome" "marine_metagenome" "soil_metagenome" "human_oral_metagenome" "human_skin_metagenome" "microbial_mat_metagenome" "freshwater_sediment_metagenome" "human_gut_metagenome")


counter=0

for distance in "${distances[@]}"
do

  python /Users/williamrshoemaker/GitHub/macroeco_phylo/scripts/make_diversity_vs_diversity_taxon_carbonate.py ${distance}

done
