#!/bin/bash

muscle -in /Users/williamrshoemaker/GitHub/macroeco_phylo/data/dalbello/seqtab-nochim.fna -out /Users/williamrshoemaker/GitHub/macroeco_phylo/data/dalbello/seqtab-nochim_muscle.fna
muscle -in /Users/williamrshoemaker/GitHub/macroeco_phylo/data/dalbello/seqtab-nochim_with_outgroup.fna -out /Users/williamrshoemaker/GitHub/macroeco_phylo/data/dalbello/seqtab-nochim_with_outgroup_muscle.fna


asv_muscle=/Users/williamrshoemaker/GitHub/macroeco_phylo/data/dalbello/seqtab-nochim_muscle.fna
asv_muscle_with_outgroup=/Users/williamrshoemaker/GitHub/macroeco_phylo/data/dalbello/seqtab-nochim_with_outgroup_muscle.fna

#raxml-ng --all --msa ${asv_muscle} --msa-format FASTA --data-type DNA --seed 123456789 --model GTR+G --bs-trees autoMRE

#~/raxml-ng_v1.1.0_macos_x86_64/raxml-ng --redo
~/raxml-ng_v1.1.0_macos_x86_64/raxml-ng --redo  --all --msa ${asv_muscle} --msa-format FASTA --data-type DNA --seed 123456789 --model GTR+G --bs-trees autoMRE
#~/raxml-ng_v1.1.0_macos_x86_64/raxml-ng --redo  --all --msa ${asv_muscle_with_outgroup} --msa-format FASTA --data-type DNA --seed 123456789 --model GTR+G --bs-trees autoMRE -outgroup NC_005042_1_353331_354795_Prochlorococcus_marinus_subsp_marinus_str_CCMP1375_complete_genome
