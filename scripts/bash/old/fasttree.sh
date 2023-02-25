#!/bin/bash

asv_muscle=/Users/williamrshoemaker/GitHub/macroeco_phylo/data/dalbello/seqtab-nochim_muscle.fna
tree=/Users/williamrshoemaker/GitHub/macroeco_phylo/data/dalbello/seqtab-nochim_muscle_fasttree.tre

#-nt = nucleotide
# -gtr = general time reversible
#-gamma = distribution 

FastTree -nt -gtr -gamma ${asv_muscle} > ${tree}
