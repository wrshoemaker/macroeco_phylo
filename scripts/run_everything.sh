#!/bin/bash

cd ~/GitHub/macroeco_phylo

# create conda environment
conda env create -f environment.yml

# activate your conda environment
source activate macroeco_phylo


# run analyses
python run_analysis.py



# Generate main manuscript plots
# Fig. 2
python plot_coarse_grained_dist.py


# Fig. 3
python plot_gamma_example.py


# Fig. 4
# predict mean estimator
python plot_mean_summary.py


# Fig. 5
# predict variance of estimator
python plot_var_summary.py


# Fig. 6
# predict DBD slope
python plot_dbd_summary.py

# Fig. 7
# predict DBD diversity slope using simulation
python plot_dbd_diversity_slm_simulation.py


###############################
# Generate supplemental plots #
###############################

# plot coarse-grained AFDs
python plot_coarse_afd_taxa.py
python plot_coarse_afd_phylo.py


# precict gamma occupancy
python plot_gamma_occupancy_all_taxon.py
python plot_gamma_occupancy_all_phylo.py


# plot gamma occupancy error across ranks
python plot_gamma_error_summary.py


# MAD vs. occupancy
python plot_abundance_vs_occupancy_all_phylo.py
python plot_abundance_vs_occupancy_all_taxon.py


# plot fine vs. coarse variance
python plot_fine_vs_coarse_variance_taxon.py
python plot_fine_vs_coarse_variance_phylo.py


# plot proof that analytic mean and var for richness and diversity matches simulations
# one plot
python plot_prediction_vs_simulation.py


# plot DBD richness slope
python plot_dbd_richness_slm_taxon_all.py
python plot_dbd_richness_slm_phylo_all.py


# plot DBD diversity slope
python plot_dbd_diversity_slm_taxon_all.py
python plot_dbd_diversity_slm_phylo_all.py


# plot DBD diversity slope simulation with correlation
python plot_dbd_diversity_slm_simulation_taxon_all.py
python plot_dbd_diversity_slm_simulation_phylo_all.py






# things to consider plotting
# plot predict mean richness all
# plot predict var richness all 
# plot predict mean diversity all
# plot predict var diversity all





