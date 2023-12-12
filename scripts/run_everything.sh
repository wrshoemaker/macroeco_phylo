#!/bin/bash

cd ~/GitHub/macroeco_phylo

# create conda environment
#conda env create -f environment.yml

# activate your conda environment
source activate macroeco_phylo


# run analyses
#python run_analysis.py



# Generate main manuscript plots

# Fig. 2
python plot_figure-2-gamma_summary.py


# Fig. 3
# predict mean estimator
python plot_figure-3-mean_summary_phylo.py


# Fig. 4
# predict variance of estimator
python plot_figure-4-var_summary_phylo.py


# Fig. 5
# predict DBD slope
# richness_dbd_dict.pickle ==> figure-5-source-data.pickle
python plot_figure-5-dbd_summary_phylo.py

# Fig. 6
# predict DBD diversity slope using simulation
# diversity_dbd_simulation_taxon_dict.pickle ==> figure-6-source-data.pickle
python plot_figure-6-dbd_diversity_slm_simulation_mean_phylo.py




###############################
# Generate supplemental plots #
###############################


# fraction OTUs vs. distance
# Figure 1–figure supplement 2
python plot_figure-1-figure-supplement-2-coarse_grained_dist.py



# plot coarse-grained AFDs
# Figure 2–figure supplement 1
python plot_figure-2-figure-supplement-1-coarse_afd_taxon.py
# Figure 2–figure supplement 2
python plot_figure-2-figure-supplement-2-coarse_afd_phylo.py


# precict gamma occupancy
# Figure 2–figure supplement 3
python plot_figure-2-figure-supplement-3-gamma_occupancy_all_taxon.py
# Figure 2–figure supplement 4
python plot_figure-2-figure-supplement-4-gamma_occupancy_all_phylo.py


# MAD vs. occupancy
# Figure 2–figure supplement 5
python plot_figure-2-figure-supplement-5-abundance_vs_occupancy_all_taxon.py
# Figure 2–figure supplement 6
python plot_figure-2-figure-supplement-6-abundance_vs_occupancy_all_phylo.py


# predict variance of occupancy
# Figure 2–figure supplement 7
python plot_figure-2-figure-supplement-7-gamma_occupancy_var_all_taxon.py
# Figure 2–figure supplement 8
python plot_figure-2-figure-supplement-8-gamma_occupancy_var_all_phylo.py


# plot gamma occupancy error across ranks
# Figure 2–figure supplement 9
python plot_figure-2-figure-supplement-9-gamma_error_summary.py


# plot fine vs. coarse variance
# Figure 2–figure supplement 10
python plot_figure-2-figure-supplement-10-fine_vs_coarse_variance_taxon.py
# Figure 2–figure supplement 11
python plot_figure-2-figure-supplement-11-fine_vs_coarse_variance_phylo.py

# variance ratio
# Figure 2–figure supplement 12
python plot_figure-2-figure-supplement-12-fine_vs_coarse_variance_ratio_taxon.py
# Figure 2–figure supplement 13
python plot_figure-2-figure-supplement-13-fine_vs_coarse_variance_ratio_phylo.py



# plot proof that analytic mean and var for richness and diversity matches simulations
# one plot
# Figure 3–figure supplement 1
python plot_figure-3-figure-supplement-1-prediction_vs_simulation.py

# mean taxon
# Figure 3–figure supplement 2
# richness_diversity_prediction_taxon_dict.pickle ==> figure-3-figure-supplement-2-source-data.pickle
python plot_figure-3-figure-supplement-2-mean_summary_taxon.py

# UNTB prediction
# Figure 3–figure supplement 3
# richness_dbd_neutral_dict.pickle ==> figure-3-figure-supplement-3-source-data.pickle
python plot_figure-3-figure-supplement-3-richness_neutral_mean.py



# var taxon
# Figure 4–figure supplement 1
python plot_figure-4-figure-supplement-1-var_summary_taxon.py

# DBD summary taxon
# Figure 5–figure supplement 2
python plot_figure-5-figure-supplement-2-dbd_summary_taxon.py

# plot DBD richness slope

# Figure 5–figure supplement 3
python plot_figure-5-figure-supplement-3-dbd_slope_slm_richness_taxon_all.py
# Figure 5–figure supplement 4
python plot_figure-5-figure-supplement-4-dbd_slope_slm_richness_phylo_all.py

# UNTB
# Figure 5–figure supplement 5
python plot_figure-5-figure-supplement-5-dbd_slope_neutral_richness_taxon_all.py
# Figure 5–figure supplement 6
python plot_figure-5-figure-supplement-6-dbd_slope_neutral_richness_phylo_all.py
# Figure 5–figure supplement 7
python plot_figure-5-figure-supplement-7-dbd_richness_neutral_mean.py
# Figure 5–figure supplement 8
python plot_figure-5-figure-supplement-8-gamma_error_summary_neutral_taxon.py
# Figure 5–figure supplement 9
python plot_figure-5-figure-supplement-9-gamma_error_summary_neutral_phylo.py


# plot DBD diversity slope
# Figure 5–figure supplement 10
python plot_figure-5-figure-supplement-10-dbd_slope_slm_diversity_taxon_all.py
# Figure 5–figure supplement 11
python plot_figure-5-figure-supplement-11-dbd_slope_slm_diversity_phylo_all.py


# plot DBD diversity slope simulation with correlation
# Figure 6–figure supplement 1
python plot_figure-6-figure-supplement-1-dbd_slope_slm_simulation_diversity_taxon_all.py
# Figure 6–figure supplement 2
python plot_figure-6-figure-supplement-2-dbd_slope_slm_simulation_diversity_phylo_all.py
# Figure 6–figure supplement 3
python plot_figure-6-figure-supplement-3-dbd_diversity_slm_simulation_mean_taxon.py









