import os
#from Bio import Phylo
import random
import copy
import sys
import numpy
import random
import pickle
import scipy.stats as stats
import matplotlib.pyplot as plt
import functools
import operator

import dbd_utils
import plot_utils

import diversity_utils
import config





richness_dict = dbd_utils.load_richness_slm_dbd_dict()
diversity_dict = dbd_utils.load_diversity_slm_dbd_dict()


fig = plt.figure(figsize = (16, 8))
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

# (#rows, # columns)
ax_richness_dist_taxon = plt.subplot2grid((2, 4), (0, 0))
ax_richness_dist_phylo = plt.subplot2grid((2, 4), (0, 1))
ax_diversity_dist_taxon = plt.subplot2grid((2, 4), (0, 2))
ax_diversity_dist_phylo = plt.subplot2grid((2, 4), (0, 3))

ax_richness_taxon = plt.subplot2grid((2, 4), (1, 0))
ax_richness_phylo = plt.subplot2grid((2, 4), (1, 1))
ax_diversity_taxon = plt.subplot2grid((2, 4), (1, 2))
ax_diversity_phylo = plt.subplot2grid((2, 4), (1, 3))


ax_richness_dist_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_richness_dist_taxon.transAxes)
ax_richness_dist_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_richness_dist_phylo.transAxes)
ax_diversity_dist_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[2], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_diversity_dist_taxon.transAxes)
ax_diversity_dist_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[3], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_diversity_dist_phylo.transAxes)

ax_richness_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[4], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_richness_taxon.transAxes)
ax_richness_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[5], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_richness_phylo.transAxes)
ax_diversity_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[6], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_diversity_taxon.transAxes)
ax_diversity_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[7], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_diversity_phylo.transAxes)



ax_richness_dist_taxon.set_title("Taxonomic coarse-graining", fontsize=11)
ax_richness_dist_phylo.set_title("Phylogenetic coarse-graining", fontsize=11)
ax_diversity_dist_taxon.set_title("Taxonomic coarse-graining", fontsize=11)
ax_diversity_dist_phylo.set_title("Phylogenetic coarse-graining", fontsize=11)




observed_taxon_all_environments = []
predicted_taxon_all_environments = []
observed_phylo_all_environments = []
predicted_phylo_all_environments = []

richness_slope_taxon_all = []
richness_slope_phylo_all = []
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    distances  = list(richness_dict[environment]['phylo'].keys())
    distances.sort()

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)
    color_environment_phylo = plot_utils.get_custom_cmap_phylo(environment, len(distances))

    taxa_ranks = list(richness_dict[environment]['taxon'].keys())

    observed_taxon_all = []
    predicted_taxon_all = []
    for coarse_rank_idx, coarse_rank in enumerate(taxa_ranks):

        if 'slope_all' in richness_dict[environment]['taxon'][coarse_rank]:
        
            observed = richness_dict[environment]['taxon'][coarse_rank]['slope_all']
            observed_taxon_all.extend(observed)

            predicted = richness_dict[environment]['taxon'][coarse_rank]['slope_slm_all']
            predicted_taxon_all.extend(predicted)

            mean_slope = richness_dict[environment]['taxon'][coarse_rank]['mean_slope']
            mean_slope_slm = richness_dict[environment]['taxon'][coarse_rank]['mean_slope_slm']

            color_taxon = color_environment_taxon[coarse_rank_idx+1]
            ax_richness_taxon.scatter(mean_slope, mean_slope_slm, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)
            richness_slope_taxon_all.extend([mean_slope, mean_slope_slm])



    observed_phylo_all = []
    predicted_phylo_all = []
    for distance_idx, distance in enumerate(distances):

        if 'slope_all' in richness_dict[environment]['phylo'][distance]:
        
            observed = richness_dict[environment]['phylo'][distance]['slope_all']
            observed_phylo_all.extend(observed)

            predicted = richness_dict[environment]['phylo'][distance]['slope_slm_all']
            predicted_phylo_all.extend(predicted)

            mean_slope = richness_dict[environment]['phylo'][distance]['mean_slope']
            mean_slope_slm = richness_dict[environment]['phylo'][distance]['mean_slope_slm']

            color_phylo = color_environment_phylo[distance_idx]
            ax_richness_phylo.scatter(mean_slope, mean_slope_slm, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)
            richness_slope_phylo_all.extend([mean_slope, mean_slope_slm])


    hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(observed_taxon_all)
    ax_richness_dist_taxon.scatter(bins_mean_to_plot, hist_to_plot, color=color_environment_taxon[-2], alpha=0.9, s=30, linewidth=0.8, edgecolors='k', label=diversity_utils.format_environment_label(environment))

    hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(observed_phylo_all)
    ax_richness_dist_phylo.scatter(bins_mean_to_plot, hist_to_plot, color=color_environment_taxon[-2], alpha=0.9, s=30, linewidth=0.8, edgecolors='k', label=diversity_utils.format_environment_label(environment))


    observed_taxon_all_environments.extend(observed_taxon_all)
    predicted_taxon_all_environments.extend(predicted_taxon_all)
    observed_phylo_all_environments.extend(observed_phylo_all)
    predicted_phylo_all_environments.extend(predicted_phylo_all)




rho2_taxon = numpy.corrcoef(observed_taxon_all_environments, predicted_taxon_all_environments)[0,1]**2
rho2_phylo = numpy.corrcoef(observed_phylo_all_environments, predicted_phylo_all_environments)[0,1]**2

ax_richness_taxon.text(0.2, 0.9, r'$\rho^{2} = $' + str(round(rho2_taxon, 3)), fontsize=10, ha='center', va='center', transform=ax_richness_taxon.transAxes)
ax_richness_phylo.text(0.2, 0.9, r'$\rho^{2} = $' + str(round(rho2_phylo, 3)), fontsize=10, ha='center', va='center', transform=ax_richness_phylo.transAxes)




observed_taxon_all_environments = []
predicted_taxon_all_environments = []
observed_phylo_all_environments = []
predicted_phylo_all_environments = []
# diversity 
diversity_slope_taxon_all = []
diversity_slope_phylo_all = []
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    distances  = list(diversity_dict[environment]['phylo'].keys())
    distances.sort()

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)
    color_environment_phylo = plot_utils.get_custom_cmap_phylo(environment, len(distances))

    taxa_ranks = list(diversity_dict[environment]['taxon'].keys())

    observed_taxon_all = []
    predicted_taxon_all = []
    for coarse_rank_idx, coarse_rank in enumerate(taxa_ranks):

        if 'slope_all' in diversity_dict[environment]['taxon'][coarse_rank]:
        
            observed = diversity_dict[environment]['taxon'][coarse_rank]['slope_all']
            observed_taxon_all.extend(observed)

            predicted = diversity_dict[environment]['taxon'][coarse_rank]['slope_slm_all']
            predicted_taxon_all.extend(predicted)

            mean_slope = diversity_dict[environment]['taxon'][coarse_rank]['mean_slope']
            mean_slope_slm = diversity_dict[environment]['taxon'][coarse_rank]['mean_slope_slm']

            color_taxon = color_environment_taxon[coarse_rank_idx+1]
            ax_diversity_taxon.scatter(mean_slope, mean_slope_slm, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)
            diversity_slope_taxon_all.extend([mean_slope, mean_slope_slm])


    observed_phylo_all = []
    predicted_phylo_all = []
    for distance_idx, distance in enumerate(distances):

        if 'slope_all' in diversity_dict[environment]['phylo'][distance]:
        
            observed = diversity_dict[environment]['phylo'][distance]['slope_all']
            observed_phylo_all.extend(observed)

            predicted = diversity_dict[environment]['phylo'][distance]['slope_slm_all']
            predicted_phylo_all.extend(predicted)
            
            mean_slope = diversity_dict[environment]['phylo'][distance]['mean_slope']
            mean_slope_slm = diversity_dict[environment]['phylo'][distance]['mean_slope_slm']

            color_phylo = color_environment_phylo[distance_idx]
            ax_diversity_phylo.scatter(mean_slope, mean_slope_slm, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)
            diversity_slope_phylo_all.extend([mean_slope, mean_slope_slm])


    hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(observed_taxon_all)
    ax_diversity_dist_taxon.scatter(bins_mean_to_plot, hist_to_plot, color=color_environment_taxon[-2], alpha=0.9, s=30, linewidth=0.8, edgecolors='k', label=diversity_utils.format_environment_label(environment))

    hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(observed_phylo_all)
    ax_diversity_dist_phylo.scatter(bins_mean_to_plot, hist_to_plot, color=color_environment_taxon[-2], alpha=0.9, s=30, linewidth=0.8, edgecolors='k', label=diversity_utils.format_environment_label(environment))


    observed_taxon_all_environments.extend(observed_taxon_all)
    predicted_taxon_all_environments.extend(predicted_taxon_all)
    observed_phylo_all_environments.extend(observed_phylo_all)
    predicted_phylo_all_environments.extend(predicted_phylo_all)



# plot R2
rho2_taxon = numpy.corrcoef(observed_taxon_all_environments, predicted_taxon_all_environments)[0,1]**2
rho2_phylo = numpy.corrcoef(observed_phylo_all_environments, predicted_phylo_all_environments)[0,1]**2

ax_diversity_taxon.text(0.2, 0.9, r'$\rho^{2} = $' + str(round(rho2_taxon, 3)), fontsize=10, ha='center', va='center', transform=ax_diversity_taxon.transAxes)
ax_diversity_phylo.text(0.2, 0.9, r'$\rho^{2} = $' + str(round(rho2_phylo, 3)), fontsize=10, ha='center', va='center', transform=ax_diversity_phylo.transAxes)



#ax_richness_dist_taxon.set_xscale('log', base=10)
ax_richness_dist_taxon.set_yscale('log', base=10)
ax_richness_dist_taxon.set_xlabel("Observed intra vs. inter-group slope, richness", fontsize = 10)
ax_richness_dist_taxon.set_ylabel("Probability density", fontsize = 10)

ax_richness_dist_phylo.set_yscale('log', base=10)
ax_richness_dist_phylo.set_xlabel("Observed intra vs. inter-group slope, richness", fontsize = 10)
ax_richness_dist_phylo.set_ylabel("Probability density", fontsize = 10)

ax_diversity_dist_taxon.set_yscale('log', base=10)
ax_diversity_dist_taxon.set_xlabel("Observed intra vs. inter-group slope, diversity", fontsize = 10)
ax_diversity_dist_taxon.set_ylabel("Probability density", fontsize = 10)

ax_diversity_dist_phylo.set_yscale('log', base=10)
ax_diversity_dist_phylo.set_xlabel("Observed intra vs. inter-group slope, diversity", fontsize = 10)
ax_diversity_dist_phylo.set_ylabel("Probability density", fontsize = 10)



ax_richness_taxon.set_xscale('log', base=10)
ax_richness_taxon.set_yscale('log', base=10)
ax_richness_taxon.set_xlabel("Observed mean slope, richness", fontsize = 10)
ax_richness_taxon.set_ylabel("Predicted mean slope, richness", fontsize = 10)
richness_slope_taxon_min = min(richness_slope_taxon_all)
richness_slope_taxon_max = max(richness_slope_taxon_all)
ax_richness_taxon.plot([richness_slope_taxon_min*0.8, richness_slope_taxon_max*1.2], [richness_slope_taxon_min*0.8, richness_slope_taxon_max*1.2], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax_richness_taxon.set_xlim([richness_slope_taxon_min*0.8,richness_slope_taxon_max*1.2])
ax_richness_taxon.set_ylim([richness_slope_taxon_min*0.8,richness_slope_taxon_max*1.2])
ax_richness_taxon.legend(loc="lower right", fontsize=7)


ax_richness_phylo.set_xscale('log', base=10)
ax_richness_phylo.set_yscale('log', base=10)
ax_richness_phylo.set_xlabel("Observed mean slope, richness", fontsize = 10)
ax_richness_phylo.set_ylabel("Predicted mean slope, richness", fontsize = 10)
richness_slope_phylo_min = min(richness_slope_phylo_all)
richness_slope_phylo_max = max(richness_slope_phylo_all)
ax_richness_phylo.plot([richness_slope_phylo_min*0.8, richness_slope_phylo_max*1.2], [richness_slope_phylo_min*0.8, richness_slope_phylo_max*1.2], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax_richness_phylo.set_xlim([richness_slope_phylo_min*0.8,richness_slope_phylo_max*1.2])
ax_richness_phylo.set_ylim([richness_slope_phylo_min*0.8,richness_slope_phylo_max*1.2])




ax_diversity_taxon.set_xlabel("Observed mean slope, diversity", fontsize = 10)
ax_diversity_taxon.set_ylabel("Predicted mean slope, diversity", fontsize = 10)
diversity_slope_taxon_min = min(diversity_slope_taxon_all)
diversity_slope_taxon_max = max(diversity_slope_taxon_all)
ax_diversity_taxon.plot([diversity_slope_taxon_min, diversity_slope_taxon_max], [diversity_slope_taxon_min, diversity_slope_taxon_max], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax_diversity_taxon.set_xlim([diversity_slope_taxon_min,diversity_slope_taxon_max])
ax_diversity_taxon.set_ylim([diversity_slope_taxon_min,diversity_slope_taxon_max])



ax_diversity_phylo.set_xlabel("Observed mean slope, diversity", fontsize = 10)
ax_diversity_phylo.set_ylabel("Predicted mean slope, diversity", fontsize = 10)
diversity_slope_phylo_min = min(diversity_slope_phylo_all)
diversity_slope_phylo_max = max(diversity_slope_phylo_all)
ax_diversity_phylo.plot([diversity_slope_phylo_min, diversity_slope_phylo_max], [diversity_slope_phylo_min, diversity_slope_phylo_max], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax_diversity_phylo.set_xlim([diversity_slope_phylo_min,diversity_slope_phylo_max])
ax_diversity_phylo.set_ylim([diversity_slope_phylo_min,diversity_slope_phylo_max])






fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sdbd_summary.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
