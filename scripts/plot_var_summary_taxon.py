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
import math

import diversity_utils
import config



# load taxon dict
taxon_dict = dbd_utils.load_richness_diversity_prediction_taxon_dict()

# load phylo dict
phylo_dict = dbd_utils.load_richness_diversity_prediction_phylo_dict()


example_environment = 'human gut metagenome'

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label




fig = plt.figure(figsize = (8, 12)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


# richness taxon
ax_rank_vs_richness_taxon = plt.subplot2grid((3, 2), (0, 0))
ax_richness_taxon = plt.subplot2grid((3, 2), (1, 0))
ax_richness_cov_taxon = plt.subplot2grid((3, 2), (2, 0))


# diversity taxon
ax_rank_vs_diversity_taxon = plt.subplot2grid((3, 2), (0, 1))
ax_diversity_taxon = plt.subplot2grid((3, 2), (1, 1))
ax_diversity_cov_taxon = plt.subplot2grid((3, 2), (2, 1))


ax_rank_vs_richness_taxon.set_title("Richness predictions", fontsize=11, fontweight='bold')
ax_rank_vs_diversity_taxon.set_title("Diversity predictions", fontsize=11, fontweight='bold')



ax_rank_vs_richness_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_rank_vs_richness_taxon.transAxes)
ax_rank_vs_diversity_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_rank_vs_diversity_taxon.transAxes)

ax_richness_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[2], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_richness_taxon.transAxes)
ax_diversity_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[3], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_diversity_taxon.transAxes)

ax_richness_cov_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[4], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_richness_cov_taxon.transAxes)
ax_diversity_cov_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[5], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_diversity_cov_taxon.transAxes)


# make example plots
mean_richness_observed_example_taxon = []
mean_richness_predicted_example_taxon = []
mean_richness_predicted_example_corr_taxon = []
mean_diversity_observed_example_taxon = []
mean_diversity_predicted_example_taxon = []
mean_diversity_predicted_example_corr_taxon = []
for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

    if rank in taxon_dict[example_environment]['taxon']:

        mean_richness_observed = taxon_dict[example_environment]['taxon'][rank]['var_richness_observed']
        mean_richness_observed_example_taxon.append(mean_richness_observed)

        mean_richness_predicted = taxon_dict[example_environment]['taxon'][rank]['var_richness_predicted']
        mean_richness_predicted_example_taxon.append(mean_richness_predicted)

        mean_richness_predicted_corr = taxon_dict[example_environment]['taxon'][rank]['var_richness_predicted_plus_covariance']
        mean_richness_predicted_example_corr_taxon.append(mean_richness_predicted_corr)


        mean_diversity_observed = taxon_dict[example_environment]['taxon'][rank]['var_diversity_observed']
        mean_diversity_observed_example_taxon.append(mean_diversity_observed)

        mean_diversity_predicted = taxon_dict[example_environment]['taxon'][rank]['var_diversity_predicted']
        mean_diversity_predicted_example_taxon.append(mean_diversity_predicted)

        mean_diversity_predicted_corr = taxon_dict[example_environment]['taxon'][rank]['var_diversity_predicted_plus_covariance']
        mean_diversity_predicted_example_corr_taxon.append(mean_diversity_predicted_corr)


richness_taxon_all = []
richness_phylo_all = []
diversity_taxon_all = []
diversity_phylo_all = []

richness_cov_taxon_all = []
richness_cov_phylo_all = []
diversity_cov_taxon_all = []
diversity_cov_phylo_all = []
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    if environment == 'soil metagenome':
        continue

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)
    
    for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

        if rank in taxon_dict[environment]['taxon']:

            color_taxon = color_environment_taxon[rank_idx+1]

            mean_richness_predicted = taxon_dict[environment]['taxon'][rank]['var_richness_predicted']
            mean_richness_cov_predicted = taxon_dict[environment]['taxon'][rank]['var_richness_predicted_plus_covariance']
            mean_richness_observed = taxon_dict[environment]['taxon'][rank]['var_richness_observed']

            mean_diversity_predicted = taxon_dict[environment]['taxon'][rank]['var_diversity_predicted']
            mean_diversity_cov_predicted = taxon_dict[environment]['taxon'][rank]['var_diversity_predicted_plus_covariance']
            mean_diversity_observed = taxon_dict[environment]['taxon'][rank]['var_diversity_observed']

            color_taxon = color_environment_taxon[rank_idx+1]

            ax_richness_taxon.scatter(mean_richness_observed, mean_richness_predicted, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)
            ax_diversity_taxon.scatter(mean_diversity_observed, mean_diversity_predicted, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)
             
            ax_richness_cov_taxon.scatter(mean_richness_observed, mean_richness_cov_predicted, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)
           

            if (math.isnan(mean_diversity_observed) == False) and (math.isnan(mean_diversity_cov_predicted) == False):       
                ax_diversity_cov_taxon.scatter(mean_diversity_observed, mean_diversity_cov_predicted, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)
                diversity_cov_taxon_all.extend([mean_diversity_cov_predicted, mean_diversity_observed])

            
            richness_taxon_all.extend([mean_richness_predicted, mean_richness_observed])
            diversity_taxon_all.extend([mean_diversity_predicted, mean_diversity_observed])
            richness_cov_taxon_all.extend([mean_richness_cov_predicted, mean_richness_observed])



# plot example richness taxon
color_taxon_example = numpy.asarray([plot_utils.get_custom_cmap_taxon(example_environment)[i+1] for i in range(len(idx_taxa_ranks))])
ax_rank_vs_richness_taxon.plot(idx_taxa_ranks, mean_richness_predicted_example_taxon, ls='-', lw=2, alpha=0.8, c='k', label='Predicted')
ax_rank_vs_richness_taxon.scatter(idx_taxa_ranks, mean_richness_predicted_example_taxon, color='k', s=45, linewidth=0.8, edgecolors='k', zorder=3)
ax_rank_vs_richness_taxon.scatter(idx_taxa_ranks, mean_richness_observed_example_taxon, color=color_taxon_example, s=45, linewidth=0.8, edgecolors='k', zorder=4, label='Observed')
# plot cor
ax_rank_vs_richness_taxon.plot(idx_taxa_ranks, mean_richness_predicted_example_corr_taxon, ls=':', lw=2, alpha=0.8, c='k', zorder=2, label='Predicted + cov.')
ax_rank_vs_richness_taxon.scatter(idx_taxa_ranks, mean_richness_predicted_example_corr_taxon, color='k', s=45, linewidth=0.8, edgecolors='k', zorder=3)


# plot example diversity taxon
ax_rank_vs_diversity_taxon.plot(idx_taxa_ranks, mean_diversity_predicted_example_taxon, ls='-', lw=2, alpha=0.8, c='k')
ax_rank_vs_diversity_taxon.scatter(idx_taxa_ranks, mean_diversity_predicted_example_taxon, color='k', s=45, linewidth=0.8, edgecolors='k', zorder=3)
ax_rank_vs_diversity_taxon.scatter(idx_taxa_ranks, mean_diversity_observed_example_taxon, color=color_taxon_example, s=45, linewidth=0.8, edgecolors='k', zorder=4)
ax_rank_vs_diversity_taxon.plot(idx_taxa_ranks, mean_diversity_predicted_example_corr_taxon, ls=':', lw=2, alpha=0.8, c='k', zorder=2)
ax_rank_vs_diversity_taxon.scatter(idx_taxa_ranks, mean_diversity_predicted_example_corr_taxon, color='k', s=45, linewidth=0.8, edgecolors='k', facecolors='white', zorder=3)


# format subplots 
ax_rank_vs_richness_taxon.set_xticks(idx_taxa_ranks)
ax_rank_vs_richness_taxon.set_xticklabels(taxa_ranks_label, fontsize=10)
ax_rank_vs_richness_taxon.set_ylabel('Variance of richness', fontsize=12)
ax_rank_vs_richness_taxon.legend(loc='upper right')
ax_rank_vs_richness_taxon.set_yscale('log', base=10)


ax_rank_vs_diversity_taxon.set_xticks(idx_taxa_ranks)
ax_rank_vs_diversity_taxon.set_xticklabels(taxa_ranks_label, fontsize=10)
ax_rank_vs_diversity_taxon.set_ylabel('Variance of diversity', fontsize=12)



########
richness_taxon_min = min(richness_taxon_all)
richness_taxon_max = max(richness_taxon_all)
ax_richness_taxon.plot([richness_taxon_min*0.8,richness_taxon_max*1.2],[richness_taxon_min*0.8,richness_taxon_max*1.2], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_richness_taxon.set_xlim([richness_taxon_min*0.8,richness_taxon_max*1.2])
ax_richness_taxon.set_ylim([richness_taxon_min*0.8,richness_taxon_max*1.2])
ax_richness_taxon.set_xscale('log', base=10)
ax_richness_taxon.set_yscale('log', base=10)
ax_richness_taxon.set_xlabel('Observed variance of richness', fontsize=12)
ax_richness_taxon.set_ylabel('Predicted variance of richness', fontsize=12)
ax_richness_taxon.legend(loc='upper left')
ax_richness_taxon.legend(handles=plot_utils.legend_elements, loc='upper left', fontsize=7)



richness_cov_taxon_min = min(richness_cov_taxon_all)
richness_cov_taxon_max = max(richness_cov_taxon_all)
ax_richness_cov_taxon.plot([richness_cov_taxon_min*0.8,richness_cov_taxon_max*1.2],[richness_cov_taxon_min*0.8,richness_cov_taxon_max*1.2], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_richness_cov_taxon.set_xlim([richness_cov_taxon_min*0.8,richness_cov_taxon_max*1.2])
ax_richness_cov_taxon.set_ylim([richness_cov_taxon_min*0.8,richness_cov_taxon_max*1.2])
ax_richness_cov_taxon.set_xscale('log', base=10)
ax_richness_cov_taxon.set_yscale('log', base=10)
ax_richness_cov_taxon.set_xlabel('Observed variance of richness', fontsize=12)
ax_richness_cov_taxon.set_ylabel('Predicted variance of\nrichness + covariance', fontsize=12)


diversity_taxon_min = min(diversity_taxon_all)
diversity_taxon_max = max(diversity_taxon_all)
ax_diversity_taxon.plot([diversity_taxon_min*0.8,diversity_taxon_max*1.2],[diversity_taxon_min*0.8,diversity_taxon_max*1.2], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_diversity_taxon.set_xlim([diversity_taxon_min*0.8,diversity_taxon_max*1.2])
ax_diversity_taxon.set_ylim([diversity_taxon_min*0.8,diversity_taxon_max*1.2])
ax_diversity_taxon.set_xlabel('Observed variance of diversity', fontsize=12)
ax_diversity_taxon.set_ylabel('Predicted variance of diversity', fontsize=12)
ax_diversity_taxon.set_xscale('log', base=10)
ax_diversity_taxon.set_yscale('log', base=10)


diversity_cov_taxon_min = min(diversity_cov_taxon_all)
diversity_cov_taxon_max = max(diversity_cov_taxon_all)
ax_diversity_cov_taxon.plot([diversity_cov_taxon_min*0.8,diversity_cov_taxon_max*1.2],[diversity_cov_taxon_min*0.8,diversity_cov_taxon_max*1.2], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_diversity_cov_taxon.set_xlim([diversity_cov_taxon_min*0.8,diversity_cov_taxon_max*1.2])
ax_diversity_cov_taxon.set_ylim([diversity_cov_taxon_min*0.8,diversity_cov_taxon_max*1.2])
ax_diversity_cov_taxon.set_xlabel('Observed variance of diversity', fontsize=12)
ax_diversity_cov_taxon.set_ylabel('Predicted variance of\ndiversity + covariance', fontsize=12)
ax_diversity_cov_taxon.set_xscale('log', base=10)
ax_diversity_cov_taxon.set_yscale('log', base=10)




fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%svar_summary_taxon.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
