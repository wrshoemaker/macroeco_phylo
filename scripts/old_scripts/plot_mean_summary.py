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



# load taxon dict
taxon_dict = dbd_utils.load_richness_diversity_prediction_taxon_dict()

# load phylo dict
phylo_dict = dbd_utils.load_richness_diversity_prediction_phylo_dict()


example_environment = 'human gut metagenome'

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label




fig = plt.figure(figsize = (16, 8)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

# (#rows, # columns)
ax_rank_vs_richness_taxon = plt.subplot2grid((2, 4), (0, 0))
ax_rank_vs_richness_phylo = plt.subplot2grid((2, 4), (0, 1))
ax_rank_vs_diversity_taxon = plt.subplot2grid((2, 4), (0, 2))
ax_rank_vs_diversity_phylo = plt.subplot2grid((2, 4), (0, 3))

ax_richness_taxon = plt.subplot2grid((2, 4), (1, 0))
ax_richness_phylo = plt.subplot2grid((2, 4), (1, 1))
ax_diversity_taxon = plt.subplot2grid((2, 4), (1, 2))
ax_diversity_phylo = plt.subplot2grid((2, 4), (1, 3))


ax_rank_vs_richness_taxon.set_title("Taxonomic coarse-graining", fontsize=11)
ax_rank_vs_richness_phylo.set_title("Phylogenetic coarse-graining", fontsize=11)
ax_rank_vs_diversity_taxon.set_title("Taxonomic coarse-graining", fontsize=11)
ax_rank_vs_diversity_phylo.set_title("Phylogenetic coarse-graining", fontsize=11)


ax_rank_vs_richness_taxon.text(1.2, 1.2, "Richness predictions", fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_rank_vs_richness_taxon.transAxes)
ax_rank_vs_diversity_taxon.text(1.2, 1.2, "Diversity predictions", fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_rank_vs_diversity_taxon.transAxes)






ax_rank_vs_richness_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_rank_vs_richness_taxon.transAxes)
ax_rank_vs_richness_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_rank_vs_richness_phylo.transAxes)
ax_rank_vs_diversity_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[2], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_rank_vs_diversity_taxon.transAxes)
ax_rank_vs_diversity_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[3], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_rank_vs_diversity_phylo.transAxes)

ax_richness_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[4], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_richness_taxon.transAxes)
ax_richness_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[5], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_richness_phylo.transAxes)
ax_diversity_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[6], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_diversity_taxon.transAxes)
ax_diversity_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[7], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_diversity_phylo.transAxes)




# make example plots
mean_richness_observed_example_taxon = []
mean_richness_predicted_example_taxon = []
mean_diversity_observed_example_taxon = []
mean_diversity_predicted_example_taxon = []
for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

    if rank in taxon_dict[example_environment]['taxon']:

        mean_richness_observed = taxon_dict[example_environment]['taxon'][rank]['mean_richness_observed']
        mean_richness_observed_example_taxon.append(mean_richness_observed)

        mean_richness_predicted = taxon_dict[example_environment]['taxon'][rank]['mean_richness_predicted']
        mean_richness_predicted_example_taxon.append(mean_richness_predicted)

        mean_diversity_observed = taxon_dict[example_environment]['taxon'][rank]['mean_diversity_observed']
        mean_diversity_observed_example_taxon.append(mean_diversity_observed)

        mean_diversity_predicted = taxon_dict[example_environment]['taxon'][rank]['mean_diversity_predicted']
        mean_diversity_predicted_example_taxon.append(mean_diversity_predicted)



distance_example = list(phylo_dict[example_environment]['phylo'].keys())
distance_example.sort()
distance_example = numpy.asarray(distance_example)
mean_richness_observed_example_phylo = []
mean_richness_predicted_example_phylo = []
mean_diversity_observed_example_phylo = []
mean_diversity_predicted_example_phylo = []
for d_idx, d in enumerate(distance_example):

    mean_richness_observed = phylo_dict[example_environment]['phylo'][d]['mean_richness_observed']
    mean_richness_observed_example_phylo.append(mean_richness_observed)

    mean_richness_predicted = phylo_dict[example_environment]['phylo'][d]['mean_richness_predicted']
    mean_richness_predicted_example_phylo.append(mean_richness_predicted)

    mean_diversity_observed = phylo_dict[example_environment]['phylo'][d]['mean_diversity_observed']
    mean_diversity_observed_example_phylo.append(mean_diversity_observed)

    mean_diversity_predicted = phylo_dict[example_environment]['phylo'][d]['mean_diversity_predicted']
    mean_diversity_predicted_example_phylo.append(mean_diversity_predicted)



richness_taxon_all = []
richness_phylo_all = []
diversity_taxon_all = []
diversity_phylo_all = []
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    if environment == 'soil metagenome':
        continue

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)
    
    for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

        if rank in taxon_dict[environment]['taxon']:

            color_taxon = color_environment_taxon[rank_idx+1]

            mean_richness_predicted = taxon_dict[environment]['taxon'][rank]['mean_richness_predicted']
            mean_richness_observed = taxon_dict[environment]['taxon'][rank]['mean_richness_observed']

            mean_diversity_predicted = taxon_dict[environment]['taxon'][rank]['mean_diversity_predicted']
            mean_diversity_observed = taxon_dict[environment]['taxon'][rank]['mean_diversity_observed']

            color_taxon = color_environment_taxon[rank_idx+1]

            ax_richness_taxon.scatter(mean_richness_observed, mean_richness_predicted, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)
            ax_diversity_taxon.scatter(mean_diversity_observed, mean_diversity_predicted, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)
            
            richness_taxon_all.extend([mean_richness_predicted, mean_richness_observed])
            diversity_taxon_all.extend([mean_diversity_predicted, mean_diversity_observed])

    
    distances = list(phylo_dict[environment]['phylo'].keys())
    distances.sort()
    color_environment_phylo = plot_utils.get_custom_cmap_phylo(environment, len(distances))
    for distance_idx, distance in enumerate(distances):

        mean_richness_predicted = phylo_dict[environment]['phylo'][distance]['mean_richness_predicted']
        mean_richness_observed = phylo_dict[environment]['phylo'][distance]['mean_richness_observed']

        mean_diversity_predicted = phylo_dict[environment]['phylo'][distance]['mean_diversity_predicted']
        mean_diversity_observed = phylo_dict[environment]['phylo'][distance]['mean_diversity_observed']
        
        color_phylo = color_environment_phylo[distance_idx]
        ax_richness_phylo.scatter(mean_richness_observed, mean_richness_predicted, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)
        ax_diversity_phylo.scatter(mean_diversity_observed, mean_diversity_predicted, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)
        
        richness_phylo_all.extend([mean_richness_predicted, mean_richness_observed])
        diversity_phylo_all.extend([mean_diversity_predicted, mean_diversity_observed])

# plot example richness taxon
color_taxon_example = numpy.asarray([plot_utils.get_custom_cmap_taxon(example_environment)[i+1] for i in range(len(idx_taxa_ranks))])
ax_rank_vs_richness_taxon.plot(idx_taxa_ranks, mean_richness_predicted_example_taxon, ls='-', lw=2, alpha=0.8, c='k')
ax_rank_vs_richness_taxon.scatter(idx_taxa_ranks, mean_richness_predicted_example_taxon, color='k', s=45, linewidth=0.8, edgecolors='k', zorder=3, label='Predicted')
ax_rank_vs_richness_taxon.scatter(idx_taxa_ranks, mean_richness_observed_example_taxon, color=color_taxon_example, s=45, linewidth=0.8, edgecolors='k', zorder=3, label='Observed')
ax_rank_vs_richness_taxon.legend(loc='upper right')


# plot example richness phylo
color_phylo_example = plot_utils.get_custom_cmap_phylo(example_environment, n=len(distance_example))[1:,:]
ax_rank_vs_richness_phylo.plot(distance_example, mean_richness_predicted_example_phylo, ls='-', lw=2, alpha=0.8, c='k')
ax_rank_vs_richness_phylo.scatter(distance_example, mean_richness_predicted_example_phylo, color='k', s=45, linewidth=0.8, edgecolors='k', zorder=3)
ax_rank_vs_richness_phylo.scatter(distance_example, mean_richness_observed_example_phylo, color=color_phylo_example, s=45, linewidth=0.8, edgecolors='k', zorder=3)


# plot example diversity taxon
ax_rank_vs_diversity_taxon.plot(idx_taxa_ranks, mean_diversity_predicted_example_taxon, ls='-', lw=2, alpha=0.8, c='k')
ax_rank_vs_diversity_taxon.scatter(idx_taxa_ranks, mean_diversity_predicted_example_taxon, color='k', s=45, linewidth=0.8, edgecolors='k', zorder=3)
ax_rank_vs_diversity_taxon.scatter(idx_taxa_ranks, mean_diversity_observed_example_taxon, color=color_taxon_example, s=45, linewidth=0.8, edgecolors='k', zorder=3)


# plot example diversity phylo
ax_rank_vs_diversity_phylo.plot(distance_example, mean_diversity_predicted_example_phylo, ls='-', lw=2, alpha=0.8, c='k')
ax_rank_vs_diversity_phylo.scatter(distance_example, mean_diversity_predicted_example_phylo, color='k', s=45, linewidth=0.8, edgecolors='k', zorder=3)
ax_rank_vs_diversity_phylo.scatter(distance_example, mean_diversity_observed_example_phylo, color=color_phylo_example, s=45, linewidth=0.8, edgecolors='k', zorder=3)



# plot example phylo



# format subplots 
ax_rank_vs_richness_taxon.set_xticks(idx_taxa_ranks)
ax_rank_vs_richness_taxon.set_xticklabels(taxa_ranks_label, fontsize=10)
ax_rank_vs_richness_taxon.set_ylabel('Mean richness', fontsize=12)

ax_rank_vs_richness_phylo.set_ylabel('Mean richness', fontsize=12)
ax_rank_vs_richness_phylo.set_xscale('log', base=10)
ax_rank_vs_richness_phylo.set_yscale('log', base=10)
ax_rank_vs_richness_phylo.set_xlabel('Phylogenetic distance', fontsize=12)

ax_rank_vs_diversity_taxon.set_xticks(idx_taxa_ranks)
ax_rank_vs_diversity_taxon.set_xticklabels(taxa_ranks_label, fontsize=10)
ax_rank_vs_diversity_taxon.set_ylabel('Mean diversity', fontsize=12)

ax_rank_vs_diversity_phylo.set_ylabel('Mean richness', fontsize=12)
ax_rank_vs_diversity_phylo.set_xscale('log', base=10)
ax_rank_vs_diversity_phylo.set_xlabel('Phylogenetic distance', fontsize=12)


########



richness_taxon_min = min(richness_taxon_all)
richness_taxon_max = max(richness_taxon_all)
ax_richness_taxon.plot([richness_taxon_min*0.8,richness_taxon_max*1.2],[richness_taxon_min*0.8,richness_taxon_max*1.2], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_richness_taxon.set_xlim([richness_taxon_min*0.8,richness_taxon_max*1.2])
ax_richness_taxon.set_ylim([richness_taxon_min*0.8,richness_taxon_max*1.2])
ax_richness_taxon.set_xscale('log', base=10)
ax_richness_taxon.set_yscale('log', base=10)
ax_richness_taxon.set_xlabel('Predicted mean richness', fontsize=12)
ax_richness_taxon.set_ylabel('Observed mean richness', fontsize=12)
ax_richness_taxon.legend(loc='upper left')


richness_phylo_min = min(richness_phylo_all)
richness_phylo_max = max(richness_phylo_all)
ax_richness_phylo.plot([richness_phylo_min*0.8,richness_phylo_max*1.2],[richness_phylo_min*0.8,richness_phylo_max*1.2], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_richness_phylo.set_xlim([richness_phylo_min*0.8,richness_phylo_max*1.2])
ax_richness_phylo.set_ylim([richness_phylo_min*0.8,richness_phylo_max*1.2])
ax_richness_phylo.set_xscale('log', base=10)
ax_richness_phylo.set_yscale('log', base=10)
ax_richness_phylo.set_xlabel('Predicted mean richness', fontsize=12)
ax_richness_phylo.set_ylabel('Observed mean richness', fontsize=12)





diversity_taxon_min = min(diversity_taxon_all)
diversity_taxon_max = max(diversity_taxon_all)
ax_diversity_taxon.plot([diversity_taxon_min*0.8,diversity_taxon_max*1.2],[diversity_taxon_min*0.8,diversity_taxon_max*1.2], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_diversity_taxon.set_xlim([diversity_taxon_min*0.8,diversity_taxon_max*1.2])
ax_diversity_taxon.set_ylim([diversity_taxon_min*0.8,diversity_taxon_max*1.2])
ax_diversity_taxon.set_xlabel('Predicted mean diversity', fontsize=12)
ax_diversity_taxon.set_ylabel('Observed mean diversity', fontsize=12)


diversity_phylo_min = min(diversity_phylo_all)
diversity_phylo_max = max(diversity_phylo_all)
ax_diversity_phylo.plot([diversity_phylo_min*0.8,diversity_phylo_max*1.2],[diversity_phylo_min*0.8,diversity_phylo_max*1.2], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_diversity_phylo.set_xlim([diversity_phylo_min*0.8,diversity_phylo_max*1.2])
ax_diversity_phylo.set_ylim([diversity_phylo_min*0.8,diversity_phylo_max*1.2])
ax_diversity_phylo.set_xlabel('Predicted mean diversity', fontsize=12)
ax_diversity_phylo.set_ylabel('Observed mean diversity', fontsize=12)









fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%smean_summary.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
