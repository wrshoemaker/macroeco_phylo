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



#dbd_utils.make_richness_diversity_prediction_taxon_dict(iter_=100)

# load taxon dict
taxon_dict = dbd_utils.load_richness_diversity_prediction_taxon_dict()

# load phylo dict
phylo_dict = dbd_utils.load_richness_diversity_prediction_phylo_dict()



fig = plt.figure(figsize = (8, 16)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

# (#rows, # columns)
ax_richness_mean_taxon = plt.subplot2grid((4, 2), (0, 0))
ax_richness_var_taxon = plt.subplot2grid((4, 2), (1, 0))
ax_diversity_mean_taxon = plt.subplot2grid((4, 2), (2, 0))
ax_diversity_var_taxon = plt.subplot2grid((4, 2), (3, 0))

ax_richness_mean_phylo = plt.subplot2grid((4, 2), (0, 1))
ax_richness_var_phylo = plt.subplot2grid((4, 2), (1, 1))
ax_diversity_mean_phylo = plt.subplot2grid((4, 2), (2, 1))
ax_diversity_var_phylo = plt.subplot2grid((4, 2), (3, 1))


ax_richness_mean_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_richness_mean_taxon.transAxes)
ax_richness_mean_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_richness_mean_phylo.transAxes)

ax_richness_var_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[2], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_richness_var_taxon.transAxes)
ax_richness_var_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[3], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_richness_var_phylo.transAxes)

ax_diversity_mean_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[4], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_diversity_mean_taxon.transAxes)
ax_diversity_mean_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[5], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_diversity_mean_phylo.transAxes)

ax_diversity_var_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[6], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_diversity_var_taxon.transAxes)
ax_diversity_var_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[7], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_diversity_var_phylo.transAxes)


# title 
ax_richness_mean_taxon.set_title("Taxonomic coarse-graining", fontsize=11)
ax_richness_mean_phylo.set_title("Phylogenetic coarse-graining", fontsize=11)


for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    if environment == 'soil metagenome':
        continue

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

    mean_richness_taxon_all = []
    var_richness_taxon_all = []
    mean_diversity_taxon_all = []
    var_diversity_taxon_all = []
    # go through taxonomic
    for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

        if rank in taxon_dict[environment]['taxon']:

            mean_richness_predicted = taxon_dict[environment]['taxon'][rank]['mean_richness_predicted']
            mean_richness_predicted_gamma = taxon_dict[environment]['taxon'][rank]['mean_richness_predicted_gamma']

            var_richness_predicted = taxon_dict[environment]['taxon'][rank]['var_richness_predicted']
            var_richness_predicted_gamma = taxon_dict[environment]['taxon'][rank]['var_richness_predicted_gamma']

            mean_diversity_predicted = taxon_dict[environment]['taxon'][rank]['mean_diversity_predicted']
            mean_diversity_predicted_gamma = taxon_dict[environment]['taxon'][rank]['mean_diversity_predicted_gamma']

            mean_diversity_observed = taxon_dict[environment]['taxon'][rank]['mean_diversity_observed']

            print(mean_diversity_observed, mean_diversity_predicted, mean_diversity_predicted_gamma)

            var_diversity_predicted = taxon_dict[environment]['taxon'][rank]['var_diversity_predicted']
            var_diversity_predicted_gamma = taxon_dict[environment]['taxon'][rank]['var_diversity_predicted_gamma']

            color_taxon = color_environment_taxon[rank_idx+1]

            ax_richness_mean_taxon.scatter(mean_richness_predicted, mean_richness_predicted_gamma, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)
            ax_richness_var_taxon.scatter(var_richness_predicted, var_richness_predicted_gamma, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)
            ax_diversity_mean_taxon.scatter(mean_diversity_predicted, mean_diversity_predicted_gamma, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)
            
            if (numpy.isnan(var_diversity_predicted) == False) and (numpy.isnan(var_diversity_predicted_gamma) == False):
            
                ax_diversity_var_taxon.scatter(var_diversity_predicted, var_diversity_predicted_gamma, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)
                var_diversity_taxon_all.append(var_diversity_predicted)
                var_diversity_taxon_all.append(var_diversity_predicted_gamma)

                rel_error = numpy.absolute(var_diversity_predicted - var_diversity_predicted_gamma)/var_diversity_predicted_gamma

                print(environment, rank, var_diversity_predicted, var_diversity_predicted_gamma, rel_error)
            
            # for xlim and ylim
            mean_richness_taxon_all.append(mean_richness_predicted)
            mean_richness_taxon_all.append(mean_richness_predicted_gamma)

            var_richness_taxon_all.append(var_richness_predicted)
            var_richness_taxon_all.append(var_richness_predicted_gamma)

            mean_diversity_taxon_all.append(mean_diversity_predicted)
            mean_diversity_taxon_all.append(mean_diversity_predicted_gamma)

            


    distances = list(phylo_dict[environment]['phylo'].keys())
    distances.sort()
    color_environment_phylo = plot_utils.get_custom_cmap_phylo(environment, len(distances))

    mean_richness_phylo_all = []
    var_richness_phylo_all = []
    mean_diversity_phylo_all = []
    var_diversity_phylo_all = []

    for distance_idx, distance in enumerate(distances):

        mean_richness_predicted = phylo_dict[environment]['phylo'][distance]['mean_richness_predicted']
        mean_richness_predicted_gamma = phylo_dict[environment]['phylo'][distance]['mean_richness_predicted_gamma']

        var_richness_predicted = phylo_dict[environment]['phylo'][distance]['var_richness_predicted']
        var_richness_predicted_gamma = phylo_dict[environment]['phylo'][distance]['var_richness_predicted_gamma']

        mean_diversity_predicted = phylo_dict[environment]['phylo'][distance]['mean_diversity_predicted']
        mean_diversity_predicted_gamma = phylo_dict[environment]['phylo'][distance]['mean_diversity_predicted_gamma']

        var_diversity_predicted = phylo_dict[environment]['phylo'][distance]['var_diversity_predicted']
        var_diversity_predicted_gamma = phylo_dict[environment]['phylo'][distance]['var_diversity_predicted_gamma']

        
        color_phylo = color_environment_phylo[distance_idx]
        ax_richness_mean_phylo.scatter(mean_richness_predicted, mean_richness_predicted_gamma, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)
        ax_richness_var_phylo.scatter(var_richness_predicted, var_richness_predicted_gamma, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)
        ax_diversity_mean_phylo.scatter(mean_diversity_predicted, mean_diversity_predicted_gamma, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)
        
        if (numpy.isnan(var_diversity_predicted) == False) and (numpy.isnan(var_diversity_predicted_gamma) == False):
            ax_diversity_var_phylo.scatter(var_diversity_predicted, var_diversity_predicted_gamma, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)
            var_diversity_phylo_all.append(var_diversity_predicted)
            var_diversity_phylo_all.append(var_diversity_predicted_gamma)


        mean_richness_phylo_all.append(mean_richness_predicted)
        mean_richness_phylo_all.append(mean_richness_predicted_gamma)

        var_richness_phylo_all.append(var_richness_predicted)
        var_richness_phylo_all.append(var_richness_predicted_gamma)

        mean_diversity_phylo_all.append(mean_diversity_predicted)
        mean_diversity_phylo_all.append(mean_diversity_predicted_gamma)





richness_mean_taxon_min = min(mean_richness_taxon_all)
richness_mean_taxon_max = max(mean_richness_taxon_all)
ax_richness_mean_taxon.plot([richness_mean_taxon_min*0.5,richness_mean_taxon_max*1.1],[richness_mean_taxon_min*0.5,richness_mean_taxon_max*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_richness_mean_taxon.set_xlim([richness_mean_taxon_min*0.5,richness_mean_taxon_max*1.1])
ax_richness_mean_taxon.set_ylim([richness_mean_taxon_min*0.5,richness_mean_taxon_max*1.1])
ax_richness_mean_taxon.set_xlabel("Predicted mean richness, analytic", fontsize = 10)
ax_richness_mean_taxon.set_ylabel("Predicted mean richness, simulation", fontsize = 10)
ax_richness_mean_taxon.set_xscale('log', base=10)
ax_richness_mean_taxon.set_yscale('log', base=10)


richness_var_taxon_min = min(var_richness_taxon_all)
richness_var_taxon_max = max(var_richness_taxon_all)
ax_richness_var_taxon.plot([richness_var_taxon_min*0.5,richness_var_taxon_max*1.1],[richness_var_taxon_min*0.5,richness_var_taxon_max*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_richness_var_taxon.set_xlim([richness_var_taxon_min*0.5,richness_var_taxon_max*1.1])
ax_richness_var_taxon.set_ylim([richness_var_taxon_min*0.5,richness_var_taxon_max*1.1])
ax_richness_var_taxon.set_xlabel("Predicted variance of richness, analytic", fontsize = 10)
ax_richness_var_taxon.set_ylabel("Predicted variance of richness, simulation", fontsize = 10)
ax_richness_var_taxon.set_xscale('log', base=10)
ax_richness_var_taxon.set_yscale('log', base=10)


mean_diversity_taxon_min = min(mean_diversity_taxon_all)
mean_diversity_taxon_max = max(mean_diversity_taxon_all)
ax_diversity_mean_taxon.plot([mean_diversity_taxon_min*0.5,mean_diversity_taxon_max*1.1],[mean_diversity_taxon_min*0.5,mean_diversity_taxon_max*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_diversity_mean_taxon.set_xlim([mean_diversity_taxon_min*0.5,mean_diversity_taxon_max*1.1])
ax_diversity_mean_taxon.set_ylim([mean_diversity_taxon_min*0.5,mean_diversity_taxon_max*1.1])
ax_diversity_mean_taxon.set_xlabel("Predicted mean diversity, analytic", fontsize = 10)
ax_diversity_mean_taxon.set_ylabel("Predicted mean diversity, simulation", fontsize = 10)


var_diversity_taxon_min = min(var_diversity_taxon_all)
var_diversity_taxon_max = max(var_diversity_taxon_all)
ax_diversity_var_taxon.plot([var_diversity_taxon_min*0.3,var_diversity_taxon_max*1.4],[var_diversity_taxon_min*0.3,var_diversity_taxon_max*1.4], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_diversity_var_taxon.set_xlim([var_diversity_taxon_min*0.3,var_diversity_taxon_max*1.4])
ax_diversity_var_taxon.set_ylim([var_diversity_taxon_min*0.3,var_diversity_taxon_max*1.4])
ax_diversity_var_taxon.set_xlabel("Predicted variance of diversity, analytic", fontsize = 10)
ax_diversity_var_taxon.set_ylabel("Predicted variance of diversity, simulation", fontsize = 10)
#ax_diversity_var_taxon.set_xscale('log', base=10)
#ax_diversity_var_taxon.set_yscale('log', base=10)


# formatting for phylo
richness_mean_phylo_min = min(mean_richness_phylo_all)
richness_mean_phylo_max = max(mean_richness_phylo_all)
ax_richness_mean_phylo.plot([richness_mean_phylo_min*0.5,richness_mean_phylo_max*1.1],[richness_mean_phylo_min*0.5,richness_mean_phylo_max*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_richness_mean_phylo.set_xlim([richness_mean_phylo_min*0.5,richness_mean_phylo_max*1.1])
ax_richness_mean_phylo.set_ylim([richness_mean_phylo_min*0.5,richness_mean_phylo_max*1.1])
ax_richness_mean_phylo.set_xlabel("Predicted mean richness, analytic", fontsize = 10)
ax_richness_mean_phylo.set_ylabel("Predicted mean richness, simulation", fontsize = 10)
ax_richness_mean_phylo.set_xscale('log', base=10)
ax_richness_mean_phylo.set_yscale('log', base=10)


richness_var_phylo_min = min(var_richness_phylo_all)
richness_var_phylo_max = max(var_richness_phylo_all)
ax_richness_var_phylo.plot([richness_var_phylo_min*0.5,richness_var_phylo_max*1.1],[richness_var_phylo_min*0.5,richness_var_phylo_max*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_richness_var_phylo.set_xlim([richness_var_phylo_min*0.5,richness_var_phylo_max*1.1])
ax_richness_var_phylo.set_ylim([richness_var_phylo_min*0.5,richness_var_phylo_max*1.1])
ax_richness_var_phylo.set_xlabel("Predicted variance of richness, analytic", fontsize = 10)
ax_richness_var_phylo.set_ylabel("Predicted variance of richness, simulation", fontsize = 10)
ax_richness_var_phylo.set_xscale('log', base=10)
ax_richness_var_phylo.set_yscale('log', base=10)


mean_diversity_phylo_min = min(mean_diversity_phylo_all)
mean_diversity_phylo_max = max(mean_diversity_phylo_all)
ax_diversity_mean_phylo.plot([mean_diversity_phylo_min*0.5,mean_diversity_phylo_max*1.1],[mean_diversity_phylo_min*0.5,mean_diversity_phylo_max*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_diversity_mean_phylo.set_xlim([mean_diversity_phylo_min*0.5,mean_diversity_phylo_max*1.1])
ax_diversity_mean_phylo.set_ylim([mean_diversity_phylo_min*0.5,mean_diversity_phylo_max*1.1])
ax_diversity_mean_phylo.set_xlabel("Predicted mean diversity, analytic", fontsize = 10)
ax_diversity_mean_phylo.set_ylabel("Predicted mean diversity, simulation", fontsize = 10)


var_diversity_phylo_min = min(var_diversity_phylo_all)
var_diversity_phylo_max = max(var_diversity_phylo_all)
ax_diversity_var_phylo.plot([var_diversity_phylo_min*0.5,var_diversity_phylo_max*1.1],[var_diversity_phylo_min*0.5,var_diversity_phylo_max*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_diversity_var_phylo.set_xlim([var_diversity_phylo_min*0.5,var_diversity_phylo_max*1.1])
ax_diversity_var_phylo.set_ylim([var_diversity_phylo_min*0.5,var_diversity_phylo_max*1.1])
ax_diversity_var_phylo.set_xlabel("Predicted variance of diversity, analytic", fontsize = 10)
ax_diversity_var_phylo.set_ylabel("Predicted variance of diversity, simulation", fontsize = 10)






fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sfigure-3-figure-supplement-1-prediction_vs_simulation.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()




# load 