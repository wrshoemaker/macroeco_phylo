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



fig = plt.figure(figsize = (12, 4)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

# (#rows, # columns)
ax_observed_vs_predicted = plt.subplot2grid((1, 3), (0, 0))
ax_observed_vs_predicted_gamma = plt.subplot2grid((1, 3), (0, 1))
ax_predicted_vs_predicted_gamma = plt.subplot2grid((1, 3), (0, 2))


predicted_gamma_all = []
predicted_all = []
observed_all = []
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    if environment == 'soil metagenome':
        continue

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)



    distances = list(phylo_dict[environment]['phylo'].keys())
    distances.sort()
    color_environment_phylo = plot_utils.get_custom_cmap_phylo(environment, len(distances))

    

    for distance_idx, distance in enumerate(distances):


        mean_diversity_predicted = phylo_dict[environment]['phylo'][distance]['mean_diversity_predicted']
        mean_diversity_predicted_gamma = phylo_dict[environment]['phylo'][distance]['mean_diversity_predicted_gamma']
        mean_diversity_observed = phylo_dict[environment]['phylo'][distance]['mean_diversity_observed']

        
        color_phylo = color_environment_phylo[distance_idx]

        ax_observed_vs_predicted.scatter(mean_diversity_observed, mean_diversity_predicted, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)
        ax_observed_vs_predicted_gamma.scatter(mean_diversity_observed, mean_diversity_predicted_gamma, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)
        ax_predicted_vs_predicted_gamma.scatter(mean_diversity_predicted, mean_diversity_predicted_gamma, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)
        
        predicted_gamma_all.append(mean_diversity_predicted_gamma)
        predicted_all.append(mean_diversity_predicted)
        observed_all.append(mean_diversity_observed)







predicted_gamma_min = min(predicted_gamma_all)
predicted_gamma_max = max(predicted_gamma_all)

predicted_min = min(predicted_all)
predicted_max = max(predicted_all)

observed_min = min(observed_all)
observed_max = max(observed_all)



ax_observed_vs_predicted.set_xlabel("Observed mean diversity", fontsize = 10)
ax_observed_vs_predicted.set_ylabel("Predicted mean diversity, simulation", fontsize = 10)
ax_observed_vs_predicted.plot([observed_min*0.5,observed_max*1.1],[predicted_min*0.5,predicted_max*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_observed_vs_predicted.set_xlim([observed_min*0.5,observed_max*1.1])
ax_observed_vs_predicted.set_ylim([predicted_min*0.5,predicted_max*1.1])



ax_observed_vs_predicted_gamma.set_xlabel("Observed mean diversity", fontsize = 10)
ax_observed_vs_predicted_gamma.set_ylabel("Predicted mean diversity, simulation", fontsize = 10)
ax_observed_vs_predicted_gamma.plot([observed_min*0.5,observed_max*1.1],[predicted_gamma_min*0.5,predicted_gamma_max*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_observed_vs_predicted_gamma.set_xlim([observed_min*0.5,observed_max*1.1])
ax_observed_vs_predicted_gamma.set_ylim([predicted_gamma_min*0.5,predicted_gamma_max*1.1])

ax_predicted_vs_predicted_gamma.set_xlabel("Predicted mean diversity, analytic", fontsize = 10)
ax_predicted_vs_predicted_gamma.set_ylabel("Predicted mean diversity, simulation", fontsize = 10)
ax_predicted_vs_predicted_gamma.plot([predicted_min*0.5,predicted_max*1.1],[predicted_gamma_min*0.5,predicted_gamma_max*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_predicted_vs_predicted_gamma.set_xlim([predicted_min*0.5,predicted_max*1.1])
ax_predicted_vs_predicted_gamma.set_ylim([predicted_gamma_min*0.5,predicted_gamma_max*1.1])








fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sprediction_vs_simulation_vs_observed.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()




# load 