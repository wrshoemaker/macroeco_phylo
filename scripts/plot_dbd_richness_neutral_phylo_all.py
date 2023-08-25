from __future__ import division
import config
import os
import sys
import subprocess
import random
import re
import itertools
import pickle
from collections import Counter

import numpy
import diversity_utils
import scipy.stats as stats

import ete3
import tree_utils
import random
import dbd_utils
import plot_utils

import matplotlib.pyplot as plt
import dbd_richness_neutral

#dbd_utils.make_dbd_slm_dict()

dbd_dict = dbd_richness_neutral.load_richness_neutral_dbd_dict()




environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label
taxa_ranks.insert(0, 'ASV')



fig = plt.figure(figsize = (12.5, 20)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


phylo_distances_all = []
for environment_idx, environment in enumerate(environments_to_keep):

    phylo_distances = list(dbd_dict[environment]['phylo'].keys())
    phylo_distances_all.extend(phylo_distances)



phylo_distances_count =  Counter(phylo_distances_all)
phylo_distances_set = list(phylo_distances_count.keys())
phylo_distances_set.sort()

phylo_distances_to_keep = []
for key in phylo_distances_set:
    if phylo_distances_count[key] == 9:
        phylo_distances_to_keep.append(key)
        

phylo_distances_to_keep = [0.02767612370754228, 0.04604239376758782, 0.07659678234751838, 0.1274274985703134, 0.3526699214174661]
   
for environment_idx, environment in enumerate(environments_to_keep):

    color_environment_phylo = plot_utils.get_custom_cmap_phylo(environment, n=len(phylo_distances_to_keep))
    
    for coarse_rank_idx, coarse_rank in enumerate(phylo_distances_to_keep):

        observed = dbd_dict[environment]['phylo'][coarse_rank]['slope_all']
        predicted = dbd_dict[environment]['phylo'][coarse_rank]['slope_neutral_all']

        observed = numpy.asarray(observed)
        predicted = numpy.asarray(predicted)

        idx_to_keep = (observed>0) & (predicted>0)
        observed = observed[idx_to_keep]
        predicted = predicted[idx_to_keep]

        color_phylo = color_environment_phylo[coarse_rank_idx+1]


        ax = plt.subplot2grid((9, 6), (environment_idx, coarse_rank_idx))

        ax.scatter(observed, predicted, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)

        merged_ = numpy.concatenate((observed, predicted),axis=0)
        min_ = min(merged_)
        max_ = max(merged_)
        ax.set_xlim([min_*0.2,max_*4]) 
        ax.set_ylim([min_*0.2,max_*4])
        ax.plot([min_*0.2,max_*4], [min_*0.2,max_*4], lw=2,ls='--',c='k',zorder=1, label='1:1')

        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)
        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)


        if coarse_rank_idx == 0:
            ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize = 12)


        if environment_idx == 0:
            #ax.set_title(diversity_utils.taxa_ranks_label_capital_dict[coarse_rank], fontsize=12)
            ax.set_title('Phylo. dist. = %s' % str(round(coarse_rank, 3)), fontsize=11)







fig.text(0.3, 0.06, "Observed richness slope", va='center', fontsize=28)
fig.text(0.02, 0.5, "Predicted richness slope, UNTB", va='center', rotation='vertical', fontsize=28)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sdbd_slope_neutral_richness_phylo_all.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()




