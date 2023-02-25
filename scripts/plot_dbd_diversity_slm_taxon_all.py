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


#dbd_utils.make_dbd_slm_dict()

dbd_dict = dbd_utils.load_diversity_slm_dbd_dict()


environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label
taxa_ranks.insert(0, 'ASV')



fig = plt.figure(figsize = (12, 20)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


for environment_idx, environment in enumerate(environments_to_keep):

    taxa_ranks = list(dbd_dict[environment]['taxon'].keys())

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

    for coarse_rank_idx, coarse_rank in enumerate(diversity_utils.taxa_ranks[1:]):

        if coarse_rank in taxa_ranks:

            if 'slope_all' not in dbd_dict[environment]['taxon'][coarse_rank]:
                continue
        
            observed = dbd_dict[environment]['taxon'][coarse_rank]['slope_all']
            predicted = dbd_dict[environment]['taxon'][coarse_rank]['slope_slm_all']

            observed = numpy.asarray(observed)
            predicted = numpy.asarray(predicted)

            #idx_to_keep = (observed!=0) & (predicted!=0)

            #if sum(idx_to_keep) == 0:
            #    continue

            #observed = observed[idx_to_keep]
            #predicted = predicted[idx_to_keep]

            color_taxon = color_environment_taxon[coarse_rank_idx+1]

            ax = plt.subplot2grid((9, 6), (environment_idx, coarse_rank_idx))

            ax.scatter(observed, predicted, color=color_taxon, s=40, alpha=0.8, linewidth=0.8, edgecolors='k', zorder=2)

            merged_ = numpy.concatenate((observed, predicted),axis=0)
            min_ = min(merged_)
            max_ = max(merged_)
            print(min_)
            ax.set_xlim([min_,max_*1.2]) 
            ax.set_ylim([min_,max_*1.2])
            ax.plot([min_,max_*1.2], [min_,max_*1.2], lw=2,ls='--',c='k',zorder=1, label='1:1')

            #ax.set_xscale('log', base=10)
            #ax.set_yscale('log', base=10)
            ax.tick_params(axis='x', labelsize=8)
            ax.tick_params(axis='y', labelsize=8)

            #if coarse_rank_idx == 0:
            #    #ax.set_ylabel("Predicted DBD slope", fontsize = 10)
            #    ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize = 12)

            if coarse_rank_idx == 0:
                ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize = 12)

#            taxa_ranks_label_capital_dict

            #if environment_idx == len(environments_to_keep)-1:
                #ax.set_xlabel("Observed DBD slope", fontsize = 10)
                #ax.set_xlabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize = 12)

            if environment_idx == 0:
                ax.set_title(diversity_utils.taxa_ranks_label_capital_dict[coarse_rank], fontsize=12)

            #if environment_idx == 0:
            #    ax.set_title(diversity_utils.taxa_ranks_label_with_asv[rank_idx], fontsize=12)


                

            #if rank_idx == 0:
            #    ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize = 12)



fig.text(0.2, 0.06, "Observed diversity slope", va='center', fontsize=28)
fig.text(0.02, 0.5, "Predicted diversity slope", va='center', rotation='vertical', fontsize=28)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sdbd_slope_slm_diversity_taxon_all.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()




