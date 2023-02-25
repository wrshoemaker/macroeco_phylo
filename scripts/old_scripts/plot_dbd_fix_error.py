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

rarefied = False


#dbd_utils.make_dbd_slm_dict()
dbd_dict = dbd_utils.load_richness_dbd_dict()

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]
environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]


taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label

taxa_ranks_idx_dict = {}
for t_idx, t in enumerate(taxa_ranks):
    taxa_ranks_idx_dict[t] = {}
    taxa_ranks_idx_dict[t]['idx'] = t_idx
    taxa_ranks_idx_dict[t]['label'] = taxa_ranks_label[t_idx]

# make the plot
fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        error_fix_beta = numpy.asarray([dbd_dict[environment]['taxon'][r]['mean_error_slope_slm_fix_beta'] for r in taxa_ranks if r in dbd_dict[environment]['taxon']])
        error_fix_mean = numpy.asarray([dbd_dict[environment]['taxon'][r]['mean_error_slope_slm_fix_mean'] for r in taxa_ranks if r in dbd_dict[environment]['taxon']])

        idx_taxa_ranks_environment = numpy.asarray([taxa_ranks_idx_dict[r]['idx'] for r in taxa_ranks if r in dbd_dict[environment]['taxon']])
        taxa_ranks_label_environment = numpy.asarray([taxa_ranks_idx_dict[r]['label'] for r in taxa_ranks if r in dbd_dict[environment]['taxon']])

        ax.plot(idx_taxa_ranks_environment, error_fix_beta, ls='-', lw=1.5, c='k',  zorder=2)
        #ax.scatter(idx_taxa_ranks, error_fix_beta, c=color_all, s=60, label='Observed', linewidth=0.8, edgecolors='k', zorder=3)
        ax.scatter(idx_taxa_ranks_environment, error_fix_beta, c='b', s=60, label='Fix beta', linewidth=0.8, edgecolors='k', zorder=3)
        
        ax.plot(idx_taxa_ranks_environment, error_fix_mean, ls='-', lw=1.5, c='k',  zorder=2)
        ax.scatter(idx_taxa_ranks_environment, error_fix_mean, c='r', s=60, label='Fix mean', linewidth=0.8, edgecolors='k', zorder=3)

        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)
        ax.set_xticks(idx_taxa_ranks_environment)
        ax.set_xticklabels(taxa_ranks_label_environment, fontsize=8)

        #ax.set_yscale('log', base=10)


        if environment_idx == 0:
            ax.set_ylabel("Mean relative errror", fontsize = 10)


        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="upper left", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%serror_slm_dbd_slope_taxon.png" % config.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()