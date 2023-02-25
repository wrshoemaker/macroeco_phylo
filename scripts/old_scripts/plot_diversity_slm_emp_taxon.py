import os
#from Bio import Phylo
import config
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
import itertools

import dbd_utils
import diversity_utils
import mle_utils
import simulation_utils

import predict_diversity_slm_emp


taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = range(len(taxa_ranks))
taxa_ranks_label = diversity_utils.taxa_ranks_label




fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environment_chunk_all = [diversity_utils.environments_to_keep[x:x+3] for x in range(0, len(diversity_utils.environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))


        diversity_dict = predict_diversity_slm_emp.load_predict_diversity_taxon_dict(environment)

        mean_diveristy = [diversity_dict[rank]['observed'] for rank in diversity_utils.taxa_ranks]
        #mean_diveristy_null = [diversity_dict[environment][rank]['mean_diveristy_null'] for rank in diversity_utils.taxa_ranks]
        upper_ci = [diversity_dict[rank]['upper_ci'] for rank in diversity_utils.taxa_ranks]
        lower_ci = [diversity_dict[rank]['lower_ci'] for rank in diversity_utils.taxa_ranks]


        ax.plot(idx_taxa_ranks, mean_diveristy, ls='-', lw=1.5, c='k',  zorder=2)
        ax.scatter(idx_taxa_ranks, mean_diveristy, c='k', label='Observed',  zorder=3)
        ax.fill_between(idx_taxa_ranks, lower_ci, upper_ci, color='darkgrey', alpha=0.7, label='95% Prediction Interval', zorder=1)


        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)

        #ax.set_yscale('log', base=10)

        if environment_idx == 0:
            ax.set_ylabel("Mean Shannon's diversity", fontsize = 10)

        #if environment_chunk_idx == len(environment_chunk_all)-1:
        #    ax.set_xlabel("Phylogenetic distance", fontsize = 10)

        ax.set_xticks(idx_taxa_ranks)
        ax.set_xticklabels(taxa_ranks_label, fontsize=8)

        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="upper right", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sdiversity_slm_emp_taxon.png" % config.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
