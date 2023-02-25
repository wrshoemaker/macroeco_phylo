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


#to_plot = sys.argv[1]


#if to_plot == 'mean_diversity_null':
#    label_ = 'Permute ASV labels'
#    file_label = 'constrain_none'

#elif to_plot == 'mean_diversity_null_constrain_mad':
#    label_ = 'Permute ASV labels while constraining mean abundances'
#    file_label = 'constrain_mad'

#elif to_plot == 'mean_diversity_null_constrain_y':
#    label_ = 'Randomize '  +  r'$\bar{x}_{i}$'  + ', constrain ' + r'$x_{i,j}/\bar{x}_{i}$'
#    file_label = 'constrain_x_div_mad'

#else:
#    print('Argument not recognized!')





taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = range(len(taxa_ranks))
taxa_ranks_label = diversity_utils.taxa_ranks_label


#dbd_utils.make_diversity_vs_diversity_taxon_dict(rarefied=False)


diversity_vs_diversity_taxon_dict = dbd_utils.load_diversity_vs_diversity_taxon_dict(rarefied=False)


fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environment_chunk_all = [diversity_utils.environments_to_keep[x:x+3] for x in range(0, len(diversity_utils.environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        mean_diveristy = [diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity'] for rank in diversity_utils.taxa_ranks]
        #mean_diveristy_null = [diversity_dict[environment][rank]['mean_diveristy_null'] for rank in diversity_utils.taxa_ranks]


        upper_ci = [diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_null_upper_ci'] for rank in diversity_utils.taxa_ranks]
        lower_ci = [diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_null_lower_ci'] for rank in diversity_utils.taxa_ranks]

        upper_ci_constrain_mad = [diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_null_constrain_mad_upper_ci'] for rank in diversity_utils.taxa_ranks]
        lower_ci_constrain_mad = [diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_null_constrain_mad_lower_ci'] for rank in diversity_utils.taxa_ranks]


        ax.plot(idx_taxa_ranks, mean_diveristy, ls='-', lw=1.5, c='k',  zorder=2)
        ax.scatter(idx_taxa_ranks, mean_diveristy, c='k', label='Observed',  zorder=3)

        ax.fill_between(idx_taxa_ranks, lower_ci, upper_ci, color='dodgerblue', alpha=0.5, label='95% PI', zorder=1)
        ax.fill_between(idx_taxa_ranks, lower_ci_constrain_mad, upper_ci_constrain_mad, color='orangered', alpha=0.5, label='95% PI, constrain ' + r'$\bar{x}$', zorder=1)


        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)

        #ax.set_yscale('log', base=10)

        if environment_idx == 0:
            ax.set_ylabel("Mean Shannon's diversity", fontsize = 10)

        #if environment_chunk_idx == len(environment_chunk_all)-1:
        #    ax.set_xlabel("Phylogenetic distance", fontsize = 10)

        ax.set_xticks(idx_taxa_ranks)
        ax.set_xticklabels(taxa_ranks_label, fontsize=8)

        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="lower left", fontsize=7)



#fig.suptitle(label_, fontsize=14, y=0.95)


fig.subplots_adjust(hspace=0.25,wspace=0.3)
#fig_name = "%sdiversity_null_%s_emp_taxon.png" % (config.analysis_directory, file_label)
fig_name = "%sdiversity_null_emp_taxon.png" % (config.analysis_directory)

fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
