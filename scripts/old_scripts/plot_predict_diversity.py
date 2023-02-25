import os
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
import predict_diversity_slm_emp

import diversity_utils
import config
import plot_utils


rarefied = False
iter = 1000

fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label


color_all = numpy.asarray([plot_utils.rgb_blue_taxon(r) for r in idx_taxa_ranks])
environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        diversity_dict = predict_diversity_slm_emp.load_diversity_gamma(environment)
        observed = diversity_dict['observed']
        predicted = diversity_dict['mean_null']

        observed = numpy.asarray(observed)
        predicted = numpy.asarray(predicted)

        predicted_and_observed = numpy.concatenate((observed,predicted),axis=0)
        min_ = min(predicted_and_observed)
        max_ = max(predicted_and_observed)

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        # occupancy
        ax.plot([min_,max_],[min_,max_], lw=2,ls='--',c='k',zorder=2, label='1:1')

        ax.scatter(observed, predicted, alpha=1, c='k', s=12, zorder=1)#, linewidth=0.8, edgecolors='k')


        ax.set_xlim([min_,max_])
        ax.set_ylim([min_,max_])


        ax.set_xlabel('Observed diversity', fontsize=11)
        ax.set_ylabel('Predicted diversity', fontsize=11)
        ax.tick_params(axis='both', which='minor', labelsize=9)
        ax.tick_params(axis='both', which='major', labelsize=9)

        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)


fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%spredict_diversity.png" % config.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
