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
        #samples = diversity_utils.subset_observations(environment=environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s, taxa)

        idx_to_keep = (occupancies>0) & (predicted_occupancies > 0)
        occupancies = occupancies[idx_to_keep]
        predicted_occupancies = predicted_occupancies[idx_to_keep]

        predicted_and_observed_occupancies = numpy.concatenate((occupancies,predicted_occupancies),axis=0)
        min_ = min(predicted_and_observed_occupancies)
        max_ = max(predicted_and_observed_occupancies)

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        # occupancy
        ax.plot([min_*0.5,1.08],[min_*0.5,1.08], lw=2,ls='--',c='k',zorder=2, label='1:1')

        ax.scatter(occupancies, predicted_occupancies, alpha=0.08, c='k', s=12, zorder=1)#, linewidth=0.8, edgecolors='k')


        ax.set_xlim([min_*0.5, 1.08])
        ax.set_ylim([min_*0.5, 1.08])

        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)
        ax.set_xlabel('Observed occupancy', fontsize=11)
        ax.set_ylabel('Predicted occupancy', fontsize=11)
        ax.tick_params(axis='both', which='minor', labelsize=9)
        ax.tick_params(axis='both', which='major', labelsize=9)

        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)


        #ax_occupancy.legend(loc="upper left", fontsize=8)





fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%soccupancy_no_subsampling.png" % config.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
