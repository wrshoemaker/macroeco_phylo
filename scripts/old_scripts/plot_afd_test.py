import os
#from Bio import Phylo
import random
import copy
import config
import sys
import numpy
import random
import pickle
import scipy.stats as stats
import matplotlib.pyplot as plt
import functools
import operator
import matplotlib as mpl
import dbd_utils

import diversity_utils
import plot_utils

import mle_utils


rarefied = False


fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environments_to_keep = diversity_utils.environments_to_keep

environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #samples = diversity_utils.subset_observations(environment=environment)
        pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)
        samples = pres_abs_dict['samples']

        taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
        s_by_s_taxon = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s_taxon = (s_by_s_taxon/s_by_s_taxon.sum(axis=0))


        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        for afd in rel_s_by_s_taxon:

            afd = afd[afd>0]
            afd_log10 = numpy.log10(afd)

            if len(afd_log10) < 80:
                continue

            afd_log10_rescaled = (afd_log10 - numpy.mean(afd_log10))/numpy.std(afd_log10)

            hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(afd_log10_rescaled, bins=15)

            #ax.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=plot_utils.rgb_blue_taxon(0), alpha=0.9, lw=2, label='ASV')
            ax.plot(bins_mean_to_plot, hist_to_plot, color=plot_utils.rgb_blue_taxon(4), alpha=0.2, lw=1, label='ASV')


        ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)
        ax.set_yscale('log', base=10)

        if environment_idx == 0:
            ax.set_ylabel("Probability density", fontsize = 10)


        if environment_chunk_idx == len(environment_chunk_all)-1:
            ax.set_xlabel("Rescaled log relative abundance", fontsize = 10)

        #if (environment_chunk_idx ==0) and (environment_idx == 0):
         #   ax.legend(loc="upper left", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%stest_afd_taxon%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
