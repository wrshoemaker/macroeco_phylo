
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

from scipy.special import rel_entr



rarefied = False
n_iter = 1000

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

coarse_grained_tree_dict_all = {}
distances_all = []
for environment in environments_to_keep:
    sys.stderr.write("Loading tree dict for %s...\n" % environment)
    coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)
    coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict
    distances = list(coarse_grained_tree_dict.keys())
    distances_all.extend(distances)

distances_all = list(set(distances_all))
distances_all.sort()
color_all = plot_utils.make_blue_cmap(len(distances_all))


distances = numpy.asarray([0.01663614249384222, 0.02767612370754228, 0.05938601867590266, 0.16435747993726982, 0.3526699214174661])



fig = plt.figure(figsize = (12, 20)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


for environment_idx, environment in enumerate(environments_to_keep):

    coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    #samples = diversity_utils.subset_observations(environment=environment)
    pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)

    samples = pres_abs_dict['samples']
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
    sys.stderr.write("Getting site-by-species matrix...\n")
    #s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)

    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

    rel_s_by_s = s_by_s/numpy.sum(s_by_s, axis=0)

    rel_s_by_s_copy = numpy.copy(rel_s_by_s)

    for distance_idx, distance in enumerate(distances):

        sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
        coarse_grained_list = coarse_grained_tree_dict[distance]
        coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

        # get indexes for each clade
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)

        rel_s_by_s_all_clades = s_by_s_all_clades/s_by_s_all_clades.sum(axis=0)
        mad = numpy.mean(rel_s_by_s_all_clades, axis=1)

        coarse_grain_n_idx = numpy.append([0], numpy.cumsum(coarse_grained_n))[:-1]
        # make null
        mad_null_all = []
        for i in range(n_iter):

            numpy.random.shuffle(rel_s_by_s_copy)
            rel_s_by_s_all_clades_null = numpy.add.reduceat(rel_s_by_s_copy, coarse_grain_n_idx, axis=0)
            mad_null_i = numpy.mean(rel_s_by_s_all_clades_null, axis=1)
            mad_null_all.extend(mad_null_i.tolist())
            

        ax = plt.subplot2grid((9, 6), (environment_idx, distance_idx))


        rescaled_mad = numpy.log10(mad) - numpy.log10(1/len(mad))
        rescaled_mad_null_all = numpy.log10(mad_null_all) - numpy.log10(1/len(mad))

        
        hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(rescaled_mad)
        ax.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=color_all(distances_all.index(distance)), alpha=0.9, lw=2)

        hist_to_plot_null, bins_mean_to_plot_null = diversity_utils.get_hist_and_bins(rescaled_mad_null_all)
        ax.scatter(bins_mean_to_plot_null, hist_to_plot_null, s=10, color='darkgrey', alpha=0.9, lw=2)


        ax.set_ylim([min(hist_to_plot), max(hist_to_plot)])
        ax.set_yscale('log', base=10)
        #ax.set_xlabel('Rescaled MAD', fontsize=11)
        #ax.set_ylabel('Probability dens', fontsize=11)
        ax.tick_params(axis='both', which='minor', labelsize=8)
        ax.tick_params(axis='both', which='major', labelsize=8)
        

        if environment_idx == 0:
            ax.set_title('Phylo. dist. = %s' % str(round(distance, 3)), fontsize=11)

        if distance_idx == 0:
            ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize = 12)

        #if (environment_idx == 0) and (distance_idx == 0):
        #    ax.legend(loc="upper left", fontsize=7)




fig.text(0.25, 0.06, "Rescaled MAD", va='center', fontsize=28)
fig.text(0.02, 0.5, "Probability density", va='center', rotation='vertical', fontsize=28)


fig.subplots_adjust(hspace=0.37, wspace=0.37)
fig_name = "%scoarse_mad_with_null_phylo.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
