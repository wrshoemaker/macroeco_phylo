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

    #distances = list(coarse_grained_tree_dict.keys())
    #distances.sort()

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    #samples = diversity_utils.subset_observations(environment=environment)
    pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)

    samples = pres_abs_dict['samples']
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
    sys.stderr.write("Getting site-by-species matrix...\n")
    #s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)

    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

    for distance_idx, distance in enumerate(distances):

        sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
        coarse_grained_list = coarse_grained_tree_dict[distance]
        coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

        #if len(coarse_grained_n) == len(taxonomy_names):
        #    continue

        # get indexes for each clade
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)


        occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s_all_clades, s_by_s_all_clades.shape[0])

        idx_to_keep = (occupancies>0) & (predicted_occupancies > 0) & (mad > 0)
        mad = mad[idx_to_keep]
        occupancies = occupancies[idx_to_keep]
        predicted_occupancies = predicted_occupancies[idx_to_keep]

        ax = plt.subplot2grid((9, 6), (environment_idx, distance_idx))

        mad_occupancy_joint = numpy.concatenate((mad, occupancies),axis=0)
        min_ = min(mad_occupancy_joint)
        max_ = max(mad_occupancy_joint)

        sorted_plot_data = plot_utils.plot_color_by_pt_dens(mad, occupancies, radius=plot_utils.color_radius, loglog=1)
        x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]


        ax.scatter(x, y, c=numpy.sqrt(z), cmap='Blues', s=70, alpha=0.9, edgecolors='none', zorder=1)
        all_ = numpy.concatenate([x, y])


        # mad vs occupancy
        mad_log10 = numpy.log10(mad)
        occupancies_log10 = numpy.log10(occupancies)
        predicted_occupancies_log10 = numpy.log10(predicted_occupancies)
        hist_all, bin_edges_all = numpy.histogram(mad_log10, density=True, bins=25)
        bins_mean_all = [0.5 * (bin_edges_all[i] + bin_edges_all[i+1]) for i in range(0, len(bin_edges_all)-1 )]
        bins_mean_all_to_keep = []
        bins_occupancies = []
        for i in range(0, len(bin_edges_all)-1 ):
            predicted_occupancies_log10_i = predicted_occupancies_log10[(mad_log10>=bin_edges_all[i]) & (mad_log10<bin_edges_all[i+1])]
            #bins_mean_all_to_keep.append(bins_mean_all[i])
            bins_mean_all_to_keep.append(bin_edges_all[i])
            bins_occupancies.append(numpy.mean(predicted_occupancies_log10_i))


        bins_mean_all_to_keep = numpy.asarray(bins_mean_all_to_keep)
        bins_occupancies = numpy.asarray(bins_occupancies)

        bins_mean_all_to_keep_no_nan = bins_mean_all_to_keep[(~numpy.isnan(bins_mean_all_to_keep)) & (~numpy.isnan(bins_occupancies))]
        bins_occupancies_no_nan = bins_occupancies[(~numpy.isnan(bins_mean_all_to_keep)) & (~numpy.isnan(bins_occupancies))]

        ax.plot(10**bins_mean_all_to_keep_no_nan, 10**bins_occupancies_no_nan, lw=3, ls='--',c='k', zorder=2, label='Gamma prediction')

        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)

        if environment_idx == 0:
            ax.set_title('Phylo. dist. = %s' % str(round(distance, 3)), fontsize=11)

        if distance_idx == 0:
            ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize = 12)

        if (environment_idx == 0) and (distance_idx == 0):
            ax.legend(loc="upper left", fontsize=6)


        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)



fig.text(0.3, 0.06, "Mean relative abundance", va='center', fontsize=28)
fig.text(0.02, 0.5, "Occupancy", va='center', rotation='vertical', fontsize=28)


fig.subplots_adjust(hspace=0.37, wspace=0.37)
fig_name = "%sabundance_vs_occupancy_all_phylo%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
