
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

        # get indexes for each clade
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)

        rel_s_by_s_all_clades = s_by_s_all_clades/s_by_s_all_clades.sum(axis=0)
        mean_all = numpy.mean(rel_s_by_s_all_clades, axis=1)
        var_all = numpy.var(rel_s_by_s_all_clades, axis=1)

        ax = plt.subplot2grid((9, 6), (environment_idx, distance_idx))

        sorted_plot_data = plot_utils.plot_color_by_pt_dens(mean_all, var_all, radius=plot_utils.color_radius, loglog=1)
        x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]
        ax.scatter(x, y, c=numpy.sqrt(z), cmap='Blues', s=70, alpha=0.9, edgecolors='none', zorder=1)

        slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(mean_all), numpy.log10(var_all))

        x_log10_range =  numpy.linspace(min(numpy.log10(mean_all)) , max(numpy.log10(mean_all)) , 10000)
        y_log10_fit_range = 10 ** (slope*x_log10_range + intercept)
        #y_log10_null_range = 10 ** (slope_null*x_log10_range + intercept)

        ax.plot(10**x_log10_range, y_log10_fit_range, c='k', lw=2.5, linestyle='--', zorder=2, label="OLS regression")

        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)
        #ax.set_xlabel('Mean relative abundance', fontsize=11)
        #ax.set_ylabel('Occupancy', fontsize=11)
        ax.tick_params(axis='both', which='minor', labelsize=9)
        ax.tick_params(axis='both', which='major', labelsize=9)


        if environment_idx == 0:
            ax.set_title('Phylo. dist. = %s' % str(round(distance, 3)), fontsize=11)

        if distance_idx == 0:
            ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize = 12)

        if (environment_idx == 0) and (distance_idx == 0):
            ax.legend(loc="upper left", fontsize=7)




fig.text(0.25, 0.06, "Mean relative abundance", va='center', fontsize=28)
fig.text(0.02, 0.5, "Variance of relative abundance", va='center', rotation='vertical', fontsize=28)


fig.subplots_adjust(hspace=0.37, wspace=0.37)
fig_name = "%staylors_law_phylo.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
