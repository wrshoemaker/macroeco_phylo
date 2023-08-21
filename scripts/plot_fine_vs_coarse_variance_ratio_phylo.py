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
import plot_utils
from collections import Counter



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

#distances_all = list(set(distances_all))
#distances_all.sort()

distances_dict = Counter(distances_all)
distances = []
for key, value in distances_dict.items():
    if value == 9:
        distances.append(key)

distances.sort()
distances = numpy.asarray(distances)
distances = distances[1:]
idx = numpy.round(numpy.linspace(0, len(distances) - 1, 50)).astype(int)
distances = distances[idx]


distances = numpy.asarray([0.01663614249384222, 0.02767612370754228, 0.05938601867590266, 0.16435747993726982, 0.3526699214174661])

color_all = plot_utils.make_blue_cmap(len(distances_all))



fig = plt.figure(figsize = (10, 20)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


for environment_idx, environment in enumerate(environments_to_keep):

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]
    pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)

    samples = pres_abs_dict['samples']
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))

    sys.stderr.write("Getting site-by-species matrix...\n")
    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

    #distances = list(coarse_grained_tree_dict.keys())
    #distances.sort()

    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)

    sys.stderr.write("Running phylogenetic coarse-graining.....\n")
    afd_dict = {}

    for distance_idx, distance in enumerate(distances):

        sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))

        coarse_grained_list = coarse_grained_tree_dict[distance]
        coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

        # get indexes for each clade
        coarse_grained_all = [numpy.asarray([i for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list if len(coarse_grained_list_i) >= 5]
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list if len(coarse_grained_list_i) >= 5]

        coarse_grained_var_all = [numpy.var(rel_s_by_s[coarse_grained_idx_i,:].sum(axis=0)) for coarse_grained_idx_i in coarse_grained_idx_all]
        sum_var_all = [sum(var_rel_s_by_s[coarse_grained_idx_i]) for coarse_grained_idx_i in coarse_grained_idx_all]

        if (len(coarse_grained_var_all) == 0) or (len(sum_var_all) == 0):
            continue

        ax = plt.subplot2grid((9, 5), (environment_idx, distance_idx))

        sum_var_all = numpy.asarray(sum_var_all)
        coarse_grained_var_all = numpy.asarray(coarse_grained_var_all)

        #plot_utils.get_scatter_density_arrays_for_loglog()
        

        ratio_ = coarse_grained_var_all/sum_var_all

        x, y, z = plot_utils.get_scatter_density_arrays_for_loglog(sum_var_all, ratio_)
        ax.scatter(x, y, c=numpy.sqrt(z), cmap='Blues', s=70, alpha=0.9, edgecolors='none', zorder=1)

        #ax.scatter(sum_var_all, coarse_grained_var_all, s=11, color='k', alpha=0.35, lw=2, zorder=2)

        min_ = min(sum_var_all)
        max_ = max(sum_var_all)
        #ax.plot([min_,max_],[min(ratio_),max(ratio_)], lw=3,ls=':',c='k',zorder=1)
        max_y = max([min(ratio_), max(ratio_)])
        ax.set_xlim([min_, max_])
        ax.set_ylim([1/max_y, max_y])

        ax.axhline(y=1, lw=3, ls=':',c='k', zorder=1)

        if environment_idx == 0:
            ax.set_title('Phylo. dist. = %s' % str(round(distance, 3)), fontsize=11)

        if distance_idx == 0:
            ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize = 12)


        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        
        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)


        #ax.set_xlabel("Sum of ASV variances", fontsize = 9)
        #ax.set_ylabel("Coarse-grained variance", fontsize = 9)


fig.text(0.3, 0.06, "Sum of OTU variances", va='center', fontsize=28)
fig.text(0.02, 0.5, "Ratio of coarse-grained variance and sum of variances", va='center', rotation='vertical', fontsize=28)


fig.subplots_adjust(hspace=0.37,wspace=0.37)
fig_name = "%sfine_vs_coarse_variance_ratio_phylo%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
