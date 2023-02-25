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
import matplotlib as mpl
import dbd_utils

import diversity_utils
import config



rarefied = True


fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environments_to_keep = diversity_utils.environments_to_keep

coarse_grained_tree_dict_all = {}
distances_all = []
for environment in environments_to_keep:
    sys.stderr.write("Loading tree dict for %s...\n" % environment)
    coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_dict(environment=environment, rarefied=rarefied)
    coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict
    distances = list(coarse_grained_tree_dict.keys())
    distances_all.extend(distances)

distances_all = list(set(distances_all))
distances_all.sort()
color_all = diversity_utils.make_blue_cmap(len(distances_all))


environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]

        distances = list(coarse_grained_tree_dict.keys())
        distances.sort()

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #samples = diversity_utils.subset_observations(environment=environment)
        pres_abs_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
        samples = pres_abs_dict['samples']

        print(environment)

        diversity_vs_diversity_dict = dbd_utils.load_diversity_vs_diversity_phylo_emp_dict(environment, rarefied)

        print(diversity_vs_diversity_dict.keys())

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        for distance_idx, distance in enumerate(distances):

            observed = diversity_vs_diversity_dict[distance]['max_mean_fine_vs_mean_coarse']['observed']
            null = diversity_vs_diversity_dict[distance]['max_mean_fine_vs_mean_coarse']['null_mean']

            ax.scatter(observed, null, s=10, color=color_all(distances_all.index(distance)), alpha=0.9, lw=2)


        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)
        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)

        if environment_idx == 0:
            ax.set_ylabel("Probability density", fontsize = 10)


        if environment_chunk_idx == len(environment_chunk_all)-1:
            ax.set_xlabel("Rescaled log relative abundance", fontsize = 10)

        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="upper left", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%scoarse_dominance_taxon%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
