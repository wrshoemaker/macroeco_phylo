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

import tree_utils

tree = tree_utils.get_emp_tree()
tree_copy = tree.copy()
# DQ798457.1.1383
leaf_names = tree_copy.get_leaf_names()






#Tree node  (-0x7ffff807e52b84ee), Tree node 'GQ358410.1.1239' (-0x7ffff807e52b84e7)

rarefied = False

environments_to_keep = diversity_utils.environments_to_keep


fig = plt.figure(figsize = (20, 20)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):
    for environment_idx, environment in enumerate(environment_chunk):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)

        samples = pres_abs_dict['samples']
        taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
        sys.stderr.write("Getting site-by-species matrix...\n")

        s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s = s_by_s/numpy.sum(s_by_s, axis=0)
        mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
        log10_mean_rel_s_by_s = numpy.log10(mean_rel_s_by_s)

        idx_to_keep = numpy.random.choice(len(rel_s_by_s), size=100, replace=False)
        mean_rel_s_by_s = mean_rel_s_by_s[idx_to_keep]
        taxa = taxa[idx_to_keep]

        n_otu = len(mean_rel_s_by_s)

        sys.stderr.write("Calculating pairwise distances...\n")
        taxa_for_tree = [str(t) for t in taxa]
        tree_copy_environment = tree_utils.subset_tree(taxa_for_tree, tree_copy)

        distance_all = []
        log_ratio_mean_all = []
        ratio_mean_all = []
        for i in range(n_otu):
            for j in range(i):

                #distance_ij = tree_copy.get_distance(str(taxa[i]), str(taxa[j]))
                distance_ij = tree_copy_environment.get_distance(taxa_for_tree[i], taxa_for_tree[j])
                max_min_ij = [mean_rel_s_by_s[i], mean_rel_s_by_s[j]]
                ratio_max_min_ij = max(max_min_ij)/min(max_min_ij)
                distance_all.append(distance_ij)
                #log_ratio_mean_all.append(numpy.absolute(log10_mean_rel_s_by_s[i] - log10_mean_rel_s_by_s[j]))
                ratio_mean_all.append(ratio_max_min_ij)


        distance_all = numpy.asarray(distance_all)
        ratio_mean_all = numpy.asarray(ratio_mean_all)

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        x, y, z = plot_utils.get_scatter_density_arrays_for_loglog(distance_all, ratio_mean_all)
        ax.scatter(x, y, c=numpy.sqrt(z), cmap='Blues', s=70, alpha=0.9, edgecolors='none', zorder=1)
        #ax.scatter(x, y, c=numpy.sqrt(z), cmap='Blues', s=70, alpha=0.9, edgecolors='none', zorder=1)

        ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)
        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)



fig.text(0.3, 0.06, "Pairwise phylogenetic distance", va='center', fontsize=28)
fig.text(0.02, 0.5, "Max. mean abundance / min. mean abundance", va='center', rotation='vertical', fontsize=28)


fig.subplots_adjust(hspace=0.37, wspace=0.37)
fig_name = "%sphylo_distance_vs_mean_ratio.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
    

