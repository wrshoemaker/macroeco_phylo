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

import diversity_utils
import config
import plot_utils

from scipy.special import rel_entr



if len(sys.argv) > 1:
    x_axis_n = True
else:
    x_axis_n = False


rarefied = False




fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)



environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]


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
color_all = plot_utils.make_blue_cmap(len(distances_all))



environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]

        distances = list(coarse_grained_tree_dict.keys())
        distances.sort()

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        pres_abs_dict = dbd_utils.load_sad_annotated_phylo_dict(environment, rarefied = rarefied)

        samples = pres_abs_dict['samples']
        taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))

        sys.stderr.write("Getting site-by-species matrix...\n")
        s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

        distances = list(coarse_grained_tree_dict.keys())
        distances.sort()

        rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

        afd_dict = {}
        for t_idx, t in enumerate(taxa):
            afd_dict[t] = {}
            afd_dict[t]['asv'] = {}
            afd_dict[t]['asv']['name'] = t
            afd_dict[t]['asv']['afd'] = rel_s_by_s[t_idx,:]
            afd_dict[t]['distances'] = {}


        sys.stderr.write("Running phylogenetic coarse-graining.....\n")
        for distance_idx, distance in enumerate(distances):

            sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
            coarse_grained_list = coarse_grained_tree_dict[distance]
            coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

            print(coarse_grained_list[0])

            # get indexes for each clade
            coarse_grained_all = [numpy.asarray([i for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            # coarse grain s-by-s for all clades
            s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
            rel_s_by_s_all_clades = (s_by_s_all_clades/s_by_s_all_clades.sum(axis=0))

            for coarse_grained_list_i_idx, coarse_grained_list_i in enumerate(coarse_grained_all):

                for taxon_i in coarse_grained_list_i:

                    if distance not in afd_dict[taxon_i]:
                        afd_dict[taxon_i]['distances'][distance] = {}

                    afd_dict[taxon_i]['distances'][distance]['afd'] = rel_s_by_s_all_clades[coarse_grained_list_i_idx,:]
                    afd_dict[taxon_i]['distances'][distance]['n'] = len(coarse_grained_idx_all[coarse_grained_list_i_idx])

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        asv_all = list(afd_dict.keys())
        for asv in asv_all:

            afd_asv = afd_dict[asv]['asv']['afd']

            if sum(afd_asv>0) <20:
                continue

            # run ranks
            distances_to_plot = []
            divergence_to_plot = []
            n_to_plot = []

            distances_asv = list(afd_dict[asv]['distances'].keys())

            for distance_idx, distance in enumerate(distances_asv):

                afd_rank = afd_dict[asv]['distances'][distance]['afd']

                idx_to_keep = (afd_rank>0) & (afd_asv>0)

                if sum(idx_to_keep) < 10:
                    continue

                afd_asv_to_keep = afd_asv[idx_to_keep]
                afd_rank_to_keep = afd_rank[idx_to_keep]

                afd_asv_to_keep = afd_asv_to_keep/sum(afd_asv_to_keep)
                afd_rank_to_keep = afd_rank_to_keep/sum(afd_rank_to_keep)

                rel_entr_sum = sum(rel_entr(afd_asv_to_keep, afd_rank_to_keep))

                n = afd_dict[asv]['distances'][distance]['n']

                if rel_entr_sum > 0:

                    distances_to_plot.append(distance)
                    divergence_to_plot.append(rel_entr_sum)
                    n_to_plot.append(n)

            if len(distances_to_plot) < 5:
                continue

            if x_axis_n == True:
                ax.plot(n_to_plot, divergence_to_plot, color='k', alpha=0.04, lw=1)

            else:
                ax.plot(distances_to_plot, divergence_to_plot, color='k', alpha=0.04, lw=1)


        ax.set_ylim([10**-4,10**1])

        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)
        ax.set_ylabel("KL divergence b/w ASV and coarse-grained AFD", fontsize = 7.5)

        if x_axis_n == True:
            ax.set_xlabel("Number coarse-grained ASVs", fontsize=9)

        else:
            ax.set_xlabel("Phylogenetic distance", fontsize = 9)

        ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)


if x_axis_n == True:
    x_axis_n_label = '_n'
else:
    x_axis_n_label = ''



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%skl_divergence_phylo%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied) + x_axis_n_label)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
