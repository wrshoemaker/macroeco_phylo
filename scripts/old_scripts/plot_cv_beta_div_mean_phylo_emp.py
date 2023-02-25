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


rarefied = False

fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

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



environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]
        pres_abs_dict = dbd_utils.load_sad_annotated_phylo_dict(environment, rarefied = rarefied)

        samples = pres_abs_dict['samples']
        taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))

        sys.stderr.write("Getting site-by-species matrix...\n")
        s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

        distances = list(coarse_grained_tree_dict.keys())
        distances.sort()

        rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
        beta_div_mean = numpy.mean(rel_s_by_s, axis=1)/numpy.var(rel_s_by_s, axis=1)

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))
        beta_div_mean_dict = {}
        sys.stderr.write("Running phylogenetic coarse-graining.....\n")
        for distance_idx, distance in enumerate(distances):

            sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
            coarse_grained_list = coarse_grained_tree_dict[distance]
            coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

            # get indexes for each clade
            coarse_grained_all = [numpy.asarray([i for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            #coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            coarse_grained_beta_div_mean_all = [numpy.asarray([beta_div_mean[numpy.where(taxa==i)[0][0]] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            # coarse grain s-by-s for all clades
            #s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
            #rel_s_by_s_all_clades = (s_by_s_all_clades/s_by_s_all_clades.sum(axis=0))

            for coarse_grained_list_i_idx, coarse_grained_list_i in enumerate(coarse_grained_all):

                #coarse_grained_idx_all_i = coarse_grained_idx_all[coarse_grained_list_i_idx]
                coarse_grained_beta_div_mean_i = coarse_grained_beta_div_mean_all[coarse_grained_list_i_idx]

                if len(coarse_grained_list_i) < 5:
                    continue

                cv_coarse_grained_beta_div_mean_i = numpy.std(coarse_grained_beta_div_mean_i)/numpy.mean(coarse_grained_beta_div_mean_i)
                for taxon_i in coarse_grained_list_i:

                    if taxon_i not in beta_div_mean_dict:
                        beta_div_mean_dict[taxon_i] = {}

                    beta_div_mean_dict[taxon_i][distance] = cv_coarse_grained_beta_div_mean_i


        species = list(beta_div_mean_dict.keys())
        for s in species:

            distances_s = list(beta_div_mean_dict[s].keys())
            distances_s.sort()
            distances_s = numpy.asarray(distances_s)

            beta_div_mean = numpy.asarray([beta_div_mean_dict[s][d] for d in distances_s])
            idx_to_keep = ~numpy.isnan(beta_div_mean)

            if sum(idx_to_keep) < 5:
                continue

            distances_s_to_plot = distances_s[idx_to_keep]
            beta_div_mean_to_plot = beta_div_mean[idx_to_keep]
            ax.plot(distances_s_to_plot, beta_div_mean_to_plot, ls='-', lw=0.75, c='k', alpha=0.008,  zorder=2)


        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)

        #ax.set_xticks(idx_taxa_ranks)
        #ax.set_xticklabels(taxa_ranks_label, fontsize=8)

        ax.axhline(y=1, color='k', linestyle=':', lw = 2, zorder=1)
        ax.set_ylim([0,1.55])
        ax.set_xscale('log', base=10)

        if environment_idx == 0:
            ax.set_ylabel("CV of gamma rate parameter", fontsize = 10)

        ax.set_xlabel("Phylogenetic distance", fontsize = 10)

        #if environment_chunk_idx == len(environment_chunk_all)-1:
        #    ax.set_xlabel("Rescaled log relative abundance", fontsize = 10)

        #if (environment_chunk_idx ==0) and (environment_idx == 0):
        #    ax.legend(loc="upper left", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%scv_beta_div_mean_phylo_emp%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
