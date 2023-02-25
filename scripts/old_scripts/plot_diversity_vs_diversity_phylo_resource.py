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


dbd_utils.make_diversity_vs_diversity_phylo_dalbello_dict()

diversity_vs_diversity_phylo_resource_dict = dbd_utils.load_diversity_vs_diversity_phylo_resource_dict()
n_resource_all = list(diversity_vs_diversity_phylo_resource_dict.keys())
n_resource_all.sort()





def plot_slopes():

    fig = plt.figure(figsize = (10, 10)) #
    fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

    n_resource_chunks = [n_resource_all[x:x+3] for x in range(0, len(n_resource_all), 3)]
    for n_resource_chunk_idx, n_resource_chunk in enumerate(n_resource_chunks):

        for n_resource_idx, n_resource in enumerate(n_resource_chunk):

            distances = list(diversity_vs_diversity_phylo_resource_dict[n_resource].keys())
            distances.sort()

            lower_ci_all = [diversity_vs_diversity_phylo_resource_dict[n_resource][d]['lower_ci'] for d in distances]
            upper_ci_all = [diversity_vs_diversity_phylo_resource_dict[n_resource][d]['upper_ci'] for d in distances]
            slope_all = [diversity_vs_diversity_phylo_resource_dict[n_resource][d]['slope'] for d in distances]

            lower_ci_all = numpy.asarray(lower_ci_all)
            upper_ci_all = numpy.asarray(upper_ci_all)
            slope_all = numpy.asarray(slope_all)

            ax = plt.subplot2grid((3, 3), (n_resource_chunk_idx, n_resource_idx))

            ax.plot(distances, slope_all, ls='-', lw=1.5, c='k',  zorder=2)
            ax.scatter(distances, slope_all, c='k', s=8, label='Observed',  zorder=3)
            ax.fill_between(distances, lower_ci_all, upper_ci_all, color='darkgrey', alpha=0.7, label='95% CI', zorder=1)

            #ax.set_xticks(idx_taxa_tanks)
            #ax.set_xticklabels(taxa_ranks_format, fontsize=9)
            ax.set_ylabel('Fine vs. coarse-grained diversity slope', fontsize=9)
            ax.set_xlabel('Phylogenetic distance', fontsize=9)
            ax.set_title('%d resources' % n_resource, fontsize=11)
            ax.set_xscale('log', base=10)
            #ax.set_xlim([0,4])

            ax.tick_params(axis='both', which='major', labelsize=7)


            if (n_resource_chunk_idx==0) and (n_resource_idx==0):
                ax.legend(loc="lower left", fontsize=7)


    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%sdiversity_vs_diversity_phylo_resource_slope_null.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()



def plot_resources_vs_slope_all_distances():

    distances = numpy.logspace(-3, 1, num=40)

    fig, ax = plt.subplots(figsize=(4,4))

    for d_idx, d in enumerate(distances):

        slopes_to_plot = []

        for n_resource in n_resource_all:

            if d in diversity_vs_diversity_phylo_resource_dict[n_resource]:
                slopes_to_plot.append(diversity_vs_diversity_phylo_resource_dict[n_resource][d]['slope'])

        if len(slopes_to_plot) < len(n_resource_all):
            continue


        ax.plot(n_resource_all, slopes_to_plot, c=diversity_utils.rgb_blue_phylo(d_idx), lw=1.5, linestyle='-', marker='o', zorder=2)


    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%sresources_vs_slope_all_distances.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()





plot_slopes()


#plot_resources_vs_slope_all_distances()
