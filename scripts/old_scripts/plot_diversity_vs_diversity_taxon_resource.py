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


dbd_utils.make_diversity_vs_diversity_taxon_dalbello_dict()


diversity_vs_diversity_taxon_dict = dbd_utils.load_diversity_vs_diversity_taxon_resource_dict()
n_resource_all = list(diversity_vs_diversity_taxon_dict.keys())
n_resource_all.sort()

taxa_ranks = diversity_utils.taxa_ranks[::-1]
taxa_ranks_format = [x.capitalize() for x in taxa_ranks]
idx_taxa_tanks = range(len(taxa_ranks))




def plot_diversity_vs_diversity():

    fig = plt.figure(figsize = (20, 10)) #
    fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

    for n_resource_idx, n_resource in enumerate(n_resource_all):


        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            ax = plt.subplot2grid((len(diversity_utils.taxa_ranks), len(n_resource_all)), (rank_idx, n_resource_idx))

            diversity_all = diversity_vs_diversity_taxon_dict[n_resource][rank]['diversity_all']
            diversity_coarse = diversity_vs_diversity_taxon_dict[n_resource][rank]['diversity_coarse']

            ax.scatter(diversity_all, diversity_coarse, alpha=0.3, c='k')

            if rank_idx == 0:
                ax.set_title('%d resources' % n_resource, fontsize=12)

            if n_resource_idx == 0:
                ax.set_ylabel(rank.capitalize(), fontsize=12)


    fig.text(0.07, 0.5, "Coarse-grained Shannon's diversity", va='center', rotation='vertical', fontsize=18)
    fig.text(0.5, 0.06, "Shannon's diversity", va='center', fontsize=18)

    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%sdiversity_vs_diversity_taxon_resource.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()




def plot_slopes():


    fig = plt.figure(figsize = (10, 10)) #
    fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

    n_resource_chunks = [n_resource_all[x:x+3] for x in range(0, len(n_resource_all), 3)]
    for n_resource_chunk_idx, n_resource_chunk in enumerate(n_resource_chunks):

        for n_resource_idx, n_resource in enumerate(n_resource_chunk):

            lower_ci_all = [diversity_vs_diversity_taxon_dict[n_resource][rank]['lower_ci'] for rank in taxa_ranks]
            upper_ci_all = [diversity_vs_diversity_taxon_dict[n_resource][rank]['upper_ci'] for rank in taxa_ranks]
            slope_all = [diversity_vs_diversity_taxon_dict[n_resource][rank]['slope'] for rank in taxa_ranks]

            lower_ci_all = numpy.asarray(lower_ci_all)
            upper_ci_all = numpy.asarray(upper_ci_all)
            slope_all = numpy.asarray(slope_all)

            ax = plt.subplot2grid((3, 3), (n_resource_chunk_idx, n_resource_idx))

            ax.plot(idx_taxa_tanks, slope_all, ls='-', lw=1.5, c='k',  zorder=2)
            ax.scatter(idx_taxa_tanks, slope_all, c='k', label='Observed',  zorder=3)
            ax.fill_between(idx_taxa_tanks, lower_ci_all, upper_ci_all, color='darkgrey', alpha=0.7, label='95% CI', zorder=1)

            ax.set_xticks(idx_taxa_tanks)
            ax.set_xticklabels(taxa_ranks_format, fontsize=9)
            ax.set_ylabel('Fine vs. coarse-grained diversity slope', fontsize=9)
            ax.set_title('%d resources' % n_resource, fontsize=11)
            ax.set_xlim([0,4])


            if (n_resource_chunk_idx==0) and (n_resource_idx==0):
                ax.legend(loc="lower left", fontsize=7)


    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%sdiversity_vs_diversity_taxon_resource_slope_null.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()



plot_diversity_vs_diversity()
plot_slopes()
