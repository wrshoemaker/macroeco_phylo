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




#dbd_utils.make_diversity_vs_diversity_taxon_dict(rarefied = False)

diversity_vs_diversity_taxon_dict = dbd_utils.load_diversity_vs_diversity_taxon_dict(rarefied=False)
taxa_ranks = diversity_utils.taxa_ranks[::-1]
taxa_ranks_format = [x.capitalize() for x in taxa_ranks]
idx_taxa_tanks = range(len(taxa_ranks))

environments_to_keep = diversity_utils.environments_to_keep


def plot_slopes():

    fig = plt.figure(figsize = (5, 8)) #
    fig.subplots_adjust(bottom= 0.15,  wspace=0.25)

    ax_slope = plt.subplot2grid((2, 1), (0, 0), colspan=1)
    ax_percentile = plt.subplot2grid((2, 1), (1, 0), colspan=1)

    for environment in diversity_utils.environments_to_keep:

        slope_all = []
        percentile_all = []
        # start at genus
        for rank in taxa_ranks:

            slope_all.append(diversity_vs_diversity_taxon_dict[environment][rank]['slope'])
            percentile_all.append(diversity_vs_diversity_taxon_dict[environment][rank]['percentile'])

        ax_slope.plot(idx_taxa_tanks, slope_all, ls='-',  lw=1.5, alpha=0.9, c='k')
        ax_percentile.plot(idx_taxa_tanks, percentile_all, ls='-',  lw=1.5, alpha=0.9, c='k')



    ax_slope.set_xticks(idx_taxa_tanks)
    ax_slope.set_xticklabels(taxa_ranks)

    #ax_slope.set_xlabel('T of resources', fontsize=12)
    ax_slope.set_ylabel('Fine vs. coarse-grained diversity slope', fontsize=12)

    ax_percentile.set_xticks(idx_taxa_tanks)
    ax_percentile.set_xticklabels(taxa_ranks)

    #ax_percentile.set_xlabel('Number of resources', fontsize=12)
    ax_percentile.set_ylabel('Position of observed slope\nrelative to the null distribution', fontsize=12)

    #ax_slope.legend(loc="upper left", fontsize=9)
    #ax_percentile.legend(loc="upper left", fontsize=9)

    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%sdiversity_vs_diversity_taxon_slope.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()



#plot_slopes()

def plot_diversity_vs_diversity():

    fig = plt.figure(figsize = (10, 20)) #
    fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        diversity_species = diversity_vs_diversity_taxon_dict[environment]['diversity_species']

        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            ax = plt.subplot2grid((9, 5), (environment_idx, rank_idx))

            diversity_rank = diversity_vs_diversity_taxon_dict[environment][rank]['diversity']



            min_ = min(numpy.concatenate([diversity_species, diversity_rank]))
            max_ = max(numpy.concatenate([diversity_species, diversity_rank]))
            ax.plot([min_,max_],[min_,max_], lw=3,ls=':',c='k',zorder=1, label='1:1')

            ax.set_xlim([min_,max_])
            ax.set_ylim([min_,max_])

            ax.scatter(diversity_species, diversity_rank, alpha=0.3, c='k', zorder=2)

            if (rank_idx == 0) and (environment_idx == 0):
                ax.legend(loc="upper left", fontsize=7)

            if rank_idx == 0:
                ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=12)
            #ax.set_ylabel(rank.capitalize(), fontsize=12)

            if environment_idx == 0:
                ax.set_title(rank.capitalize(), fontsize=12)



    fig.text(0.03, 0.5, "Coarse-grained Shannon's diversity", va='center', rotation='vertical', fontsize=18)
    fig.text(0.4, 0.08, "Shannon's diversity", va='center', fontsize=18)


    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%sdiversity_vs_diversity_taxon.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()



def plot_diversity_vs_diversity_null():


    fig = plt.figure(figsize = (10, 10)) #
    fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

    environment_chunks = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
    for environment_chunk_idx, environment_chunk in enumerate(environment_chunks):

        for environment_idx, environment in enumerate(environment_chunk):

            lower_ci_all = [diversity_vs_diversity_taxon_dict[environment][rank]['lower_ci'] for rank in taxa_ranks]
            upper_ci_all = [diversity_vs_diversity_taxon_dict[environment][rank]['upper_ci'] for rank in taxa_ranks]
            slope_all = [diversity_vs_diversity_taxon_dict[environment][rank]['slope'] for rank in taxa_ranks]

            lower_ci_all = numpy.asarray(lower_ci_all)
            upper_ci_all = numpy.asarray(upper_ci_all)
            slope_all = numpy.asarray(slope_all)

            ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

            ax.plot(idx_taxa_tanks, slope_all, ls='-', lw=1.5, c='k',  zorder=2)
            ax.scatter(idx_taxa_tanks, slope_all, c='k', label='Observed',  zorder=3)
            ax.fill_between(idx_taxa_tanks, lower_ci_all, upper_ci_all, color='darkgrey', alpha=0.7, label='95% CI', zorder=1)

            ax.set_xticks(idx_taxa_tanks)
            ax.set_xticklabels(taxa_ranks_format, fontsize=9)
            ax.set_ylabel('Fine vs. coarse-grained diversity slope', fontsize=9)
            ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)
            ax.set_xlim([0,4])


            if (environment_chunk_idx==0) and (environment_idx==0):
                ax.legend(loc="lower left", fontsize=7)


    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%sdiversity_vs_diversity_slope_null.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()




plot_diversity_vs_diversity()
