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



diversity_vs_diversity_phylo_emp_template = config.data_directory + "diversity_vs_diversity_phylo.pickle"
environments_to_keep = diversity_utils.environments_to_keep



dbd_utils.make_diversity_vs_diversity_phylo_dict(rarefied = False)




def load_diversity_vs_diversity_phylo_emp_dict():


    with open(diversity_vs_diversity_phylo_emp_template, 'rb') as handle:
        dict_ = pickle.load(handle)

    return dict_





def make_diversity_vs_diversity_phylo_emp_slope_null():

    fig = plt.figure(figsize = (10, 10)) #
    fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

    environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
    for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

        for environment_idx, environment in enumerate(environment_chunk):

            diversity_vs_diversity_phylo_emp_dict = load_diversity_vs_diversity_phylo_emp_dict(environment, rarefied=False)
            diversity_species = diversity_vs_diversity_phylo_emp_dict['diversity_species']

            distances = list(diversity_vs_diversity_phylo_emp_dict.keys())
            distances.remove('diversity_species')
            distances.sort()

            lower_ci_all = [diversity_vs_diversity_phylo_emp_dict[d]['lower_ci'] for d in distances]
            upper_ci_all = [diversity_vs_diversity_phylo_emp_dict[d]['upper_ci'] for d in distances]
            slope_all = [diversity_vs_diversity_phylo_emp_dict[d]['slope'] for d in distances]

            lower_ci_all = numpy.asarray(lower_ci_all)
            upper_ci_all = numpy.asarray(upper_ci_all)
            slope_all = numpy.asarray(slope_all)

            ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

            ax.plot(distances, slope_all, ls='-', lw=1.5, c='k',  zorder=2)
            ax.scatter(distances, slope_all, c='k', label='Observed',  zorder=3)
            ax.fill_between(distances, lower_ci_all, upper_ci_all, color='darkgrey', alpha=0.7, label='95% CI', zorder=1)

            #ax.set_xticks(idx_taxa_tanks)
            #ax.set_xticklabels(taxa_ranks_format, fontsize=9)
            ax.set_ylabel('Fine vs. coarse-grained diversity slope', fontsize=9)
            #ax.set_title('%d resources' % n_resource, fontsize=11)
            #ax.set_xlim([0,4])
            ax.set_xscale('log', base=10)
            ax.set_xlabel('Phylogenetic distance', fontsize=9)

            ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)




    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%sdiversity_vs_diversity_phylo_emp_slope_null.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()






def plot_diversity_vs_diversity():

    fig = plt.figure(figsize = (10, 20)) #
    fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

    diversity_vs_diversity_phylo_emp_dict = load_diversity_vs_diversity_phylo_emp_dict()

    distance_all = numpy.asarray([0.01663614249384222, 0.02767612370754228, 0.05938601867590266, 0.16435747993726982, 0.3526699214174661])


    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        diversity_species = diversity_vs_diversity_phylo_emp_dict[environment]['diversity_species']

        for distance_idx, distance in enumerate(distance_all):

            ax = plt.subplot2grid((9, 5), (environment_idx, distance_idx))

            diversity_rank = diversity_vs_diversity_phylo_emp_dict[environment][distance]['diversity']


            min_ = min(numpy.concatenate([diversity_species, diversity_rank]))
            max_ = max(numpy.concatenate([diversity_species, diversity_rank]))
            ax.plot([min_,max_],[min_,max_], lw=3,ls=':',c='k',zorder=1, label='1:1')

            ax.set_xlim([min_,max_])
            ax.set_ylim([min_,max_])

            ax.scatter(diversity_species, diversity_rank, alpha=0.3, c='k', zorder=2)

            if (distance_idx == 0) and (environment_idx == 0):
                ax.legend(loc="upper left", fontsize=7)

            if distance_idx == 0:
                ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=12)
            #ax.set_ylabel(rank.capitalize(), fontsize=12)

            if environment_idx == 0:
                ax.set_title('Phylo. dist. = %s' % str(round(distance, 3)), fontsize=11)


    fig.text(0.03, 0.5, "Coarse-grained Shannon's diversity", va='center', rotation='vertical', fontsize=18)
    fig.text(0.4, 0.08, "Shannon's diversity", va='center', fontsize=18)


    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%sdiversity_vs_diversity_phylo.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()









plot_diversity_vs_diversity()


#dbd_utils.make_coarse_grained_tree_dict('marine sediment metagenome', rarefied=True)



#dbd_utils.run_all_environments()

#dbd_utils.make_diversity_vs_diversity_taxon_dalbello_dict()



#print(diversity_vs_diversity_taxon_resource_dict)
