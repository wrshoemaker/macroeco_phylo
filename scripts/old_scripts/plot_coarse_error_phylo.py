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
import simulation_utils
import mle_utils

import config



rarefied = True
iter = 1000

fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]


error_dict_path_template = config.data_directory + "error_dict_phylo%s.pickle"

#dbd_utils.make_coarse_grained_tree_dict()


def make_error_dict():

    error_dict = {}
    for environment in environments_to_keep:
        error_dict[environment] = {}

        coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_dict(environment=environment, rarefied=rarefied)

        distances = list(coarse_grained_tree_dict.keys())
        distances.sort()

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        sys.stderr.write("Getting site-by-species matrix...\n")
        s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)

        occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s, taxonomy_names)
        error = numpy.absolute(occupancies - predicted_occupancies)/occupancies

        idx_to_keep =  (~numpy.isinf(error)) & (error>0)
        species = species[idx_to_keep]
        error = error[idx_to_keep]

        for s_idx, s in enumerate(species):
            error_dict[environment][s] = {}


        sys.stderr.write("Running phylogenetic coarse-graining.....\n")
        for distance_idx, distance in enumerate(distances):

            sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
            coarse_grained_list = coarse_grained_tree_dict[distance]
            coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

            if len(coarse_grained_n) == len(taxonomy_names):
                continue

            # get indexes for each clade

            #coarse_idx = numpy.append([0], numpy.cumsum(counts_coarse))[:-1]
            coarse_grained_idx_all = [numpy.asarray([numpy.where(taxonomy_names==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            # coarse grain s-by-s for all clades
            s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)

            occupancies_clades, predicted_occupancies_clades, mad_clades, beta_clades, species_clades = diversity_utils.predict_occupancy(s_by_s_all_clades, list(range(len(coarse_grained_idx_all))))
            error_clades = numpy.absolute(occupancies_clades - predicted_occupancies_clades)/occupancies_clades

            #idx_to_keep =  (~numpy.isinf(error_clades)) & (error_clades>0)
            #error_clades = error_clades[idx_to_keep]
            #species_clades = species_clades[idx_to_keep]

            for e_idx, e in enumerate(error_clades):

                if (numpy.isinf(e)==True):
                    continue

                coarse_grained_idx = coarse_grained_idx_all[e_idx]

                for t_idx in coarse_grained_idx:

                    taxonomy_name_t = taxonomy_names[t_idx]
                    if taxonomy_name_t in error_dict[environment]:
                        error_dict[environment][taxonomy_name_t][distance] = e


    error_dict_path = error_dict_path_template % (diversity_utils.get_rarefied_label(rarefied))

    with open(error_dict_path, 'wb') as handle:
        pickle.dump(error_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_error_dict(rarefied):

    error_dict_path = error_dict_path_template % (diversity_utils.get_rarefied_label(rarefied))

    with open(error_dict_path, 'rb') as handle:
        error_dict = pickle.load(handle)

    return error_dict




#make_error_dict()
error_dict = load_error_dict(rarefied)


environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        print(environment)

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        species = list(error_dict[environment].keys())
        #ranks_in_dict = ['otu', 'genus', 'family', 'order',  'class',  'phylum']

        error_distance_dict = {}
        for s_idx, s in enumerate(species):

            #if s_idx > 100:
            #    continue

            #error_to_plot = [error_dict[s][r] for r in ranks_in_dict if r in error_dict[s]]
            distances_to_plot = []
            error_to_plot = []

            distances = list(error_dict[environment][s].keys())
            distances.sort()

            for d_idx, d in enumerate(distances):

                error_d = error_dict[environment][s][d]

                if error_d == 0:
                    continue

                distances_to_plot.append(d)
                error_to_plot.append(error_d)

                if d not in error_distance_dict:
                    error_distance_dict[d] = []

                error_distance_dict[d].append(error_d)

            #if len(error_to_plot) < 10:
            #    continue

            ax.plot(distances_to_plot, error_to_plot, ls='-', lw=0.3, c='dodgerblue', alpha=0.003,  zorder=2)


        distances = numpy.asarray(list(error_distance_dict.keys()))
        distances = numpy.sort(distances)

        if environment == 'soil metagenome':
            print(distances)

        #distances.sort()

        #mean_error = numpy.asarray([numpy.mean(numpy.log10(error_distance_dict[d])) for d in distances])
        mean_error = numpy.asarray([numpy.mean(error_distance_dict[d]) for d in distances])
        n_ = numpy.asarray([numpy.mean(len(error_distance_dict[d])) for d in distances])


        #idx_to_keep = (n_ >= 20)
        #distances = distances[idx_to_keep]
        #mean_error = mean_error[idx_to_keep]

        ax.plot(distances, 10**mean_error, ls='-', lw=1.5, c='k', zorder=3)
        ax.scatter(distances, 10**mean_error, c='k', zorder=4)

        ax.set_xlim([min(distances), max(distances)])
        ax.set_ylim([min(10**mean_error)*0.1, max(10**mean_error)*10])

        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)

        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)


        if environment_idx == 0:
            ax.set_ylabel("Occupancy prediction error, gamma", fontsize = 10)

        if environment_chunk_idx == len(environment_chunk_all)-1:
            ax.set_xlabel("Phylogenetic distance", fontsize = 10)

        #if (environment_chunk_idx ==0) and (environment_idx == 0):
        #    ax.legend(loc="upper left", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%serror_coarse_emp_phylo%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
