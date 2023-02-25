import os
#from Bio import Phylo
import config
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
import itertools

import dbd_utils
import diversity_utils
import mle_utils
import plot_utils



iter = 1000
rarefied = False

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label
color_all = numpy.asarray([plot_utils.rgb_blue_taxon(r) for r in idx_taxa_ranks])


error_null_dict_path = config.data_directory + "error_null_phylo_dict.pickle"


coarse_grained_tree_dict_all = {}
distances_all = []
for environment in diversity_utils.environments_to_keep:
    sys.stderr.write("Loading tree dict for %s...\n" % environment)
    coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_dict(environment=environment, rarefied=rarefied)
    coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict


def make_error_null_dict():

    error_null_dict = {}

    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        distances = list(coarse_grained_tree_dict.keys())
        distances.sort()

        error_null_dict[environment] = {}
        for distance_idx, distance in enumerate(distances):

            sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
            coarse_grained_list = coarse_grained_tree_dict[distance]
            coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

            # get indexes for each clade
            coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

            # coarse grain s-by-s for all clades
            s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)

            occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s_all_clades, s_by_s_all_clades.shape[0])
            error = diversity_utils.calculate_error(occupancies, predicted_occupancies)
            median_log_error = numpy.mean(numpy.log10(error))

            # get indices for null genera coarse graining
            coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(coarse_grained_n))[:-1]

            s_by_s_copy = numpy.copy(s_by_s)
            median_log_error_null_all = []
            for i in range(iter):

                print(i)

                numpy.random.shuffle(s_by_s_copy)
                s_by_s_coarse_null = numpy.add.reduceat(s_by_s_copy, coarse_grain_by_genus_idx, axis=0)
                s_by_s_coarse_null = s_by_s_coarse_null[:,~(numpy.all(s_by_s_coarse_null == 0, axis=0))]

                occupancies_null, predicted_occupancies_null, mad_null, beta_null, species_null = diversity_utils.predict_occupancy(s_by_s_coarse_null, s_by_s_coarse_null.shape[0])
                error_null = diversity_utils.calculate_error(occupancies_null, predicted_occupancies_null)
                median_log_error_null = numpy.mean(numpy.log10(error_null))
                median_log_error_null_all.append(median_log_error_null)


            median_log_error_null_all = numpy.asarray(median_log_error_null_all)
            median_log_error_null_all = numpy.sort(median_log_error_null_all)

            lower_ci = median_log_error_null_all[int(iter*0.025)]
            upper_ci = median_log_error_null_all[int(iter*0.975)]

            error_null_dict[environment][distance] = {}
            error_null_dict[environment][distance]['median_log_error'] = median_log_error
            error_null_dict[environment][distance]['median_log_error_lower_ci'] = lower_ci
            error_null_dict[environment][distance]['median_log_error_upper_ci'] = upper_ci



    sys.stderr.write("Saving dictionary...\n")
    with open(error_null_dict_path, 'wb') as handle:
        pickle.dump(error_null_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




make_error_null_dict()

def load_error_null_dict():

    with open(error_null_dict_path, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_



error_null_dict = load_error_null_dict()




fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environment_chunk_all = [diversity_utils.environments_to_keep[x:x+3] for x in range(0, len(diversity_utils.environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        distances = list(error_null_dict[environment].keys())
        distances.sort()

        median_log_error = numpy.asarray([error_null_dict[environment][r]['median_log_error'] for r in distances])
        lower_ci = numpy.asarray([error_null_dict[environment][r]['median_log_error_lower_ci'] for r in distances])
        upper_ci = numpy.asarray([error_null_dict[environment][r]['median_log_error_upper_ci'] for r in distances])


        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        ax.plot(distances, median_log_error, ls='-', lw=1.5, c='k',  zorder=2)
        ax.scatter(distances, median_log_error, c=color_all, s=60, label='Observed', linewidth=0.8, edgecolors='k', zorder=3)
        ax.fill_between(distances, lower_ci, upper_ci, color='darkgrey', alpha=0.7, label='95% Confidence Interval', zorder=1)

        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)
        #ax.set_xticks(idx_taxa_ranks)
        #ax.set_xticklabels(taxa_ranks_label, fontsize=8)
        ax.set_xscale('log', base=10)


        if environment_idx == 0:
            ax.set_ylabel("Median " + r'$\mathrm{log}_{10}$' +  " occupancy error", fontsize = 10)

        if environment_chunk_idx == len(environment_chunk_all)-1:
            ax.set_xlabel("Phylogenetic distance", fontsize = 10)


        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="lower left", fontsize=7)





fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%smedian_log_error_gamma_phylo.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
