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


error_null_taxon_dict_path = config.data_directory + "error_null_taxon_dict.pickle"



def make_error_null_dict():

    error_null_dict = {}

    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        error_null_dict[environment] = {}

        sys.stderr.write("Running taxonomic coarse-graining.....\n")
        for rank_idx, rank in enumerate(taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)
            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))

            genus_to_taxa_dict = {}
            sad_genera_all = []
            for genus_idx, genus in enumerate(all_genera):
                genus_to_taxa_dict[genus] = []
                var_asv = 0
                for t in taxa:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
                s_by_s_genus = s_by_s[g_taxa_idx,:]
                sad_genera_all.append(s_by_s_genus.sum(axis=0))


            s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
            # remove sites where there are no observations
            s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]

            # get indices for null genera coarse graining
            unique_genera, counts_genera = numpy.unique([sad_annotated_dict['taxa'][t][rank] for t in taxa], return_counts=True)
            coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_genera))[:-1]

            occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s_genera, s_by_s_genera.shape[0])
            error = diversity_utils.calculate_error(occupancies, predicted_occupancies)
            median_log_error = numpy.mean(numpy.log10(error))

            s_by_s_copy = numpy.copy(s_by_s)
            median_log_error_null_all = []
            for i in range(iter):

                numpy.random.shuffle(s_by_s_copy)
                #print(coarse_grain_by_genus_idx)
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

            error_null_dict[environment][rank] = {}
            error_null_dict[environment][rank]['median_log_error'] = median_log_error
            error_null_dict[environment][rank]['median_log_error_lower_ci'] = lower_ci
            error_null_dict[environment][rank]['median_log_error_upper_ci'] = upper_ci



    sys.stderr.write("Saving dictionary...\n")
    with open(error_null_taxon_dict_path, 'wb') as handle:
        pickle.dump(error_null_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)





def load_error_null_dict():

    with open(error_null_taxon_dict_path, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_




make_error_null_dict()
error_null_dict = load_error_null_dict()




fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environment_chunk_all = [diversity_utils.environments_to_keep[x:x+3] for x in range(0, len(diversity_utils.environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        median_log_error = numpy.asarray([error_null_dict[environment][r]['median_log_error'] for r in taxa_ranks])
        lower_ci = numpy.asarray([error_null_dict[environment][r]['median_log_error_lower_ci'] for r in taxa_ranks])
        upper_ci = numpy.asarray([error_null_dict[environment][r]['median_log_error_upper_ci'] for r in taxa_ranks])


        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))


        ax.plot(idx_taxa_ranks, median_log_error, ls='-', lw=1.5, c='k',  zorder=2)
        ax.scatter(idx_taxa_ranks, median_log_error, c=color_all, s=60, label='Observed', linewidth=0.8, edgecolors='k', zorder=3)
        ax.fill_between(idx_taxa_ranks, lower_ci, upper_ci, color='darkgrey', alpha=0.7, label='95% Confidence Interval', zorder=1)

        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)
        ax.set_xticks(idx_taxa_ranks)
        ax.set_xticklabels(taxa_ranks_label, fontsize=8)


        if environment_idx == 0:
            ax.set_ylabel("Median " + r'$\mathrm{log}_{10}$' +  " occupancy error", fontsize = 10)


        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="lower left", fontsize=7)





fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%smedian_log_error_gamma_taxon.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
