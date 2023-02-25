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
import mle_utils

import predict_diversity_slm_emp

from scipy.special import erf


import simulation_utils



diversity_slm_taxon_dict_template = config.data_directory + "diversity_slm_taxon%s.pickle"


rarefied = False
iter = 1000
#environment = 'human gut metagenome'

def make_diversity_slm_prediction_taxon_dict():

    diversity_dict = {}
    for environment in [diversity_utils.environments_to_keep[0]]:

        parameter_dict = predict_diversity_slm_emp.load_slm_parameter_dict(environment)

        # get best parameter
        idx_ = parameter_dict['mean_euc_distance'].index(min(parameter_dict['mean_euc_distance']))
        c = parameter_dict['c'][idx_]
        sigma = parameter_dict['sigma'][idx_]

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #samples = diversity_utils.subset_observations(environment=environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        sys.stderr.write("Getting site-by-species matrix...\n")
        #s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        diversity_species = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s)

        diversity_dict[environment] = {}
        diversity_dict[environment]['asv'] = {}
        diversity_dict[environment]['asv']['observed'] = numpy.mean(diversity_species)
        diversity_dict[environment]['asv']['null'] = []

        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))

            genus_to_taxa_dict = {}
            sad_genera_all = []
            for genus in all_genera:
                genus_to_taxa_dict[genus] = []
                for t in taxa:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
                sad_genera_all.append(s_by_s[g_taxa_idx,:].sum(axis=0))


            s_by_s_genera = numpy.stack(sad_genera_all, axis=0)

            diversity_species_rank = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_genera)

            unique_genera, counts_genera = numpy.unique([sad_annotated_dict['taxa'][t][rank] for t in taxa], return_counts=True)
            coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_genera))[:-1]

            diversity_dict[environment][rank] = {}
            diversity_dict[environment][rank]['observed'] = numpy.mean(diversity_species_rank)
            diversity_dict[environment][rank]['null'] = []
            diversity_dict[environment][rank]['idx'] = coarse_grain_by_genus_idx


        #mu=-16
        #s=3
        rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
        mean = numpy.mean(rel_s_by_s, axis=1)
        #c = 0.00005
        mu, s = diversity_utils.Klogn(mean, c = c)

        N = s_by_s.sum(axis=0)
        #S = s_by_s.shape[0]
        S_obs = s_by_s.shape[0]
        S_tot = int(2*S_obs / erf((numpy.log(c)-mu) / numpy.sqrt(2)*s))

        count = 0
        #for i in range(iter):
        while count < iter:

            #if count %100 == 0:
            #    print(count)

            s_by_s_null = simulation_utils.generate_community(mu, s, S_tot, N, 'constant', sigma, len(N))
            # we cant coarse grain more species than we observed

            n_non_zero = sum(~numpy.all(s_by_s_null == 0, axis=1))
            if n_non_zero > S_obs:
                continue

            # add "missing" species back in
            if n_non_zero < S_obs:

                s_by_s_null_no_zeros = s_by_s_null[~(numpy.all(s_by_s_null == 0, axis=1))]
                n_zero = S_obs - s_by_s_null_no_zeros.shape[0]
                # put zeros back in to simulated number of species is equal to observed number
                s_by_s_zeros = numpy.zeros((n_zero, s_by_s_null_no_zeros.shape[1]))
                s_by_s_null = numpy.vstack((s_by_s_null_no_zeros, s_by_s_zeros))


            numpy.random.shuffle(s_by_s_null)

            diversity_species_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_null)
            diversity_dict[environment]['asv']['null'].append(numpy.mean(diversity_species_null))

            for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

                s_by_s_coarse_null = numpy.add.reduceat(s_by_s_null, diversity_dict[environment][rank]['idx'], axis=0)

                diversity_species_coarse_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_coarse_null)
                diversity_dict[environment][rank]['null'].append(numpy.mean(diversity_species_coarse_null))

            count += 1


        for rank in diversity_dict[environment].keys():

            null = numpy.asarray(diversity_dict[environment][rank]['null'])
            null = numpy.sort(null)

            lower_ci = null[int(iter*0.025)]
            upper_ci = null[int(iter*0.975)]

            diversity_dict[environment][rank]['lower_ci'] = lower_ci
            diversity_dict[environment][rank]['upper_ci'] = upper_ci

            del diversity_dict[environment][rank]['null']

            if rank != 'asv':

                del diversity_dict[environment][rank]['idx']

            print(diversity_dict[environment][rank]['observed'], lower_ci, upper_ci)



    diversity_slm_taxon_dict_path = diversity_slm_taxon_dict_template % diversity_utils.get_rarefied_label(rarefied)
    sys.stderr.write("Saving diversity dictionary...\n")
    with open(diversity_slm_taxon_dict_path, 'wb') as handle:
        pickle.dump(diversity_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




make_diversity_slm_prediction_taxon_dict()

def load_diversity_slope_phylo_dict():

    diversity_slm_taxon_dict_path = diversity_slm_taxon_dict_template % diversity_utils.get_rarefied_label(rarefied)

    with open(diversity_slm_taxon_dict_path, 'rb') as handle:
        diversity_slm_taxon_dict = pickle.load(handle)
    return diversity_slm_taxon_dict







def plot_diversity_prediction():


    diversity_slm_taxon_dict = load_diversity_slope_phylo_dict()

    #make_diversity_slm_taxon_dict()

    fig, ax = plt.subplots(figsize=(4,4))

    environment = 'human gut metagenome'

    taxa_ranks = diversity_utils.taxa_ranks
    idx_taxa_ranks = list(range(len(taxa_ranks)))
    idx_taxa_ranks.append(len(taxa_ranks))
    taxa_ranks_label = diversity_utils.taxa_ranks_label
    taxa_ranks_label.insert(0,'ASV')


    ranks_in_dict = ['asv', 'genus', 'family', 'order',  'class',  'phylum']

    observed_all = []
    lower_ci_all = []
    upper_ci_all = []
    for r in ranks_in_dict:
        observed_all.append(diversity_slm_taxon_dict[environment][r]['observed'])
        lower_ci_all.append(diversity_slm_taxon_dict[environment][r]['lower_ci'])
        upper_ci_all.append(diversity_slm_taxon_dict[environment][r]['upper_ci'])



    ax.plot(idx_taxa_ranks, observed_all, ls='-', lw=1.5, c='k',  zorder=2)
    ax.scatter(idx_taxa_ranks, observed_all, c='k', label='Observed',  zorder=3)
    ax.fill_between(idx_taxa_ranks, lower_ci_all, upper_ci_all, color='darkgrey', alpha=0.7, label='SLM 95% CI', zorder=1)

    ax.set_xticks(idx_taxa_ranks)
    ax.set_xticklabels(taxa_ranks_label, fontsize=9)
    ax.set_ylabel("Mean shannon's diversity", fontsize=9)
    #ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)
    #ax.set_xlim([0,4])


    #ax.hist(diversity_species, 20, density=True, histtype='step', color='b', alpha=0.75, label = 'Observed')
    #ax.hist(diversity_species_sim, 20, density=True, histtype='step', color='r', alpha=0.75, label = 'Simulated')


    #ax.set_xlabel('Diversity', fontsize=9)

    #ax.set_ylabel('Probability density', fontsize=9)

    ax.legend(loc="upper right", fontsize=7)

    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%sdiversity_slm_taxon.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()




    #sys.stderr.write("Running taxonomic coarse-graining.....\n")
    #for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

    #    sys.stderr.write("Starting %s level analysis...\n" % rank)

    #    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    #    all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
    #    all_genera_idx = numpy.arange(len(all_genera))

    #    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
    #    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

    #    genus_to_taxa_dict = {}
    #    sad_genera_all = []

    #    beta_div_mean_ratio_all = []
    #    for genus in all_genera:
    #        genus_to_taxa_dict[genus] = []
    #        for t in taxa:
    #            if sad_annotated_dict['taxa'][t][rank] == genus:
    #                genus_to_taxa_dict[genus].append(t)

    #        g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
    #        sad_genera_all.append(s_by_s[g_taxa_idx,:].sum(axis=0))


    #    s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
    #    # remove sites where there are no observations
    #    s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]

    #    diversity_genera = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_genera)
