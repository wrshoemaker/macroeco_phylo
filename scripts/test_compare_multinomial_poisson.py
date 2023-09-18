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
import plot_utils

import diversity_utils
import config
import simulation_utils




def predict_second_moment_diversity(s_by_s, iter_=100):

    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

    n_reads = s_by_s.sum(axis=0)

    mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
    beta_rel_s_by_s = (mean_rel_s_by_s**2)/var_rel_s_by_s


    u = numpy.identity(len(mean_rel_s_by_s))
    s = numpy.asarray([1]*len(mean_rel_s_by_s))
    v = numpy.identity(len(mean_rel_s_by_s))
    svd = (u, s, v)

    expected_diversity_second_moment_rvs_all = []
    product_pairs_null = []
    var_diversity_rvs_multinomial_all = []
    var_diversity_rvs_poisson_all = []
    for i in range(iter_):

        rel_abundances_gamma, read_counts_gamma_multinomial, read_counts_gamma_poisson = simulation_utils.genrate_community_from_mean_and_var_svd_multinomial_and_poisson(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads), svd)

        rel_read_counts_gamma_multinomial = read_counts_gamma_multinomial/numpy.sum(read_counts_gamma_multinomial, axis=0)
        rel_read_counts_gamma_poisson = read_counts_gamma_poisson/numpy.sum(read_counts_gamma_poisson, axis=0)

        var_diversity_rvs_multinomial = numpy.var(numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_read_counts_gamma_multinomial))
        var_diversity_rvs_multinomial_all.append(var_diversity_rvs_multinomial)

        var_diversity_rvs_poisson = numpy.var(numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_read_counts_gamma_poisson))
        var_diversity_rvs_poisson_all.append(var_diversity_rvs_poisson)

    print(rank, numpy.mean(var_diversity_rvs_multinomial_all), numpy.mean(var_diversity_rvs_poisson_all))





environment =  'marine metagenome'

sys.stderr.write("Subsetting samples for %s...\n" % environment)
sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = False)

samples = sad_annotated_dict['samples']
taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
sys.stderr.write("Getting site-by-species matrix...\n")

s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

#taxa_ranks = ['genus', 'family', 'order', 'class']
taxa_ranks = ['genus', 'family', 'order']


for rank_idx, rank in enumerate(taxa_ranks):

    if rank == 'OTU':

        s_by_s_genera = s_by_s

    else:

        sys.stderr.write("Running taxonomic coarse-graining.....\n")

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

        predict_second_moment_diversity(s_by_s_genera)
    
