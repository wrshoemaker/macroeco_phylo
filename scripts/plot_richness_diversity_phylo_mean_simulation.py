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

from scipy.special import polygamma
from scipy.stats import norm
from scipy.stats import gamma
#import dbd_utils
#import plot_utils

#import diversity_utils
import config
import scipy.integrate as integrate
import scipy.special as special

import diversity_utils
import dbd_utils
import simulation_utils


#dbd_utils.make_svd_phylo_dict()

#for environment in diversity_utils.environments_to_keep:
#dbd_utils.richness_diversity_phylo_mean_simulation('soil metagenome')

#    dbd_utils.make_diversity_slm_phlyo_integral_otu_level_dict(environment)


def test_():
    environment =  'marine metagenome'

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = False)

    samples = sad_annotated_dict['samples']
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    sys.stderr.write("Getting site-by-species matrix...\n")

    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
    rel_s_by_s = s_by_s/numpy.sum(s_by_s, axis=0)
    mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
    n_reads = numpy.sum(s_by_s, axis=0)


    #a, reads_svd_ = simulation_utils.genrate_community_from_mean_and_var_svd(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads), svd_null)
    #print(f"Run time {monotonic() - start_time} seconds")
    #print(time.clock() - start_time, "seconds")


    # rvs
    #start_time = time.clock()
    #start_time = monotonic()
    a, reads_rv_ = simulation_utils.genrate_community_from_mean_and_var_using_rvs(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads))
    #print(f"Run time {monotonic() - start_time} seconds")
    #print(time.clock() - start_time, "seconds")

    #diversity_svd = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, reads_svd_)
    diversity_rv = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, reads_rv_)
    diversity_observed = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s)


    print(numpy.mean(diversity_observed), numpy.mean(diversity_rv))




    taxa_ranks = ['family']

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


        rel_s_by_s_genera = s_by_s_genera/numpy.sum(s_by_s_genera, axis=0)
        mean_rel_s_by_s_genera = numpy.mean(rel_s_by_s_genera, axis=1)
        var_rel_s_by_s_genera = numpy.var(rel_s_by_s_genera, axis=1)
        n_reads = numpy.sum(s_by_s_genera, axis=0)

        diversity_observed = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_genera)


        # test inverse transform sampling
        u_null = numpy.identity(len(mean_rel_s_by_s_genera))
        s_null = numpy.array([1]*len(mean_rel_s_by_s_genera))
        v_null = numpy.identity(len(mean_rel_s_by_s_genera))
        svd_null = (u_null, s_null, v_null)

        import time

        for i in range(10):

            start_time = time.clock()
            #start_time = monotonic()
            #print(mean_rel_s_by_s_genera.shape, var_rel_s_by_s_genera.shape, n_reads.shape)
            a, reads_svd = simulation_utils.genrate_community_from_mean_and_var_svd(mean_rel_s_by_s_genera, var_rel_s_by_s_genera, n_reads, len(n_reads), svd_null)
            #print(f"Run time {monotonic() - start_time} seconds")
            #print(time.clock() - start_time, "seconds")


            # rvs
            start_time = time.clock()
            #start_time = monotonic()
            a, reads_rv = simulation_utils.genrate_community_from_mean_and_var_using_rvs(mean_rel_s_by_s_genera, var_rel_s_by_s_genera, n_reads, len(n_reads))
            #print(f"Run time {monotonic() - start_time} seconds")
            #print(time.clock() - start_time, "seconds")

            rel_reads_svd = reads_svd/numpy.sum(reads_svd, axis=0)
            rel_reads_rv = reads_rv/numpy.sum(reads_rv, axis=0)

            diversity_svd = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_reads_svd)
            diversity_rv = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_reads_rv)


            print(numpy.mean(diversity_observed), numpy.mean(diversity_svd), numpy.mean(diversity_rv))



test_()
