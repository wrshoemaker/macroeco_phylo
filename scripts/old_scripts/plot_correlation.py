import os
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
import plot_utils


rarefied = False

environments_to_keep = diversity_utils.environments_to_keep




def test_symmetry_ks(array_, iter=1000):

    neg_array_ = -1*array_

    D, p = stats.ks_2samp(array_, neg_array_)

    array_merged = numpy.concatenate((array_, neg_array_))

    D_null_all = []
    for i in range(iter):
        numpy.random.shuffle(array_merged)

        array_null = array_merged[:len(array_)]
        neg_array_null = array_merged[len(array_):]

        D_null, p_null = stats.ks_2samp(array_null, neg_array_null)
        D_null_all.append(D_null)

    D_null_all = numpy.asarray(D_null_all)
    D_null_all = numpy.sort(D_null_all)

    p_perm = sum(D_null_all > D)/iter

    #lower_ci = D_null_all[int(iter*0.025)]
    #upper_ci = D_null_all[int(iter*0.975)]

    return D, p_perm


def symmetry_ks_bootstrap(array_, iter=10000):

    neg_array_ = -1*array_

    D, p = stats.ks_2samp(array_, neg_array_)

    D_bs_all = []
    for i in range(iter):

        array_boostrap = numpy.random.choice(array_, replace=True, size=len(array_))
        neg_array_bootstrap = -1*array_boostrap

        D_bs, p_bs = stats.ks_2samp(array_null, neg_array_null)
        D_bs_all.append(D_bs)

    D_bs_all = numpy.asarray(D_bs_all)
    D_bs_all = numpy.sort(D_bs_all)

    lower_ci = D_bs_all[int(iter*0.025)]
    upper_ci = D_bs_all[int(iter*0.975)]

    return D, lower_ci, upper_ci







def make_symmetry_rho_ks_dict():

    for environment in environments_to_keep:

        sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

        rho = numpy.cov(rel_s_by_s)
        rho_flat = rho[numpy.triu_indices(rho.shape[0], k = 1)]

        #D, p_perm = test_symmetry_ks(rho_flat)

        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)
            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))

            genus_to_taxa_dict = {}
            sad_genera_all = []
            for genus_idx, genus in enumerate(all_genera):
                genus_to_taxa_dict[genus] = []
                for t in taxa:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])

                if len(g_taxa_idx) < 5:
                    continue


                rho_g = rho[numpy.ix_(g_taxa_idx, g_taxa_idx)]


                rho_g_flat = rho_g[numpy.triu_indices(rho_g.shape[0], k = 1)]

                #rho_g_flat_neg = -1*rho_g_flat


                D, p_perm = test_symmetry_ks(rho_g_flat)

                #print(rho[g_taxa_idx,:])

                print(D, p_perm)

                #print(g_taxa_idx)


make_symmetry_rho_ks_dict()

#make_symmetry_rho_ks_dict()
