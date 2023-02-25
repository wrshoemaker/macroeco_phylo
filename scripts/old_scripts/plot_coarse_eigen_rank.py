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






iter = 1000
environment = 'freshwater  metagenome'
rarefied=True

sys.stderr.write("Subsetting samples for %s...\n" % environment)
sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
samples = sad_annotated_dict['samples']

taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
#rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))


for rank in diversity_utils.taxa_ranks:

    sys.stderr.write("Starting %s level analysis...\n" % rank)
    all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
    all_genera_idx = numpy.arange(len(all_genera))

    genus_to_taxa_dict = {}
    sad_genera_all = []
    g_taxa_idx_all = []
    for genus in all_genera:
        genus_to_taxa_dict[genus] = []
        for t in taxa:
            if sad_annotated_dict['taxa'][t][rank] == genus:
                genus_to_taxa_dict[genus].append(t)

        g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
        g_taxa_idx_all.append(g_taxa_idx)
        sad_genera_all.append(s_by_s[g_taxa_idx,:].sum(axis=0))

    s_by_s_all_clades = numpy.stack(sad_genera_all, axis=0)
    # remove sites where there are no observations
    s_by_s_all_clades = s_by_s_all_clades[:,~(numpy.all(s_by_s_all_clades == 0, axis=0))]
    rel_s_by_s_all_clades = (s_by_s_all_clades/s_by_s_all_clades.sum(axis=0))

    rho_rel_s_by_s_all_clades = numpy.corrcoef(rel_s_by_s_all_clades)
    # eigendecomposition on a real symmetric matrix, returns real eigenvalues
    lambda_rho, v_rho = numpy.linalg.eigh(rho_rel_s_by_s_all_clades)
    lambda_rho = numpy.sort(lambda_rho)[::-1]

    # get indices for null genera coarse graining
    counts_genera = [len(n) for n in g_taxa_idx_all]
    coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_genera))[:-1]

    #idx_all = numpy.arange(len(coarse_grained_idx_all))
    s_by_s_copy = numpy.copy(s_by_s)
    max_all = []
    for i in range(iter):

        # diversity vs diversity null
        numpy.random.shuffle(s_by_s_copy)
        s_by_s_coarse_null = numpy.add.reduceat(s_by_s_copy, all_genera_idx, axis=0)
        s_by_s_coarse_null = s_by_s_coarse_null[:,~(numpy.all(s_by_s_coarse_null == 0, axis=0))]

        rel_s_by_s_coarse_null = (s_by_s_coarse_null/s_by_s_coarse_null.sum(axis=0))

        rho_rel_s_by_s_all_clades_null = numpy.corrcoef(rel_s_by_s_coarse_null)
        # eigendecomposition on a real symmetric matrix, returns real eigenvalues
        lambda_rho_null, v_rho_null = numpy.linalg.eigh(rho_rel_s_by_s_all_clades_null)
        lambda_rho_null = numpy.sort(lambda_rho_null)[::-1]

        max_all.append(lambda_rho_null[0])

    max_all = numpy.asarray(max_all)


    p_value = sum(max_all > lambda_rho[0])/iter

    print(rank, lambda_rho[0], p_value)
