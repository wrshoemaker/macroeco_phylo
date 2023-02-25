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



iter_ = 100


for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = False)

    samples = sad_annotated_dict['samples']
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    sys.stderr.write("Getting site-by-species matrix...\n")
    #s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)

    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

    sys.stderr.write("Running taxonomic coarse-graining.....\n")


    coarse_graining_dict = {}

    for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

        print(rank)

        if rank == 'OTU':

            s_by_s_genera = s_by_s

        else:

            sys.stderr.write("Starting %s level analysis...\n" % rank)
            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))

            genus_to_taxa_dict = {}
            sad_genera_all = []
            g_taxa_idx_all = []
            n_taxa_all = []
            for genus_idx, genus in enumerate(all_genera):
                genus_to_taxa_dict[genus] = []
                var_asv = 0
                for t in taxa:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)


                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
                s_by_s_genus = s_by_s[g_taxa_idx,:]
                sad_genera_all.append(s_by_s_genus.sum(axis=0))

                g_taxa_idx_all.extend(g_taxa_idx.tolist())
                n_taxa_all.append(len(g_taxa_idx))


            s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
            # remove sites where there are no observations
            s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]

            g_taxa_idx_all = numpy.asarray(g_taxa_idx_all)
            n_taxa_all = numpy.asarray(n_taxa_all)

            coarse_graining_dict[rank] = {}
            coarse_graining_dict[rank]['s_by_s_sorted'] = s_by_s[g_taxa_idx_all,:]
            coarse_graining_dict[rank]['taxa_idx_all'] = g_taxa_idx_all
            coarse_graining_dict[rank]['n_taxa_all'] = n_taxa_all

            print(len(n_taxa_all), s_by_s.shape)


    
    # generate null matrices at OTU level
    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

    mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
    n_reads = numpy.sum(s_by_s, axis=0)

    read_counts_all = []

    for i in range(iter_):


        rel_abundances, read_counts = simulation_utils.genrate_community_from_mean_and_var(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads))
        read_counts_all.append(read_counts)

        print(read_counts)


        #if rank == 'OTU':




