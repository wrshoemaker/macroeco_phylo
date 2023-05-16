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



coarse_eigen_taxon_dict_path = "%scoarse_eigen_taxon_dict.pickle" % config.data_directory 
coarse_eigen_phylo_dict_path = "%scoarse_eigen_phylo_dict.pickle" % config.data_directory 


n_iter = 1000
rarefied=False

def make_coarse_eigen_taxon_dict():

    coarse_eigen_taxon_dict = {}

    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        coarse_eigen_taxon_dict[environment] = {}

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
            max_lambda = lambda_rho[0]

            # get indices for null genera coarse graining
            #counts_genera = [len(n) for n in g_taxa_idx_all]
            #coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_genera))[:-1]

            #idx_all = numpy.arange(len(coarse_grained_idx_all))
            s_by_s_copy = numpy.copy(s_by_s)
            max_all = []
            for i in range(n_iter):

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

            max_all = numpy.sort(max_all)
            p_value = sum(max_all > max_lambda)/n_iter
            lower_ci = max_all[int(n_iter*0.025)]
            upper_ci = max_all[int(n_iter*0.975)]

            coarse_eigen_taxon_dict[environment][rank] = {}
            coarse_eigen_taxon_dict[environment][rank]['max_eigenvalue'] = max_lambda
            coarse_eigen_taxon_dict[environment][rank]['p_value'] = p_value
            coarse_eigen_taxon_dict[environment][rank]['lower_ci'] = lower_ci
            coarse_eigen_taxon_dict[environment][rank]['upper_ci'] = upper_ci

            print(environment, rank, max_lambda, p_value)



    with open(coarse_eigen_taxon_dict_path, 'wb') as handle:
        pickle.dump(coarse_eigen_taxon_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


def make_coarse_eigen_phylo_dict():

    coarse_grained_tree_dict_all = {}
    distances_all = []
    for environment in diversity_utils.environments_to_keep:
        sys.stderr.write("Loading tree dict for %s...\n" % environment)
        coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)
        coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict
        distances = list(coarse_grained_tree_dict.keys())
        distances_all.extend(distances)


    coarse_eigen_phylo_dict = {}
    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]
        pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)

        samples = pres_abs_dict['samples']
        taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))

        sys.stderr.write("Getting site-by-species matrix...\n")
        s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

        distances = list(coarse_grained_tree_dict.keys())
        distances.sort()

        coarse_eigen_phylo_dict[environment] = {}
        sys.stderr.write("Running phylogenetic coarse-graining.....\n")
        for distance_idx, distance in enumerate(distances):

            sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))

            coarse_grained_list = coarse_grained_tree_dict[distance]

            coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            # re-sort the data
            coarse_grained_idx_all_flat = []
            for c in coarse_grained_idx_all:
                coarse_grained_idx_all_flat.extend(c.tolist())

            coarse_grained_idx_all_flat = numpy.asarray(coarse_grained_idx_all_flat)
            n_coarse_grained_all = numpy.asarray([len(c) for c in coarse_grained_idx_all])
            coarse_grain_idx = numpy.append([0], numpy.cumsum(n_coarse_grained_all))[:-1]

            s_by_s_distance = s_by_s[coarse_grained_idx_all_flat,:]
            s_by_s_distance_copy = numpy.copy(s_by_s_distance)

            s_by_s_distance_coarse = numpy.add.reduceat(s_by_s_distance_copy, coarse_grain_idx, axis=0)
            rel_s_by_s_distance_coarse = s_by_s_distance_coarse/numpy.sum(s_by_s_distance_coarse, axis=0)
            print(rel_s_by_s_distance_coarse.shape)
            rho_rel_s_by_s_all_clades = numpy.corrcoef(rel_s_by_s_distance_coarse)
            print(rho_rel_s_by_s_all_clades.shape)
            # eigendecomposition on a real symmetric matrix, returns real eigenvalues
            lambda_rho, v_rho = numpy.linalg.eigh(rho_rel_s_by_s_all_clades)
            lambda_rho = numpy.sort(lambda_rho)[::-1]
            max_lambda = lambda_rho[0]
            print(max_lambda)

            max_all = []
            for i in range(n_iter):

                numpy.random.shuffle(s_by_s_distance_copy)
                s_by_s_coarse_null = numpy.add.reduceat(s_by_s_distance_copy, coarse_grain_idx, axis=0)
                s_by_s_coarse_null = s_by_s_coarse_null[:,~(numpy.all(s_by_s_coarse_null == 0, axis=0))]
                rel_s_by_s_coarse_null = s_by_s_coarse_null/numpy.sum(s_by_s_coarse_null, axis=0)
                rho_rel_s_by_s_all_clades_null = numpy.corrcoef(rel_s_by_s_coarse_null)
                # eigendecomposition on a real symmetric matrix, returns real eigenvalues
                lambda_rho_null, v_rho_null = numpy.linalg.eigh(rho_rel_s_by_s_all_clades_null)
                lambda_rho_null = numpy.sort(lambda_rho_null)[::-1]
                print(lambda_rho_null[0])
                max_all.append(lambda_rho_null[0])


            max_all = numpy.sort(max_all)
            p_value = sum(max_all > max_lambda)/n_iter
            lower_ci = max_all[int(n_iter*0.025)]
            upper_ci = max_all[int(n_iter*0.975)]

            coarse_eigen_phylo_dict[environment][distance] = {}
            coarse_eigen_phylo_dict[environment][distance]['max_eigenvalue'] = max_lambda
            coarse_eigen_phylo_dict[environment][distance]['p_value'] = p_value
            coarse_eigen_phylo_dict[environment][distance]['lower_ci'] = lower_ci
            coarse_eigen_phylo_dict[environment][distance]['upper_ci'] = upper_ci

            print(environment, d, max_lambda, p_value)


    with open(coarse_eigen_phylo_dict_path, 'wb') as handle:
        pickle.dump(coarse_eigen_phylo_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


make_coarse_eigen_phylo_dict()

    
