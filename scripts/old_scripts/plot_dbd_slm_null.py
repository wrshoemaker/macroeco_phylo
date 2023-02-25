from __future__ import division
import config
import os
import sys
import subprocess
import random
import re
import itertools
import pickle
from collections import Counter

import numpy
import diversity_utils
import scipy.stats as stats

import ete3
import tree_utils
import random
import dbd_utils
import plot_utils



#coarse_grained_tree_dict = dbd_utils.get_coarse_grained_tree_dict()
taxa_ranks = diversity_utils.taxa_ranks


richness_dbd_null_dict_path = config.data_directory + "richness_dbd_null_dict.pickle"
environments_to_keep = diversity_utils.environments_to_keep
environments_to_keep = [environments_to_keep[0]]


dbd_dict = {}
coarse_idx_dict = {}
coarse_idx_dict['taxon'] = {}
coarse_idx_dict['phylo'] = {}


def make_null_make_dbd_slm_dict_old():

    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        sys.stderr.write("Loading taxon dict for %s...\n" % environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied=False)
        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s = s_by_s/s_by_s.sum(axis=0)

        coarse_idx_dict[environment] = {}
        coarse_idx_dict[environment]['taxon'] = {}
        coarse_idx_dict[environment]['phylo'] = {}


        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)
            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))

            genus_to_taxa_dict = {}
            coarse_grained_idx_all = []
            for genus_idx, genus in enumerate(all_genera):
                genus_to_taxa_dict[genus] = []
                for t in taxa:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
                coarse_grained_idx_all.append(g_taxa_idx)


            #coarse_grained_n = numpy.asarray(coarse_grained_n)
            # get indexes for each clade
            # re-sort the data
            coarse_grained_idx_all_flat = []
            for c in coarse_grained_idx_all:
                coarse_grained_idx_all_flat.extend(c.tolist())

            coarse_grained_idx_all_flat = numpy.asarray(coarse_grained_idx_all_flat)
            n_coarse_grained_all = numpy.asarray([len(c) for c in coarse_grained_idx_all])
            coarse_grain_idx = numpy.append([0], numpy.cumsum(n_coarse_grained_all))[:-1]

            coarse_idx_dict[environment]['taxon'][rank] = coarse_grain_idx


        rank_ratio_idx_dict = {}
        # create map for each pair of ranks to identify the number 


        s_by_s_copy = numpy.copy(s_by_s)
        for i in range(1):

            numpy.random.shuffle(s_by_s_copy)

            for rank_idx in range(len(taxa_ranks)-1):

                idx_fine = coarse_idx_dict[environment]['taxon'][taxa_ranks[rank_idx]]
                idx_coarse = coarse_idx_dict[environment]['taxon'][taxa_ranks[rank_idx+1]]

                print(idx_fine, idx_coarse)

                s_by_s_fine = numpy.add.reduceat(s_by_s_copy, idx_fine, axis=0)
                s_by_s_coarse = numpy.add.reduceat(s_by_s_copy, idx_coarse, axis=0)



                #print(s_by_s_fine.shape, s_by_s_coarse.shape)

            #s_by_s_non_focal_coarse_null = numpy.add.reduceat(s_by_s_non_focal_null, coarse_grain_by_genus_idx, axis=0)




def make_null_make_dbd_slm_dict(iter=1000):

    dbd_dict = {}

    for environment_idx, environment in enumerate(environments_to_keep):

        sys.stderr.write("Loading taxon dict for %s...\n" % environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied=False)
        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
        s_by_s = numpy.copy(s_by_s)

        dbd_dict[environment] = {}
        dbd_dict[environment]['taxon'] = {}
        dbd_dict[environment]['phylo'] = {}

        for coarse_rank_idx, coarse_rank in enumerate(diversity_utils.taxa_ranks):
            dbd_dict[environment]['taxon'][coarse_rank] = {}
            dbd_dict[environment]['taxon'][coarse_rank]['mean_slope'] = []
            dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm'] = []

        for i in range(iter):

            numpy.random.shuffle(s_by_s)

            for coarse_rank_idx, coarse_rank in enumerate(diversity_utils.taxa_ranks):

                #sys.stderr.write("Starting %s level analysis...\n" % coarse_rank)

                if coarse_rank == 'genus':
                    fine_rank = 'asv'
                    all_fine = numpy.copy(taxa).tolist()

                #elif coarse_rank == 'phylum':
                #    continue

                else:
                    fine_rank = diversity_utils.taxa_ranks[coarse_rank_idx-1]
                    all_fine = list([sad_annotated_dict['taxa'][t][fine_rank] for t in taxa])

                sys.stderr.write("Estimating dbd slope...\n")
                set_fine = list(set(all_fine))
                set_fine =  numpy.asarray(set_fine)

                fine_coarse_dict = {}
                for t in taxa:

                    if fine_rank == 'asv':
                        fine_coarse_dict[t] = sad_annotated_dict['taxa'][t][coarse_rank]

                    else:
                        fine_rank_t = sad_annotated_dict['taxa'][t][fine_rank]
                        if fine_rank_t not in fine_coarse_dict:
                            fine_coarse_dict[fine_rank_t] = sad_annotated_dict['taxa'][t][coarse_rank]


                if fine_rank == 'asv':
                    s_by_s_fine = numpy.copy(s_by_s)
                    all_coarse = []
                    for fine in set_fine:
                        all_coarse.append(fine_coarse_dict[fine])
                    all_coarse = numpy.asarray(all_coarse)


                else:
                    sad_fine_all = []
                    all_coarse = []
                    for fine in set_fine:
                        taxa_in_fine = []
                        for t in taxa:

                            if sad_annotated_dict['taxa'][t][fine_rank] == fine:
                                taxa_in_fine.append(t)

                        # numpy.where() is a major bottleneck
                        g_taxa_idx = numpy.asarray([numpy.where(taxa == taxa_in_fine_i)[0][0] for taxa_in_fine_i in taxa_in_fine])
                        sad_fine_all.append(s_by_s[g_taxa_idx,:].sum(axis=0))
                        all_coarse.append(fine_coarse_dict[fine])

                    all_coarse = numpy.asarray(all_coarse)
                    s_by_s_fine = numpy.stack(sad_fine_all, axis=0)


                # remove sites where there are no observations
                s_by_s_fine = s_by_s_fine[:,~(numpy.all(s_by_s_fine == 0, axis=0))]

                # sort s_by_s_fine
                # so we can use cumsum
                set_coarse = list(set(all_coarse))
                idx_to_sort = []
                counts_coarse = []
                for coarse in set_coarse:
                    idx_coarse_i = numpy.where(all_coarse==coarse)[0].tolist()
                    counts_coarse.append(len(idx_coarse_i))
                    idx_to_sort.append(idx_coarse_i)

                idx_to_sort_flat =  list(itertools.chain(*idx_to_sort))
                idx_to_sort_flat = numpy.asarray(idx_to_sort_flat)
                s_by_s_fine = s_by_s_fine[idx_to_sort_flat,:]

                total_reads_fine = numpy.sum(s_by_s_fine, axis=0)
                rel_s_by_s_fine =  s_by_s_fine/numpy.sum(s_by_s_fine, axis=0)
                mean_rel_s_by_s_fine = numpy.mean(rel_s_by_s_fine, axis=1)
                var_rel_s_by_s_fine = numpy.var(rel_s_by_s_fine, axis=1)

                coarse_idx = numpy.append([0], numpy.cumsum(counts_coarse))[:-1]
                coarse_s_by_s = numpy.add.reduceat(s_by_s_fine, coarse_idx, axis=0)       

                total_reads_coarse = numpy.sum(coarse_s_by_s, axis=0)
                rel_coarse_s_by_s = coarse_s_by_s/numpy.sum(coarse_s_by_s, axis=0)
                mean_rel_s_by_s_coarse = numpy.mean(rel_coarse_s_by_s, axis=1)
                var_rel_s_by_s_coarse = numpy.var(rel_coarse_s_by_s, axis=1)

                slope_all = []
                slope_slm_all = []

                for focal_coarse_idx, focal_coarse in enumerate(set_coarse):

                    # ignore coarse-grained taxa with less than five fine-grained taxa
                    if counts_coarse[focal_coarse_idx] < 5:
                        continue

                    # all the fine-scale indices for the focal
                    focal_coarse_s_by_s_idx = numpy.asarray(idx_to_sort[focal_coarse_idx])
                    focal_fine_s_by_s = s_by_s_fine[focal_coarse_s_by_s_idx,:]

                    focal_mean_fine_rel_s_by_s = mean_rel_s_by_s_fine[focal_coarse_s_by_s_idx]
                    focal_var_fine_rel_s_by_s = var_rel_s_by_s_fine[focal_coarse_s_by_s_idx]

                    nonfocal_coarse_s_by_s = numpy.delete(coarse_s_by_s, focal_coarse_idx, axis=0)
                    
                    nonfocal_mean_coarse_rel_s_by_s = numpy.delete(mean_rel_s_by_s_coarse, focal_coarse_idx, axis=0)          
                    nonfocal_var_coarse_rel_s_by_s = numpy.delete(var_rel_s_by_s_coarse, focal_coarse_idx, axis=0)      

                    focal_fine_richness_predicted, nonfocal_coarse_richness_predicted = dbd_utils.predict_richness_dbd(focal_mean_fine_rel_s_by_s, focal_var_fine_rel_s_by_s, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, nonfocal_var_coarse_rel_s_by_s, total_reads_coarse)    

                    measure_coarse = numpy.sum(nonfocal_coarse_s_by_s > 0, axis=0)
                    measure_fine = numpy.sum(focal_fine_s_by_s > 0, axis=0)

                    # remove sites where there are no observations in either focal or non-focal
                    idx_to_remove = (measure_fine>0) | (measure_coarse>0)

                    measure_coarse = measure_coarse[idx_to_remove]
                    measure_fine = measure_fine[idx_to_remove]

                    nonfocal_coarse_richness_predicted = nonfocal_coarse_richness_predicted[idx_to_remove]
                    focal_fine_richness_predicted = focal_fine_richness_predicted[idx_to_remove]

                    if len(measure_fine) < 5:
                        continue

                    slope, intercept, r_value, p_value, std_err = stats.linregress(measure_coarse, measure_fine)
                    slope_slm, intercept_slm, r_value_slm, p_value_slm, std_err_slm = stats.linregress(nonfocal_coarse_richness_predicted, focal_fine_richness_predicted)

                    slope_all.append(slope)
                    slope_slm_all.append(slope_slm)

                slope_all = numpy.asarray(slope_all)
                slope_slm_all = numpy.asarray(slope_slm_all)

                idx_to_keep = ~(numpy.isnan(slope_all) | numpy.isnan(slope_slm_all))

                if sum(idx_to_keep) < 3:
                    continue

                slope_all = slope_all[idx_to_keep]
                slope_slm_all = slope_slm_all[idx_to_keep]

                mean_slope_all = numpy.mean(slope_all)
                mean_slope_slm_all = numpy.mean(slope_slm_all)

                dbd_dict[environment]['taxon'][coarse_rank]['mean_slope'].append(mean_slope_all)
                dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm'].append(mean_slope_slm_all)



    dbd_dict_final = {}
    # make dict with 95% CIs
    for environment_idx, environment in enumerate(environments_to_keep):

        dbd_dict_final[environment] = {}
        dbd_dict_final[environment]['taxon'] = {}
        for coarse_rank_idx, coarse_rank in enumerate(diversity_utils.taxa_ranks):
            dbd_dict_final[environment]['taxon'][coarse_rank] = {}

            mean_slope = numpy.asarray(dbd_dict[environment]['taxon'][coarse_rank]['mean_slope'])
            mean_slope_slm = numpy.asarray(dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm'])

            mean_slope = numpy.sort(mean_slope)
            mean_slope_slm = numpy.sort(mean_slope_slm)

            dbd_dict_final[environment]['taxon'][coarse_rank]['mean_mean_slope'] = numpy.mean(mean_slope)
            dbd_dict_final[environment]['taxon'][coarse_rank]['mean_mean_slope_slm'] = numpy.mean(mean_slope_slm)

            dbd_dict_final[environment]['taxon'][coarse_rank]['lower_ci_mean_slope'] = mean_slope[int(0.025*iter)]
            dbd_dict_final[environment]['taxon'][coarse_rank]['upper_ci_mean_slope'] = mean_slope[int(0.975*iter)]

            dbd_dict_final[environment]['taxon'][coarse_rank]['lower_ci_mean_slope_slm'] = mean_slope_slm[int(0.025*iter)]
            dbd_dict_final[environment]['taxon'][coarse_rank]['upper_ci_mean_slope_slm'] = mean_slope_slm[int(0.975*iter)]




    sys.stderr.write("Saving diversity dictionary...\n")
    with open(richness_dbd_null_dict_path, 'wb') as handle:
        pickle.dump(dbd_dict_final, handle, protocol=pickle.HIGHEST_PROTOCOL)
    


make_null_make_dbd_slm_dict()



