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
import simulation_utils



rarefied = True
#environment = 'marine metagenome'
measure = 'richness'

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = range(len(taxa_ranks))
taxa_ranks_label = diversity_utils.taxa_ranks_label





for environment in diversity_utils.environments_to_keep:

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    #samples = diversity_utils.subset_observations(environment=environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
    samples = sad_annotated_dict['samples']

    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

    mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
    N = s_by_s.sum(axis=0)

    #div = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s)

    #s_by_s = simulation_utils.genrate_community_from_mean_and_var(mean_rel_s_by_s, var_rel_s_by_s, N, len(N))
    #rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

    #div_sim = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s)

    #print(numpy.mean(div_sim), numpy.mean(div))



    dbd_slope_dict = {}

    for coarse_rank_idx, coarse_rank in enumerate(diversity_utils.taxa_ranks):

        slope_all = []

        if coarse_rank != 'phylum':
            continue

        sys.stderr.write("Starting %s level analysis...\n" % coarse_rank)
        dbd_slope_dict[coarse_rank] = {}

        if coarse_rank == 'genus':
            fine_rank = 'asv'
            all_fine = numpy.copy(taxa).tolist()
        else:
            fine_rank = diversity_utils.taxa_ranks[coarse_rank_idx-1]
            all_fine = list([sad_annotated_dict['taxa'][t][fine_rank] for t in taxa])


        sys.stderr.write("Estimating dbd slope...\n")

        set_fine = list(set(all_fine))
        all_fine = numpy.asarray(all_fine)
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

        fine_rel_s_by_s =  s_by_s_fine/numpy.sum(s_by_s_fine, axis=0)
        coarse_idx = numpy.append([0], numpy.cumsum(counts_coarse))[:-1]
        coarse_rel_s_by_s = numpy.add.reduceat(fine_rel_s_by_s, coarse_idx, axis=0)

        # get index for fine for null
        unique_fine, counts_fine = numpy.unique(all_fine, return_counts=True)
        fine_idx = numpy.append([0], numpy.cumsum(counts_fine))[:-1]
        for focal_coarse_idx, focal_coarse in enumerate(set_coarse):

            # ignore coarse-grained taxa with less than five fine-grained taxa
            if counts_coarse[focal_coarse_idx] < 1:
                continue

            # all the fine-scale indices for the focal
            focal_coarse_s_by_s_idx = numpy.asarray(idx_to_sort[focal_coarse_idx])
            focal_fine_rel_s_by_s = fine_rel_s_by_s[focal_coarse_s_by_s_idx,:]
            nonfocal_coarse_rel_s_by_s = numpy.delete(coarse_rel_s_by_s, focal_coarse_idx, axis=0)

            n_fine_focal_coarse = len(focal_coarse_s_by_s_idx)

            if measure == 'richness':
                measure_coarse = numpy.sum(nonfocal_coarse_rel_s_by_s > 0, axis=0)
                measure_fine = numpy.sum(focal_fine_rel_s_by_s > 0, axis=0)

            else:
                measure_coarse = numpy.apply_along_axis(diversity_utils.calculate_richness, 0, nonfocal_coarse_rel_s_by_s)
                measure_fine = numpy.apply_along_axis(diversity_utils.calculate_richness, 0, focal_fine_rel_s_by_s)

            # remove sites where there are no observations in either focal or non-focal
            idx_to_remove = (measure_fine>0) | (measure_coarse>0)

            measure_coarse = measure_coarse[idx_to_remove]
            measure_fine = measure_fine[idx_to_remove]

            if len(measure_fine) < 10:
                continue

            slope, intercept, r_value, p_value, std_err = stats.linregress(measure_coarse, measure_fine)
            slope_all.append(slope)


        print(environment, len(coarse_idx), numpy.mean(slope_all))







def old__():


    for environment in diversity_utils.environments_to_keep:

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #samples = diversity_utils.subset_observations(environment=environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

        mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
        var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
        N = s_by_s.sum(axis=0)


        s_by_s_sim = simulation_utils.genrate_community_from_mean_and_var(mean_rel_s_by_s, var_rel_s_by_s, N, len(N))
        #rel_s_by_s_sim = (s_by_s_sim/s_by_s_sim.sum(axis=0))


        #sys.stderr.write("Running taxonomic coarse-graining.....\n")
        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)

            all_genera = numpy.asarray(list([sad_annotated_dict['taxa'][t][rank] for t in taxa]))
            #all_genera_idx = numpy.arange(len(all_genera))

            unique_genera, counts_genera = numpy.unique(all_genera, return_counts=True)
            coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_genera))[:-1]

            s_by_s_coarse = numpy.add.reduceat(s_by_s, coarse_grain_by_genus_idx, axis=0)

            s_by_s_sim_coarse = numpy.add.reduceat(s_by_s_sim, coarse_grain_by_genus_idx, axis=0)


            slope_all = []
            idx_all = numpy.arange(len(counts_genera))
            for i in idx_all:

                if counts_genera[i] < 10:
                    continue

                if i == idx_all[-1]:
                    s_by_s_focal = s_by_s[coarse_grain_by_genus_idx[i]:,:]

                else:
                    s_by_s_focal = s_by_s[coarse_grain_by_genus_idx[i]:coarse_grain_by_genus_idx[i+1],:]

                #s_by_s_focal = s_by_s[coarse_grained_idx_all[i],:]
                s_by_s_non_focal = s_by_s_coarse[idx_all!=i,:]

                diversity_focal = numpy.apply_along_axis(diversity_utils.calculate_richness, 0, s_by_s_focal)
                diversity_non_focal = numpy.apply_along_axis(diversity_utils.calculate_richness, 0, s_by_s_non_focal)

                idx_to_keep = (diversity_focal>0) & (diversity_non_focal>0) & (numpy.isfinite(diversity_focal)) & (numpy.isfinite(diversity_non_focal))
                diversity_focal = diversity_focal[idx_to_keep]
                diversity_non_focal = diversity_non_focal[idx_to_keep]

                # ignore cut-offs with low richness
                if len(diversity_focal) < 10:
                    continue

                ## ignore cut-offs where all diversity estimates are the same
                #if (sum(diversity_non_focal==1) == len(diversity_non_focal)) or (sum(diversity_focal==1) == len(diversity_focal)):
                #    continue

                slope, intercept, r_value, p_value, std_err = stats.linregress(diversity_non_focal, diversity_focal)

                slope_all.append(slope)

            print(environment, rank, len(unique_genera), numpy.mean(slope))



        #    genus_to_taxa_dict = {}
        #    sad_genera_all = []
        #    for genus_idx, genus in enumerate(all_genera):
        #        genus_to_taxa_dict[genus] = []
        #        for t in taxa:
        #            if sad_annotated_dict['taxa'][t][rank] == genus:
        #                genus_to_taxa_dict[genus].append(t)

        #        g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])

        #        rel_s_by_s_genus = rel_s_by_s[g_taxa_idx,:]
