import os
import random
import copy
import math
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



#environment = 'human gut metagenome'

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label

diversity_vs_diversity_taxon_emp_template = config.data_directory + "test_constrain_occupancy_taxon_%d.pickle"



block_richnes = 100

def make_taxon_dict(iter = 1000):

    taxa_ranks = diversity_utils.taxa_ranks
    idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
    taxa_ranks_label = diversity_utils.taxa_ranks_label

    diversity_dict = {}
    for environment in diversity_utils.environments_to_keep:

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied=False)
        samples = sad_annotated_dict['samples']

        sys.stderr.write("Building site-by-species matrix...\n")
        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
        richness = len(taxa)

        # figure out number of blocks
        n_blocks = math.ceil(richness / block_richnes)
        block_richness_all = [block_richnes]*(n_blocks-1)
        # add remaining species
        block_richness_all.append(richness - block_richnes*(n_blocks-1))
        block_richness_all = numpy.asarray(block_richness_all)

        rel_s_by_s = s_by_s/numpy.sum(s_by_s,axis=0)
        mean_abundance = numpy.mean(rel_s_by_s, axis=1)

        # create array of block numbers
        # sort rel_s_by_s by mean abundance
        s_idx = list(range(s_by_s.shape[0]))
        mean_abundance_s_idx_tuple_all = list(zip(mean_abundance.tolist(), s_idx))

        # sort tuples by mean abundance, order of increasing abundane
        mean_abundance_s_idx_tuple_all.sort(key=lambda tup: tup[0])
        s_sort_idx = numpy.asarray([s[1] for s in mean_abundance_s_idx_tuple_all])
        mean_abundance_sort = numpy.asarray([s[0] for s in mean_abundance_s_idx_tuple_all])

        s_by_s_sort = s_by_s[s_sort_idx,:]
        rel_s_by_s_sort = rel_s_by_s[s_sort_idx,:]
        taxa_sort = taxa[s_sort_idx]
        scaled_rel_s_by_s_sort = (rel_s_by_s_sort.T/mean_abundance_sort).T

        # get coarse_grained
        rank_coarse_grain_dict = {}
        diversity_dict[environment] = {}
        for rank_idx, rank in enumerate(taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)
            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa_sort])))
            all_genera_idx = numpy.arange(len(all_genera))

            genus_to_taxa_dict = {}
            sad_genera_all = []
            #coarse_grained_n = []
            #coarse_grained_list = []
            coarse_grained_idx_all = []
            for genus_idx, genus in enumerate(all_genera):
                genus_to_taxa_dict[genus] = []
                for t in taxa_sort:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa_sort == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
                s_by_s_sort_genus = s_by_s_sort[g_taxa_idx,:]
                #coarse_grained_n.append(len(g_taxa_idx))
                sad_genera_all.append(s_by_s_sort_genus.sum(axis=0))
                coarse_grained_idx_all.append(g_taxa_idx)


            #coarse_grained_n = numpy.asarray(coarse_grained_n)
            # get indexes for each clade
            #coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa_sort==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            # re-sort the data
            coarse_grained_idx_all_flat = []
            for c in coarse_grained_idx_all:
                coarse_grained_idx_all_flat.extend(c.tolist())

            coarse_grained_idx_all_flat = numpy.asarray(coarse_grained_idx_all_flat)
            n_coarse_grained_all = numpy.asarray([len(c) for c in coarse_grained_idx_all])
            coarse_grain_idx = numpy.append([0], numpy.cumsum(n_coarse_grained_all))[:-1]

            rank_coarse_grain_dict[rank] = coarse_grain_idx

            # get diversity
            rel_s_by_s_sort_coarse = numpy.add.reduceat(rel_s_by_s_sort, coarse_grain_idx, axis=0)
            diversity = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_s_by_s_sort_coarse)

            diversity_dict[environment][rank] = {}
            diversity_dict[environment][rank]['observed'] = numpy.mean(diversity)
            diversity_dict[environment][rank]['null_constran_mean_and_occupancy_all'] = []
            diversity_dict[environment][rank]['null_constran_mean_all'] = []
            diversity_dict[environment][rank]['null_all'] = []


        for i in range(iter):

            if (i%100 == 0) and (i>0):
                print(i)

            # shuffle unconstrained
            numpy.random.shuffle(rel_s_by_s_sort)

            # shuffle constrained on mean and occupancy
            scaled_rel_s_by_s_sort_copy = numpy.copy(scaled_rel_s_by_s_sort)

            numpy.random.shuffle(block_richness_all)
            permutation_block_idx = numpy.append([0], numpy.cumsum(block_richness_all))[:-1]
            # add final index
            permutation_block_idx = numpy.append(permutation_block_idx, richness)

            # shuffle each block
            for n in range(n_blocks):
                numpy.random.shuffle(scaled_rel_s_by_s_sort_copy[permutation_block_idx[n]:permutation_block_idx[n+1],:])

            unscaled_rel_s_by_s_sort_copy = ((scaled_rel_s_by_s_sort_copy.T)*mean_abundance_sort).T
            # renormalize in the off-change rows sum to more than one
            unscaled_rel_s_by_s_sort_copy = unscaled_rel_s_by_s_sort_copy/numpy.sum(unscaled_rel_s_by_s_sort_copy,axis=0)

            # shuffle constrained on just the mean
            scaled_rel_s_by_s_sort_just_mean_copy = numpy.copy(scaled_rel_s_by_s_sort)
            numpy.random.shuffle(scaled_rel_s_by_s_sort_just_mean_copy)
            unscaled_rel_s_by_s_sort_just_mean_copy = ((scaled_rel_s_by_s_sort_just_mean_copy.T)*mean_abundance_sort).T


            # coarse-grain
            for rank in taxa_ranks:

                unscaled_rel_s_by_s_sort_copy_coarse = numpy.add.reduceat(unscaled_rel_s_by_s_sort_copy, rank_coarse_grain_dict[rank], axis=0)
                diversity_null_constrain_mean_and_occupancy = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, unscaled_rel_s_by_s_sort_copy_coarse)
                diversity_dict[environment][rank]['null_constran_mean_and_occupancy_all'].append(numpy.mean(diversity_null_constrain_mean_and_occupancy))

                rel_s_by_s_sort_coarse = numpy.add.reduceat(rel_s_by_s_sort, rank_coarse_grain_dict[rank], axis=0)
                diversity_null_unconstrained = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_s_by_s_sort_coarse)
                diversity_dict[environment][rank]['null_all'].append(numpy.mean(diversity_null_unconstrained))

                unscaled_rel_s_by_s_sort_copy_constrain_mean_coarse = numpy.add.reduceat(unscaled_rel_s_by_s_sort_just_mean_copy, rank_coarse_grain_dict[rank], axis=0)
                diversity_null_constrain_mean = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, unscaled_rel_s_by_s_sort_copy_constrain_mean_coarse)
                diversity_dict[environment][rank]['null_constran_mean_all'].append(numpy.mean(diversity_null_constrain_mean))


        # get confidence intervals
        for rank in taxa_ranks:

            for measure in ['null_constran_mean_and_occupancy_all', 'null_all', 'null_constran_mean_all']:

                data = numpy.sort(numpy.asarray(diversity_dict[environment][rank][measure]))
                lower_ci_data = data[int((0.025*iter))]
                upper_ci_data = data[int((0.975*iter))]

                del diversity_dict[environment][rank][measure]
                diversity_dict[environment][rank]['%s_lower_ci'%  measure] = lower_ci_data
                diversity_dict[environment][rank]['%s_upper_ci'%  measure] = upper_ci_data


    diversity_vs_diversity_taxon_emp_path = diversity_vs_diversity_taxon_emp_template % block_richnes

    sys.stderr.write("Saving diversity dictionary...\n")
    with open(diversity_vs_diversity_taxon_emp_path, 'wb') as handle:
        pickle.dump(diversity_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



#make_taxon_dict()
with open(diversity_vs_diversity_taxon_emp_template % block_richnes, 'rb') as handle:
    diversity_dict = pickle.load(handle)


fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environment_chunk_all = [diversity_utils.environments_to_keep[x:x+3] for x in range(0, len(diversity_utils.environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        distances = list(diversity_dict[environment].keys())

        rank_idx = list(range(len(distances)))

        mean_diveristy = [numpy.mean(diversity_dict[environment][d]['observed']) for d in distances]

        upper_ci = [diversity_dict[environment][d]['null_all_upper_ci'] for d in distances]
        lower_ci = [diversity_dict[environment][d]['null_all_lower_ci'] for d in distances]

        upper_ci_constrain_mad = [diversity_dict[environment][d]['null_constran_mean_all_upper_ci'] for d in distances]
        lower_ci_constrain_mad = [diversity_dict[environment][d]['null_constran_mean_all_lower_ci'] for d in distances]

        upper_ci_constrain_mad_and_occupancy = [diversity_dict[environment][d]['null_constran_mean_and_occupancy_all_upper_ci'] for d in distances]
        lower_ci_constrain_mad_and_occupancy = [diversity_dict[environment][d]['null_constran_mean_and_occupancy_all_lower_ci'] for d in distances]


        ax.plot(rank_idx, mean_diveristy, ls='-', lw=1.5, c='k',  zorder=2)
        ax.scatter(rank_idx, mean_diveristy, c='k', label='Observed',  zorder=3)

        #ax.fill_between(distances, lower_ci, upper_ci, color='#87CEEB', alpha=0.5, label='95% PI', zorder=1)
        ax.fill_between(rank_idx, lower_ci_constrain_mad, upper_ci_constrain_mad, color='#FFA500', alpha=0.5, label='95% PI, constrain ' + r'$\bar{x}$', zorder=1)
        ax.fill_between(rank_idx, lower_ci_constrain_mad_and_occupancy, upper_ci_constrain_mad_and_occupancy, color='#FF6347', alpha=0.5, label='95% PI, constrain ' + r'$\bar{x}$' + ' and ' + r'$\hat{o}$', zorder=1)

        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)

        ax.set_xticks(idx_taxa_ranks)
        ax.set_xticklabels(taxa_ranks_label, fontsize=9.5)

        if environment_idx == 0:
            ax.set_ylabel("Mean Shannon's diversity", fontsize = 10)

        #if environment_chunk_idx == len(environment_chunk_all)-1:
        #    ax.set_xlabel("Phylogenetic distance", fontsize = 12)

        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="lower left", fontsize=7)



fig.subplots_adjust(hspace=0.25,wspace=0.3)
fig_name = "%sdiversity_slm_constrain_null_emp_taxon.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
