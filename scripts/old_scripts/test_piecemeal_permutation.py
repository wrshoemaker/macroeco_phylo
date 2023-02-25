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


n_iter = 1000
diversity_vs_diversity_taxon_emp_template = config.data_directory + "test_constrain_occupancy_taxon_%d.pickle"
diversity_vs_diversity_phylo_emp_template = config.data_directory + "test_constrain_occupancy_phylo_%d.pickle"




def make_phylo_dict(block_richnes = 1000):

    diversity_dict = {}
    for environment in diversity_utils.environments_to_keep:

        sys.stderr.write("Loading tree dict for %s...\n" % environment)
        coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=False)

        sys.stderr.write("Subsetting samples for %s...\n" % environment)

        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied=False)
        samples = sad_annotated_dict['samples']

        sys.stderr.write("Building site-by-species matrix...\n")
        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
        richness = len(taxa)

        distances = list(coarse_grained_tree_dict.keys())
        distances.sort()

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
        rel_s_by_s_sort = rel_s_by_s[s_sort_idx,:]
        taxa_sort = taxa[s_sort_idx]
        scaled_rel_s_by_s_sort = (rel_s_by_s_sort.T/mean_abundance_sort).T

        # get coarse_grained
        distances_coarse_grain_dict = {}
        diversity_dict[environment] = {}
        for distance in distances:

            sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
            coarse_grained_list = coarse_grained_tree_dict[distance]
            coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

            # get indexes for each clade
            coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa_sort==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            # re-sort the data
            coarse_grained_idx_all_flat = []
            for c in coarse_grained_idx_all:
                coarse_grained_idx_all_flat.extend(c.tolist())

            coarse_grained_idx_all_flat = numpy.asarray(coarse_grained_idx_all_flat)
            n_coarse_grained_all = numpy.asarray([len(c) for c in coarse_grained_idx_all])
            coarse_grain_idx = numpy.append([0], numpy.cumsum(n_coarse_grained_all))[:-1]

            distances_coarse_grain_dict[distance] = coarse_grain_idx

            # get diversity
            rel_s_by_s_sort_coarse = numpy.add.reduceat(rel_s_by_s_sort, coarse_grain_idx, axis=0)
            diversity = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_s_by_s_sort_coarse)

            diversity_dict[environment][distance] = {}
            diversity_dict[environment][distance]['observed'] = numpy.mean(diversity)
            diversity_dict[environment][distance]['null_constran_mean_and_occupancy_all'] = []
            diversity_dict[environment][distance]['null_constran_mean_all'] = []
            diversity_dict[environment][distance]['null_all'] = []


        for i in range(n_iter):

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
            for distance in distances:

                unscaled_rel_s_by_s_sort_copy_coarse = numpy.add.reduceat(unscaled_rel_s_by_s_sort_copy, distances_coarse_grain_dict[distance], axis=0)
                diversity_null_constrain_mean_and_occupancy = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, unscaled_rel_s_by_s_sort_copy_coarse)
                diversity_dict[environment][distance]['null_constran_mean_and_occupancy_all'].append(numpy.mean(diversity_null_constrain_mean_and_occupancy))

                rel_s_by_s_sort_coarse = numpy.add.reduceat(rel_s_by_s_sort, distances_coarse_grain_dict[distance], axis=0)
                diversity_null_unconstrained = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_s_by_s_sort_coarse)
                diversity_dict[environment][distance]['null_all'].append(numpy.mean(diversity_null_unconstrained))

                unscaled_rel_s_by_s_sort_copy_constrain_mean_coarse = numpy.add.reduceat(unscaled_rel_s_by_s_sort_just_mean_copy, distances_coarse_grain_dict[distance], axis=0)
                diversity_null_constrain_mean = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, unscaled_rel_s_by_s_sort_copy_constrain_mean_coarse)
                diversity_dict[environment][distance]['null_constran_mean_all'].append(numpy.mean(diversity_null_constrain_mean))



        # get confidene intervals
        for distance in distances:

            for measure in ['null_constran_mean_and_occupancy_all', 'null_all', 'null_constran_mean_all']:

                data = numpy.sort(numpy.asarray(diversity_dict[environment][distance][measure]))
                lower_ci_data = data[int((0.025*iter))]
                upper_ci_data = data[int((0.975*iter))]

                del diversity_dict[environment][distance][measure]
                diversity_dict[environment][distance]['%s_lower_ci'%  measure] = lower_ci_data
                diversity_dict[environment][distance]['%s_upper_ci'%  measure] = upper_ci_data

    diversity_vs_diversity_phylo_emp_path = diversity_vs_diversity_phylo_emp_template % block_richnes

    sys.stderr.write("Saving diversity dictionary...\n")
    with open(diversity_vs_diversity_phylo_emp_path, 'wb') as handle:
        pickle.dump(diversity_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


make_phylo_dict()



# plot figure

def plot_figure():

    with open(diversity_vs_diversity_phylo_emp_path, 'rb') as handle:
        diversity_dict = pickle.load(handle)


    fig = plt.figure(figsize = (10, 10)) #
    fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

    environment_chunk_all = [diversity_utils.environments_to_keep[x:x+3] for x in range(0, len(diversity_utils.environments_to_keep), 3)]
    for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

        for environment_idx, environment in enumerate(environment_chunk):

            ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

            distances = list(diversity_dict[environment].keys())
            #distances.remove('diversity_species')
            #distances.sort()

            mean_diveristy = [numpy.mean(diversity_dict[environment][d]['observed']) for d in distances]

            upper_ci = [diversity_dict[environment][d]['null_all_upper_ci'] for d in distances]
            lower_ci = [diversity_dict[environment][d]['null_all_lower_ci'] for d in distances]

            upper_ci_constrain_mad = [diversity_dict[environment][d]['null_constran_mean_all_upper_ci'] for d in distances]
            lower_ci_constrain_mad = [diversity_dict[environment][d]['null_constran_mean_all_lower_ci'] for d in distances]

            upper_ci_constrain_mad_and_occupancy = [diversity_dict[environment][d]['null_constran_mean_and_occupancy_all_upper_ci'] for d in distances]
            lower_ci_constrain_mad_and_occupancy = [diversity_dict[environment][d]['null_constran_mean_and_occupancy_all_lower_ci'] for d in distances]


            ax.plot(distances, mean_diveristy, ls='-', lw=1.5, c='k',  zorder=2)
            ax.scatter(distances, mean_diveristy, c='k', label='Observed',  zorder=3)

            #ax.fill_between(distances, lower_ci, upper_ci, color='#87CEEB', alpha=0.5, label='95% PI', zorder=1)
            ax.fill_between(distances, lower_ci_constrain_mad, upper_ci_constrain_mad, color='#FFA500', alpha=0.5, label='95% PI, constrain ' + r'$\bar{x}$', zorder=1)
            ax.fill_between(distances, lower_ci_constrain_mad_and_occupancy, upper_ci_constrain_mad_and_occupancy, color='#FF6347', alpha=0.5, label='95% PI, constrain ' + r'$\bar{x}$' + ' and ' + r'$\hat{o}$', zorder=1)

            ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)

            ax.set_xscale('log', base=10)

            if environment_idx == 0:
                ax.set_ylabel("Mean Shannon's diversity", fontsize = 10)

            if environment_chunk_idx == len(environment_chunk_all)-1:
                ax.set_xlabel("Phylogenetic distance", fontsize = 12)

            if (environment_chunk_idx ==0) and (environment_idx == 0):
                ax.legend(loc="lower left", fontsize=7)




    #fig.suptitle(label_, fontsize=14, y=0.95)

    fig.subplots_adjust(hspace=0.25,wspace=0.3)
    fig_name = "%sdiversity_slm_constrain_null_emp_phylo.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()
