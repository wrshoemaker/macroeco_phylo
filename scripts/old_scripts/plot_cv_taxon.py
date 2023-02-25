
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
import matplotlib.gridspec as gridspec

import dbd_utils
import plot_utils

import diversity_utils
import config

from scipy.special import rel_entr



rarefied = False

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

taxa_ranks = diversity_utils.taxa_ranks_with_asv
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label_with_asv
#taxa_ranks.insert(0, 'OTU')



mean_cv_squared_dict_path = config.data_directory + "mean_cv_squared_dict.pickle"









def make_mean_cv_squared_dict(iter=1000):

    coarse_grained_tree_dict_all = {}
    distances_all = []
    for environment in environments_to_keep:
        #dbd_utils.make_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)
        sys.stderr.write("Loading tree dict for %s...\n" % environment)
        coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)
        coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict
        distances = list(coarse_grained_tree_dict.keys())
        distances_all.extend(distances)

    distances_all = list(set(distances_all))
    distances_all.sort()

    mean_cv_square_dict = {}

    for environment_idx, environment in enumerate(environments_to_keep):

        sys.stderr.write("Running taxonomic coarse-graining.....\n")

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s = s_by_s/s_by_s.sum(axis=0)
        rel_s_by_s_copy = numpy.copy(rel_s_by_s)

        mean_cv_square_dict[environment] = {}
        mean_cv_square_dict[environment]['taxon'] = {}
        mean_cv_square_dict[environment]['phylo'] = {}

        for rank_idx, rank in enumerate(taxa_ranks):

            if rank == 'OTU':

                s_by_s_genera = s_by_s

            else:

                sys.stderr.write("Starting %s level analysis...\n" % rank)
                all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
                all_genera_idx = numpy.arange(len(all_genera))

                genus_to_taxa_dict = {}
                sad_genera_all = []
                coarse_grained_idx_all = []
                for genus_idx, genus in enumerate(all_genera):
                    genus_to_taxa_dict[genus] = []
                    for t in taxa:
                        if sad_annotated_dict['taxa'][t][rank] == genus:
                            genus_to_taxa_dict[genus].append(t)

                    g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
                    s_by_s_genus = s_by_s[g_taxa_idx,:]
                    sad_genera_all.append(s_by_s_genus.sum(axis=0))
                    coarse_grained_idx_all.append(g_taxa_idx)


                s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
                # remove sites where there are no observations
                s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]


            rel_s_by_s_genera = s_by_s_genera/s_by_s_genera.sum(axis=0)
            mean_rel_abundance = numpy.mean(rel_s_by_s_genera, axis=1)
            var_rel_abundance = numpy.var(rel_s_by_s_genera, axis=1)

            cv_rel_abundance = numpy.sqrt(var_rel_abundance)/mean_rel_abundance

            square_cv_rel_abundance = cv_rel_abundance**2
            mean_square_cv_rel_abundance = numpy.mean(square_cv_rel_abundance)

            mean_cv_square_dict[environment]['taxon'][rank] = {}
            mean_cv_square_dict[environment]['taxon'][rank]['mean_observed'] = mean_square_cv_rel_abundance

            # null distribution
            if rank == 'OTU':
                
                mean_cv_square_dict[environment]['taxon'][rank]['mean_null'] = mean_square_cv_rel_abundance
                mean_cv_square_dict[environment]['taxon'][rank]['lower_ci_null'] = mean_square_cv_rel_abundance
                mean_cv_square_dict[environment]['taxon'][rank]['upper_ci_null'] = mean_square_cv_rel_abundance

            else:

                coarse_grained_idx_all_flat = []
                for c in coarse_grained_idx_all:
                    coarse_grained_idx_all_flat.extend(c.tolist())

                coarse_grained_idx_all_flat = numpy.asarray(coarse_grained_idx_all_flat)
                n_coarse_grained_all = numpy.asarray([len(c) for c in coarse_grained_idx_all])
                coarse_grain_idx = numpy.append([0], numpy.cumsum(n_coarse_grained_all))[:-1]

                mean_cv_square_dict[environment]['taxon'][rank]['null_all'] = []
                for i in range(iter):

                    numpy.random.shuffle(rel_s_by_s_copy)
                    rel_s_by_s_coarse_null = numpy.add.reduceat(rel_s_by_s_copy, coarse_grain_idx, axis=0)

                    mean_rel_abundance_null = numpy.mean(rel_s_by_s_coarse_null, axis=1)
                    var_rel_abundance_null = numpy.var(rel_s_by_s_coarse_null, axis=1)

                    cv_rel_abundance_null = numpy.sqrt(var_rel_abundance_null)/mean_rel_abundance_null

                    square_cv_rel_abundance_null = cv_rel_abundance_null**2
                    mean_square_cv_rel_abundance_null = numpy.mean(square_cv_rel_abundance_null)
                    mean_cv_square_dict[environment]['taxon'][rank]['null_all'].append(mean_square_cv_rel_abundance_null)


                null_all = numpy.asarray(mean_cv_square_dict[environment]['taxon'][rank]['null_all'])
                null_all = numpy.sort(null_all)

                mean_cv_square_dict[environment]['taxon'][rank]['mean_null'] = numpy.mean(null_all)
                mean_cv_square_dict[environment]['taxon'][rank]['lower_ci_null'] = null_all[int(iter*0.025)]
                mean_cv_square_dict[environment]['taxon'][rank]['upper_ci_null'] = null_all[int(iter*0.975)]
                del mean_cv_square_dict[environment]['taxon'][rank]['null_all']

                sys.stderr.write("%s mean CV^2 = %f...\n" % (rank, round(mean_square_cv_rel_abundance, 3)))


        sys.stderr.write("Running phylogenetic coarse-graining.....\n")

        coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #samples = diversity_utils.subset_observations(environment=environment)
        pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)

        samples = pres_abs_dict['samples']
        taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
        sys.stderr.write("Getting site-by-species matrix...\n")
        #s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)

        s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s = s_by_s/s_by_s.sum(axis=0)
        rel_s_by_s_copy = numpy.copy(rel_s_by_s)

        distances = coarse_grained_tree_dict.keys()

        for distance_idx, distance in enumerate(distances):

            coarse_grained_list = coarse_grained_tree_dict[distance]
            coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

            # get indexes for each clade
            coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

            # coarse grain s-by-s for all clades
            s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)

            rel_s_by_s_all_clades = s_by_s_all_clades/s_by_s_all_clades.sum(axis=0)
            
            mean_rel_abundance = numpy.mean(rel_s_by_s_all_clades, axis=1)
            var_rel_abundance = numpy.var(rel_s_by_s_all_clades, axis=1)

            cv_rel_abundance = numpy.sqrt(var_rel_abundance)/mean_rel_abundance

            square_cv_rel_abundance = cv_rel_abundance**2
            mean_square_cv_rel_abundance = numpy.mean(square_cv_rel_abundance)

            
            mean_cv_square_dict[environment]['phylo'][distance] = {}
            mean_cv_square_dict[environment]['phylo'][distance]['mean_observed'] = mean_square_cv_rel_abundance
            mean_cv_square_dict[environment]['phylo'][distance]['null_all'] = []
            coarse_grain_idx = numpy.append([0], numpy.cumsum(coarse_grained_n))[:-1]

            for i in range(iter):

                numpy.random.shuffle(rel_s_by_s_copy)
                rel_s_by_s_coarse_null = numpy.add.reduceat(rel_s_by_s_copy, coarse_grain_idx, axis=0)

                mean_rel_abundance_null = numpy.mean(rel_s_by_s_coarse_null, axis=1)
                var_rel_abundance_null = numpy.var(rel_s_by_s_coarse_null, axis=1)

                cv_rel_abundance_null = numpy.sqrt(var_rel_abundance_null)/mean_rel_abundance_null

                square_cv_rel_abundance_null = cv_rel_abundance_null**2
                mean_square_cv_rel_abundance_null = numpy.mean(square_cv_rel_abundance_null)

                mean_cv_square_dict[environment]['phylo'][distance]['null_all'].append(mean_square_cv_rel_abundance_null)


            null_all = numpy.asarray(mean_cv_square_dict[environment]['phylo'][distance]['null_all'])
            null_all = numpy.sort(null_all)

            mean_cv_square_dict[environment]['phylo'][distance]['mean_null'] = numpy.mean(null_all)
            mean_cv_square_dict[environment]['phylo'][distance]['lower_ci_null'] = null_all[int(iter*0.025)]
            mean_cv_square_dict[environment]['phylo'][distance]['upper_ci_null'] = null_all[int(iter*0.975)]
            del mean_cv_square_dict[environment]['phylo'][distance]['null_all']

            sys.stderr.write("Phylo distance = %f mean CV^2 = %f...\n" % (round(distance, 3), round(mean_square_cv_rel_abundance, 3)))


    sys.stderr.write("Saving diversity dictionary...\n")
    with open(mean_cv_squared_dict_path, 'wb') as handle:
        pickle.dump(mean_cv_square_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)





def load_mean_cv_squared_dict():

    with open(mean_cv_squared_dict_path, 'rb') as handle:
        distance_collapsed_dict = pickle.load(handle)
    return distance_collapsed_dict








#make_mean_cv_squared_dict()
mean_cv_squared_dict = load_mean_cv_squared_dict()



fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

idx_taxa_ranks = numpy.asarray(list(range(len(diversity_utils.taxa_ranks_label_with_asv))))
taxa_ranks_label = diversity_utils.taxa_ranks_label_with_asv

color_all = numpy.asarray([plot_utils.rgb_blue_taxon(r) for r in idx_taxa_ranks])
environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)

        lower_ci_null = []
        for r in diversity_utils.taxa_ranks_with_asv:
            if r == 'OTU':
                lower_ci_null.append(mean_cv_squared_dict[environment]['taxon'][r]['lower_ci_null'])

            else:
                lower_ci_null.append(mean_cv_squared_dict[environment]['taxon'][r]['lower_ci_null'])

        lower_ci_null = numpy.asarray(lower_ci_null)
        #lower_ci_null = numpy.asarray([mean_cv_squared_dict[environment]['taxon'][r]['lower_ci_null'] for r in diversity_utils.taxa_ranks_with_asv])
        upper_ci_null = numpy.asarray([mean_cv_squared_dict[environment]['taxon'][r]['upper_ci_null'] for r in diversity_utils.taxa_ranks_with_asv])
        observed = numpy.asarray([mean_cv_squared_dict[environment]['taxon'][r]['mean_observed'] for r in diversity_utils.taxa_ranks_with_asv])


        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        ax.plot(idx_taxa_ranks, observed, ls='-', lw=1.5, c='k',  zorder=2)
        ax.scatter(idx_taxa_ranks, observed, c=color_all, s=60, label='Observed', linewidth=0.8, edgecolors='k', zorder=3)
        ax.fill_between(idx_taxa_ranks, lower_ci_null, upper_ci_null, color='darkgrey', alpha=0.7, label='95% PI of the null', zorder=1)

        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)
        ax.set_xticks(idx_taxa_ranks)
        ax.set_xticklabels(taxa_ranks_label, fontsize=8)


        if environment_idx == 0:
            ax.set_ylabel("Mean CV^2", fontsize = 10)

        #if environment_chunk_idx == len(environment_chunk_all)-1:
        #    ax.set_xlabel("Phylogenetic dist log relative abundance", fontsize = 10)

        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="upper right", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%scv_coarse_taxon_null.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()






fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

idx_taxa_ranks = numpy.asarray(list(range(len(diversity_utils.taxa_ranks_label_with_asv))))
taxa_ranks_label = diversity_utils.taxa_ranks_label_with_asv


environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)


        distances = list(mean_cv_squared_dict[environment]['phylo'].keys())
        distances.sort()

        lower_ci_null = numpy.asarray([mean_cv_squared_dict[environment]['phylo'][r]['lower_ci_null'] for r in distances])
        upper_ci_null = numpy.asarray([mean_cv_squared_dict[environment]['phylo'][r]['upper_ci_null'] for r in distances])
        observed = numpy.asarray([mean_cv_squared_dict[environment]['phylo'][r]['mean_observed'] for r in distances])


        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        color_all = numpy.asarray([plot_utils.rgb_blue_phylo(r) for r in distances])

        ax.plot(distances, observed, ls='-', lw=1.5, c='k',  zorder=2)
        ax.scatter(distances, observed, c=color_all, s=60, label='Observed', linewidth=0.8, edgecolors='k', zorder=3)
        ax.fill_between(distances, lower_ci_null, upper_ci_null, color='darkgrey', alpha=0.7, label='95% PI of the null', zorder=1)
        ax.set_xscale('log', base=10)
        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)

        if environment_idx == 0:
            ax.set_ylabel("Mean CV^2", fontsize = 10)

        if environment_chunk_idx == len(environment_chunk_all)-1:
            ax.set_xlabel("Phylogenetic distance", fontsize = 10)

        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="lower left", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%scv_coarse_phylo_null.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
