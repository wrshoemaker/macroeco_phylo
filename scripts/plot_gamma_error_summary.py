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



rarefied = False


taxa_ranks = diversity_utils.taxa_ranks_with_asv
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label_with_asv


#fig, ax = plt.subplots(figsize=(4,4))
fig = plt.figure(figsize = (8, 4)) #
fig.subplots_adjust(bottom= 0.15)


ax_error_taxon = plt.subplot2grid((1,2), (0,0), colspan=1)
ax_error_phylo = plt.subplot2grid((1,2), (0,1), colspan=1)

ax_error_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_error_taxon.transAxes)
ax_error_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_error_phylo.transAxes)

coarse_grained_tree_dict_all = {}
distances_all = []
for environment in diversity_utils.environments_to_keep:
    sys.stderr.write("Loading tree dict for %s...\n" % environment)
    coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)
    coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict



mean_richness_observed_taxon_all = []
mean_richness_predicted_taxon_all = []
mean_richness_observed_phylo_all = []
mean_richness_predicted_phylo_all = []

error_taxon_all = []
error_phylo_all = []
error_taxon_dict = {}
error_phylo_dict = {}
#
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)

    samples = sad_annotated_dict['samples']
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    sys.stderr.write("Getting site-by-species matrix...\n")

    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

    for rank_idx, rank in enumerate(taxa_ranks):

        if rank not in error_taxon_dict:
            error_taxon_dict[rank] = []

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


        color = color_environment_taxon[rank_idx+1]

        # error
        occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s_genera, s_by_s_genera.shape[0])
        error = diversity_utils.calculate_error(occupancies, predicted_occupancies)
        mean_log_error = numpy.mean(numpy.log10(error))
        error_taxon_dict[rank].append(mean_log_error)
        ax_error_taxon.scatter(rank_idx, mean_log_error, color=color, s=40, linewidth=0.8, edgecolors='k')


        mean_richness_observed, mean_richness_predicted = diversity_utils.predict_mean_richness(s_by_s_genera, s_by_s_genera.shape[0])
        mean_richness_observed_taxon_all.append(mean_richness_observed)
        mean_richness_predicted_taxon_all.append(mean_richness_predicted)



    # phylogenetic
    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)
    samples = sad_annotated_dict['samples']

    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    distances = list(coarse_grained_tree_dict.keys())
    distances.sort()

    color_environment_phylo = plot_utils.get_custom_cmap_phylo(environment, len(distances))

    for distance_idx, distance in enumerate(distances):


        sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
        coarse_grained_list = coarse_grained_tree_dict[distance]
        coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

        # get indexes for each clade
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)

        #if s_by_s_all_clades.shape[0] < 5:
        #    continue


        if distance not in error_phylo_dict:
            error_phylo_dict[distance] = []


        color = color_environment_phylo[distance_idx+1]

        occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s_all_clades, s_by_s_all_clades.shape[0])
        error = diversity_utils.calculate_error(occupancies, predicted_occupancies)
        mean_log_error = numpy.mean(numpy.log10(error))
        error_phylo_dict[distance].append(mean_log_error)
        ax_error_phylo.scatter(distance, mean_log_error, color=color, s=40, linewidth=0.8, edgecolors='k')

        mean_richness_observed, mean_richness_predicted = diversity_utils.predict_mean_richness(s_by_s_all_clades, s_by_s_all_clades.shape[0])
        mean_richness_observed_phylo_all.append(mean_richness_observed)
        mean_richness_predicted_phylo_all.append(mean_richness_predicted)




ax_error_taxon.set_xticks(idx_taxa_ranks)
ax_error_taxon.set_xticklabels(taxa_ranks_label, fontsize=9.5)
ax_error_taxon.set_ylabel("Mean " + r'$\mathrm{log}_{10}$' +  " occupancy error", fontsize = 12)

# plot mean error
mean_mean_error = [numpy.mean(error_taxon_dict[rank]) for rank in taxa_ranks]
ax_error_taxon.plot(idx_taxa_ranks, mean_mean_error, lw=2,ls='--',c='k',zorder=1, label='Mean across environments')
ax_error_taxon.legend(loc="lower left", fontsize=7)
ax_error_taxon.set_xlabel("Taxonomic rank", fontsize=12)


mean_richness_observed_taxon_all = numpy.asarray(mean_richness_observed_taxon_all)
mean_richness_predicted_taxon_all = numpy.asarray(mean_richness_predicted_taxon_all)

predicted_and_observed_richness_taxon = numpy.concatenate((mean_richness_observed_taxon_all,mean_richness_predicted_taxon_all),axis=0)
min_ = min(predicted_and_observed_richness_taxon)*0.5
max_ = max(predicted_and_observed_richness_taxon)*1.1






# Phylo

ax_error_phylo.set_xscale('log', base=10)
ax_error_phylo.set_xlabel("Phylogenetic distance", fontsize = 12)
ax_error_phylo.set_ylabel("Mean " + r'$\mathrm{log}_{10}$' +  " occupancy error", fontsize = 12)

# plot mean error
distance_all = list(error_phylo_dict.keys())
distance_all.sort()
mean_mean_error = [numpy.mean(error_phylo_dict[distance]) for distance in distance_all]
ax_error_phylo.plot(distance_all, mean_mean_error, lw=2,ls='--',c='k',zorder=1, label='Mean across environments')


mean_richness_observed_phylo_all = numpy.asarray(mean_richness_observed_phylo_all)
mean_richness_predicted_phylo_all = numpy.asarray(mean_richness_predicted_phylo_all)

predicted_and_observed_richness_phylo = numpy.concatenate((mean_richness_observed_phylo_all,mean_richness_predicted_phylo_all),axis=0)
min_ = min(predicted_and_observed_richness_phylo)*0.5
max_ = max(predicted_and_observed_richness_phylo)*1.1






fig.subplots_adjust(hspace=0.2,wspace=0.25)
fig_name = "%sgamma_error_summary.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
