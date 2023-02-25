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
import matplotlib as mpl
import dbd_utils

import diversity_utils
import config
import plot_utils


rarefied = False


taxa_ranks = diversity_utils.taxa_ranks_with_asv
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label_with_asv




fig = plt.figure(figsize = (8, 4)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

ax_taxon = plt.subplot2grid((1, 2), (0,0))
ax_phylo = plt.subplot2grid((1, 2), (0,1))



coarse_grained_tree_dict_all = {}
distances_all = []
for environment in diversity_utils.environments_to_keep:
    sys.stderr.write("Loading tree dict for %s...\n" % environment)
    coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)
    coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict




slope_taxon_dict = {}
slope_phylo_dict = {}
#
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)

    samples = sad_annotated_dict['samples']
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    sys.stderr.write("Getting site-by-species matrix...\n")

    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

    sys.stderr.write("Running taxonomic coarse-graining.....\n")
    for rank_idx, rank in enumerate(taxa_ranks):

        if rank not in slope_taxon_dict:
            slope_taxon_dict[rank] = []

        if rank == 'OTU':

            s_by_s_genera = s_by_s

        else:


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

        rel_s_by_s_genera = s_by_s_genera/s_by_s_genera.sum(axis=0)

        mean_rel_s_by_s_genera = numpy.mean(rel_s_by_s_genera, axis=1)
        var_rel_s_by_s_genera = numpy.var(rel_s_by_s_genera, axis=1)

        slope_taxon, intercept_taxon, r_value_taxon, p_value_taxon, std_err_taxon = stats.linregress(numpy.log10(mean_rel_s_by_s_genera), numpy.log10(var_rel_s_by_s_genera))
        ax_taxon.scatter(rank_idx, slope_taxon, color=color, s=40, linewidth=0.8, edgecolors='k',  zorder=2)

        slope_taxon_dict[rank].append(slope_taxon)



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
    sys.stderr.write("Running phylogenetic coarse-graining.....\n")
    for distance_idx, distance in enumerate(distances):

        sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
        coarse_grained_list = coarse_grained_tree_dict[distance]
        coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

        # get indexes for each clade
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)

        if distance not in slope_phylo_dict:
            slope_phylo_dict[distance] = []

        color = color_environment_phylo[distance_idx+1]

        rel_s_by_s_all_clades = s_by_s_all_clades/s_by_s_all_clades.sum(axis=0)

        mean_rel_s_by_s_all_clades = numpy.mean(rel_s_by_s_all_clades, axis=1)
        var_rel_s_by_s_all_clades = numpy.var(rel_s_by_s_all_clades, axis=1)

        slope_phylo, intercept_phylo, r_value_phylo, p_value_phylo, std_err_phylo = stats.linregress(numpy.log10(mean_rel_s_by_s_all_clades), numpy.log10(var_rel_s_by_s_all_clades))
        ax_phylo.scatter(distance, slope_phylo, color=color, zorder=2, s=40, linewidth=0.8, edgecolors='k')

        slope_phylo_dict[distance].append(slope_phylo)





ax_taxon.set_xticks(idx_taxa_ranks)
ax_taxon.set_xticklabels(taxa_ranks_label, fontsize=9.5)
ax_taxon.set_ylabel("Taylor's Law slope", fontsize = 12)

# plot mean error
mean_slope = [numpy.mean(slope_taxon_dict[rank]) for rank in taxa_ranks]
ax_taxon.plot(idx_taxa_ranks, mean_slope, lw=2,ls='--',c='k', zorder=3, label='Mean across environments')
ax_taxon.axhline(2, lw=1.5, ls=':',color='k', label=r'$\alpha =2$', zorder=1)
#ax_taxon.set_ylim([1,2.05])
legend = ax_taxon.legend(loc="upper right", fontsize=7)
legend.get_frame().set_alpha(None)

# Phylo


# plot mean error
distance_all = list(slope_phylo_dict.keys())
distance_all.sort()
mean_slope = [numpy.mean(slope_phylo_dict[distance]) for distance in distance_all]
ax_phylo.plot(distance_all, mean_slope, lw=2,ls='--',c='k', zorder=3, label='Mean across environments')

ax_phylo.set_xscale('log', base=10)
ax_phylo.set_xlabel("Phylogenetic distance", fontsize = 12)
ax_phylo.set_ylabel("Taylor's Law slope", fontsize = 12)
ax_phylo.axhline(2, lw=1.5, ls=':',color='k', label=r'$\alpha =2$', zorder=1)
#ax_phylo.set_ylim([1,2.05])





fig.subplots_adjust(hspace=0.2,wspace=0.25)
fig_name = "%staylors_law_summary.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
