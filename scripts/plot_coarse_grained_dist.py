import os
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
from collections import Counter
import matplotlib.gridspec as gridspec

import dbd_utils

import diversity_utils
import config
import plot_utils


rarefied = False
iter = 1000

#fig = plt.figure(figsize = (10, 10)) #
#fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


fig = plt.figure(figsize = (9.5, 4))
gs = gridspec.GridSpec(nrows=1, ncols=2)

ax_taxon = fig.add_subplot(gs[0, 0])
ax_phylo = fig.add_subplot(gs[0, 1])

ax_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_taxon.transAxes)
ax_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_phylo.transAxes)



environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label



coarse_grained_tree_dict_all = {}
distances_all = []
for environment in environments_to_keep:
    sys.stderr.write("Loading tree dict for %s...\n" % environment)
    coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)
    coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict
    distances = list(coarse_grained_tree_dict.keys())
    distances_all.extend(distances)



color_all = numpy.asarray([plot_utils.rgb_blue_taxon(r) for r in idx_taxa_ranks])



coarse_grained_fract = {}

# ........
environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #samples = diversity_utils.subset_observations(environment=environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

        # dictionary to map taxon names to higher order rank
        sys.stderr.write("Running taxonomic coarse-graining.....\n")

        fraction_taxa = []

        for rank_idx, rank in enumerate(taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)
            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            #all_genera_idx = numpy.arange(len(all_genera))
            fraction_remaining = len(all_genera)/ len(taxa)
            fraction_taxa.append(fraction_remaining)

            if rank not in coarse_grained_fract:
                coarse_grained_fract[rank] = []

            coarse_grained_fract[rank].append(1-fraction_remaining)

            #all_genera = [sad_annotated_dict['taxa'][t][rank] for t in taxa]
            #all_genera_counts = numpy.asarray(list(dict(Counter(all_genera)).values()))

            #all_genera_counts_sort = numpy.sort(all_genera_counts)

            #all_genera_counts_sort_range = numpy.arange(1, max(all_genera_counts_sort)+1)

            #cdf = [sum(all_genera_counts_sort[all_genera_counts_sort>=i])/sum(all_genera_counts_sort) for i in range(max(all_genera_counts_sort))]
            #cdf = numpy.asarray(cdf)

            #ax.plot(all_genera_counts_sort_range, cdf, c=plot_utils.rgb_blue_taxon(rank_idx), ls='-', lw=2, alpha=1)



        # get counts for phylo
        coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]
        pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)

        samples = pres_abs_dict['samples']
        taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))

        sys.stderr.write("Getting site-by-species matrix...\n")
        s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

        distances = list(coarse_grained_tree_dict.keys())
        distances.sort()

        sys.stderr.write("Running phylogenetic coarse-graining.....\n")
        fraction_phylo = []
        for distance_idx, distance in enumerate(distances):

            sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))

            coarse_grained_list = coarse_grained_tree_dict[distance]

            fraction_phylo.append(len(coarse_grained_list) / len(taxa))


        fraction_taxa = numpy.asarray(fraction_taxa)
        fraction_phylo = numpy.asarray(fraction_phylo)

        fraction_taxa = 1-fraction_taxa
        fraction_phylo = 1-fraction_phylo

        c = plot_utils.environment_color_map[environment]

        #ax_taxon.plot(idx_taxa_ranks, fraction_taxa, ls='-', lw=2, alpha=0.8, c=c)
        #ax_phylo.plot(distances, fraction_phylo, ls='-', lw=2, alpha=0.8, c=c, label=diversity_utils.format_environment_label(environment))

        ax_taxon.plot(idx_taxa_ranks, fraction_taxa, ls='-', lw=1.2, alpha=0.8, c='k')
        ax_phylo.plot(distances, fraction_phylo, ls='-', lw=1.2, alpha=0.8, c='k', zorder=2)

        ax_taxon.scatter(idx_taxa_ranks, fraction_taxa, color=c, s=45, linewidth=0.8, edgecolors='k', zorder=3, label=diversity_utils.format_environment_label(environment))
        ax_phylo.scatter(distances, fraction_phylo, color=c, s=45, linewidth=0.8, edgecolors='k', zorder=3)




for rank_idx, rank in enumerate(taxa_ranks):

    mean_fraction_coarseg_grained = numpy.mean(coarse_grained_fract[rank])
    color_rank = color_all[rank_idx]

    ax_phylo.axhline(y=mean_fraction_coarseg_grained, color=color_rank, linestyle=':', lw = 2, zorder=1, label=taxa_ranks_label[rank_idx])
    print(rank, mean_fraction_coarseg_grained)






ax_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_taxon.transAxes)
ax_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_phylo.transAxes)



#ax_taxon.set_yscale('log', base=10)
ax_taxon.set_xticks(idx_taxa_ranks)
ax_taxon.set_xticklabels(taxa_ranks_label, fontsize=10)

#ax_taxon.set_ylabel("# coarse-grained taxa/ # ASVs", fontsize=12)
ax_taxon.set_ylabel('Fraction of coarse-grained OTUs', fontsize=12)
ax_taxon.legend(loc="lower right", fontsize=7)
ax_taxon.set_xlabel("Taxonomic rank", fontsize=12)




ax_phylo.set_xscale('log', base=10)
#ax_phylo.set_yscale('log', base=10)

ax_phylo.set_xlabel("Phylogenetic distance", fontsize=12)
ax_phylo.set_ylabel('Fraction of coarse-grained OTUs', fontsize=12)
ax_phylo.legend(loc="lower right", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%scoarse_grained_dist%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
