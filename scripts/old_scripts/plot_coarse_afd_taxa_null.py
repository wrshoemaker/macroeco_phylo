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



rarefied = True
iter = 1000

fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environments_to_keep = diversity_utils.environments_to_keep

environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        print(environment)

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #samples = diversity_utils.subset_observations(environment=environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        sys.stderr.write("Running taxonomic coarse-graining.....\n")
        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)

            taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))

            s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
            rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

            genus_to_taxa_dict = {}
            sad_genera_all = []

            beta_div_mean_ratio_all = []
            for genus in all_genera:
                genus_to_taxa_dict[genus] = []
                for t in taxa:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
                sad_genera_all.append(s_by_s[g_taxa_idx,:].sum(axis=0))


            s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
            # remove sites where there are no observations
            s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]
            rel_s_by_s_genera = (s_by_s_genera/s_by_s_genera.sum(axis=0))

            #clade_log10_rescaled_all = []



            # generate null CIs
            # get indices for null genera coarse graining
            unique_genera, counts_genera = numpy.unique([sad_annotated_dict['taxa'][t][rank] for t in taxa], return_counts=True)
            coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_genera))[:-1]

            s_by_s_copy = numpy.copy(s_by_s)
            clade_log10_rescaled_all = []
            for i in range(100):

                numpy.random.shuffle(s_by_s_copy)
                s_by_s_coarse_null = numpy.add.reduceat(s_by_s_copy, coarse_grain_by_genus_idx, axis=0)
                s_by_s_coarse_null = s_by_s_coarse_null[:,~(numpy.all(s_by_s_coarse_null == 0, axis=0))]
                rel_s_by_s_coarse_null = (s_by_s_coarse_null/s_by_s_coarse_null.sum(axis=0))

                for clade in rel_s_by_s_coarse_null:

                    clade = clade[clade>0]
                    clade_log10 = numpy.log10(clade)

                    if len(clade_log10) < 4:
                        continue

                    clade_log10_rescaled = (clade_log10 - numpy.mean(clade_log10))/numpy.std(clade_log10)
                    clade_log10_rescaled_all.append(clade_log10_rescaled)



            clade_log10_rescaled_all_flat = numpy.concatenate(clade_log10_rescaled_all).ravel()
            hist_, bin_edges_ = numpy.histogram(clade_log10_rescaled_all_flat, density=True, bins=20)
            bins_mean_ = numpy.asarray([0.5 * (bin_edges_[i] + bin_edges_[i+1]) for i in range(0, len(bin_edges_)-1 )])
            hist_to_plot = hist_[hist_>0]
            bins_mean_to_plot = bins_mean_[hist_>0]

            ax.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=diversity_utils.rgb_blue_taxon(rank_idx), alpha=0.9, lw=2, label=rank)


        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)
        ax.set_yscale('log', base=10)

        if environment_idx == 0:
            ax.set_ylabel("Probability density", fontsize = 10)


        if environment_chunk_idx == len(environment_chunk_all)-1:
            ax.set_xlabel("Rescaled log relative abundance", fontsize = 10)

        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="upper left", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%scoarse_afd_taxon_null%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
