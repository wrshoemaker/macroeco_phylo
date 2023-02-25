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

import dbd_utils

import diversity_utils
import config
import plot_utils


rarefied = False
iter = 1000

fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label


color_all = numpy.asarray([plot_utils.rgb_blue_taxon(r) for r in idx_taxa_ranks])
environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #samples = diversity_utils.subset_observations(environment=environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

        # dictionary to map taxon names to higher order rank
        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))
        sys.stderr.write("Running taxonomic coarse-graining.....\n")
        beta_div_mean_dict = {}
        for rank_idx, rank in enumerate(taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)
            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))

            genus_to_taxa_dict = {}
            sad_genera_all = []
            for genus_idx, genus in enumerate(all_genera):
                genus_to_taxa_dict[genus] = []
                for t in taxa:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                        if t not in beta_div_mean_dict:
                            beta_div_mean_dict[t] = {}

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])

                s_by_s_genus = s_by_s[g_taxa_idx,:]
                sad_genera_all.append(s_by_s_genus.sum(axis=0))

                rel_s_by_s_genus = rel_s_by_s[g_taxa_idx,:]
                mean_genus = numpy.mean(rel_s_by_s_genus, axis=1)
                std_genus = numpy.std(rel_s_by_s_genus, axis=1)
                beta_div_mean = mean_genus/(std_genus**2)

                if len(beta_div_mean) >=5:
                    cv_beta_div_mean = numpy.std(beta_div_mean)/numpy.mean(beta_div_mean)
                else:
                    cv_beta_div_mean = float('nan')

                for t in genus_to_taxa_dict[genus]:
                    beta_div_mean_dict[t][rank] = cv_beta_div_mean


            #s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
            # remove sites where there are no observations
            #s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]
            #rel_s_by_s_genera = (s_by_s_genera/s_by_s_genera.sum(axis=0))


            #mean_genera = numpy.mean(rel_s_by_s_genera, axis=1)
            #std_genera = numpy.std(rel_s_by_s_genera, axis=1)

            #beta_div_mean_genera = mean_genera/(std_genera**2)

            #for genus_idx, genus in enumerate(all_genera):
            #    for t_idx, t in enumerate(genus_to_taxa_dict[genus]):
            #        beta_div_mean_dict[t][rank] = beta_div_mean_genera[genus_idx]


        #species = list(beta_div_mean_dict.keys())
        #ranks_in_dict = ['otu', 'genus', 'family', 'order',  'class',  'phylum']
        #for s in species:
        #    beta_div_mean_to_plot = [beta_div_mean_dict[s][r] for r in ranks_in_dict]
        #    ax.plot(idx_taxa_ranks, beta_div_mean_to_plot, ls='-', lw=0.75, c='dodgerblue', alpha=0.005,  zorder=2)


        genera_plotted = []
        species = list(beta_div_mean_dict.keys())
        for s in species:

            if beta_div_mean_dict[s]['genus'] in genera_plotted:
                continue

            beta_div_mean = numpy.asarray([beta_div_mean_dict[s][r] for r in taxa_ranks])
            idx_to_keep = ~numpy.isnan(beta_div_mean)

            if sum(idx_to_keep) < 5:
                continue

            idx_taxa_ranks_to_plot = idx_taxa_ranks[idx_to_keep]
            beta_div_mean_to_plot = beta_div_mean[idx_to_keep]
            ax.plot(idx_taxa_ranks_to_plot, beta_div_mean_to_plot, ls='-', lw=0.75, c='k', alpha=0.1,  zorder=2)

            genera_plotted.append(beta_div_mean_dict[s]['genus'])

        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)

        ax.set_xticks(idx_taxa_ranks)
        ax.set_xticklabels(taxa_ranks_label, fontsize=8)

        ax.axhline(y=1, color='k', linestyle=':', lw = 2, zorder=1)

        ax.set_ylim([0,1.55])

        if environment_idx == 0:
            ax.set_ylabel("CV of gamma rate parameter", fontsize = 10)

        #if environment_chunk_idx == len(environment_chunk_all)-1:
        #    ax.set_xlabel("Rescaled log relative abundance", fontsize = 10)

        #if (environment_chunk_idx ==0) and (environment_idx == 0):
        #    ax.legend(loc="upper left", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%scv_beta_div_mean_taxon_emp%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
