
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

from scipy.special import rel_entr



rarefied = False

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label
taxa_ranks.insert(0, 'ASV')



fig = plt.figure(figsize = (12, 20)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


for environment_idx, environment in enumerate(environments_to_keep):

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)
    samples = sad_annotated_dict['samples']

    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    for rank_idx, rank in enumerate(taxa_ranks):

        print(rank)

        if rank == 'ASV':

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


        rel_s_by_s_genera = s_by_s_genera/s_by_s_genera.sum(axis=0)
        mean_all = numpy.mean(rel_s_by_s_genera, axis=1)
        var_all = numpy.var(rel_s_by_s_genera, axis=1)

        ax = plt.subplot2grid((9, 6), (environment_idx, rank_idx))

        sorted_plot_data = plot_utils.plot_color_by_pt_dens(mean_all, var_all, radius=plot_utils.color_radius, loglog=1)
        x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]
        ax.scatter(x, y, c=numpy.sqrt(z), cmap='Blues', s=70, alpha=0.9, edgecolors='none', zorder=1)

        slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(mean_all), numpy.log10(var_all))

        x_log10_range =  numpy.linspace(min(numpy.log10(mean_all)) , max(numpy.log10(mean_all)) , 10000)
        y_log10_fit_range = 10 ** (slope*x_log10_range + intercept)
        #y_log10_null_range = 10 ** (slope_null*x_log10_range + intercept)

        ax.plot(10**x_log10_range, y_log10_fit_range, c='k', lw=2.5, linestyle='--', zorder=2, label="OLS regression")

        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)
        #ax.set_xlabel('Mean relative abundance', fontsize=11)
        #ax.set_ylabel('Occupancy', fontsize=11)
        ax.tick_params(axis='both', which='minor', labelsize=9)
        ax.tick_params(axis='both', which='major', labelsize=9)


        if environment_idx == 0:
            ax.set_title(diversity_utils.taxa_ranks_label_with_asv[rank_idx], fontsize=12)

        if rank_idx == 0:
            ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize = 12)

        if (environment_idx == 0) and (rank_idx == 0):
            ax.legend(loc="upper left", fontsize=7)

        #ax.tick_params(axis='x', labelsize=8)
        #ax.tick_params(axis='y', labelsize=8)


        #ax.set_xlabel("Sum of ASV variances", fontsize = 9)
        #ax.set_ylabel("Coarse-grained variance", fontsize = 9)



fig.text(0.3, 0.06, "Mean relative abundance", va='center', fontsize=28)
fig.text(0.02, 0.5, "Variance of relative abundance", va='center', rotation='vertical', fontsize=28)


fig.subplots_adjust(hspace=0.37, wspace=0.37)
fig_name = "%staylors_law_taxon.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
