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



        occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s_genera, s_by_s_genera.shape[0])

        idx_to_keep = (occupancies>0) & (predicted_occupancies > 0) & (mad > 0)
        mad = mad[idx_to_keep]
        occupancies = occupancies[idx_to_keep]
        predicted_occupancies = predicted_occupancies[idx_to_keep]

        ax = plt.subplot2grid((9, 6), (environment_idx, rank_idx))

        mad_occupancy_joint = numpy.concatenate((mad, occupancies),axis=0)
        min_ = min(mad_occupancy_joint)
        max_ = max(mad_occupancy_joint)

        sorted_plot_data = plot_utils.plot_color_by_pt_dens(mad, occupancies, radius=plot_utils.color_radius, loglog=1)
        x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]


        ax.scatter(x, y, c=numpy.sqrt(z), cmap='Blues', s=70, alpha=0.9, edgecolors='none', zorder=1)
        all_ = numpy.concatenate([x, y])


        # mad vs occupancy
        mad_log10 = numpy.log10(mad)
        occupancies_log10 = numpy.log10(occupancies)
        predicted_occupancies_log10 = numpy.log10(predicted_occupancies)
        hist_all, bin_edges_all = numpy.histogram(mad_log10, density=True, bins=25)
        bins_mean_all = [0.5 * (bin_edges_all[i] + bin_edges_all[i+1]) for i in range(0, len(bin_edges_all)-1 )]
        bins_mean_all_to_keep = []
        bins_occupancies = []
        for i in range(0, len(bin_edges_all)-1 ):
            predicted_occupancies_log10_i = predicted_occupancies_log10[(mad_log10>=bin_edges_all[i]) & (mad_log10<bin_edges_all[i+1])]
            #bins_mean_all_to_keep.append(bins_mean_all[i])
            bins_mean_all_to_keep.append(bin_edges_all[i])
            bins_occupancies.append(numpy.mean(predicted_occupancies_log10_i))


        bins_mean_all_to_keep = numpy.asarray(bins_mean_all_to_keep)
        bins_occupancies = numpy.asarray(bins_occupancies)

        bins_mean_all_to_keep_no_nan = bins_mean_all_to_keep[(~numpy.isnan(bins_mean_all_to_keep)) & (~numpy.isnan(bins_occupancies))]
        bins_occupancies_no_nan = bins_occupancies[(~numpy.isnan(bins_mean_all_to_keep)) & (~numpy.isnan(bins_occupancies))]

        ax.plot(10**bins_mean_all_to_keep_no_nan, 10**bins_occupancies_no_nan, lw=3, ls='--',c='k', zorder=2, label='Gamma prediction')

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
            ax.legend(loc="upper left", fontsize=6)

        #ax.tick_params(axis='x', labelsize=8)
        #ax.tick_params(axis='y', labelsize=8)


        #ax.set_xlabel("Sum of ASV variances", fontsize = 9)
        #ax.set_ylabel("Coarse-grained variance", fontsize = 9)


fig.text(0.3, 0.06, "Mean relative abundance", va='center', fontsize=28)
fig.text(0.02, 0.5, "Occupancy", va='center', rotation='vertical', fontsize=28)


fig.subplots_adjust(hspace=0.37, wspace=0.37)
fig_name = "%sfigure-2-figure-supplement-5-abundance_vs_occupancy_all_taxon%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
