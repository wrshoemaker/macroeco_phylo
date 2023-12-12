import os
#from Bio import Phylo
import random
import copy
import config
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
import plot_utils

import mle_utils


rarefied = False


fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


min_occupancy = 1.0
environments_to_keep = diversity_utils.environments_to_keep

environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #samples = diversity_utils.subset_observations(environment=environment)
        pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)
        samples = pres_abs_dict['samples']

        taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
        s_by_s_taxon = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s_taxon = (s_by_s_taxon/s_by_s_taxon.sum(axis=0))

        occupancy_taxon = (s_by_s_taxon>0).sum(axis=1)/s_by_s_taxon.shape[1]
        
        if environment == 'microbial mat metagenome':
            taxon_to_keep = (occupancy_taxon >= 0.85)
        else:
            taxon_to_keep = (occupancy_taxon >= min_occupancy)
        
        #taxon_to_keep = (occupancy_taxon >= min_occupancy)
        rel_s_by_s_taxon = rel_s_by_s_taxon[taxon_to_keep,:]


        asv_log10_rescaled_all = []
        for clade in rel_s_by_s_taxon:

            clade = clade[clade>0]
            clade_log10 = numpy.log10(clade)

            if len(clade_log10) < 5:
                continue

            clade_log10_rescaled = (clade_log10 - numpy.mean(clade_log10))/numpy.std(clade_log10)
            asv_log10_rescaled_all.extend(clade_log10_rescaled.tolist())

        hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(asv_log10_rescaled_all)

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))
        ax.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=plot_utils.rgb_blue_taxon(0), alpha=0.9, lw=2, label='OTU')


        asv_log10_rescaled_all_to_fit = []

        sys.stderr.write("Running taxonomic coarse-graining.....\n")
        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)

            all_genera = numpy.asarray(list(set([pres_abs_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))


            genus_to_taxa_dict = {}
            sad_genera_all = []

            beta_div_mean_ratio_all = []
            for genus in all_genera:
                genus_to_taxa_dict[genus] = []
                for t in taxa:
                    if pres_abs_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
                sad_genera_all.append(s_by_s_taxon[g_taxa_idx,:].sum(axis=0))


            s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
            # remove sites where there are no observations
            s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]
            rel_s_by_s_genera = (s_by_s_genera/s_by_s_genera.sum(axis=0))

            occupancy_genera = (s_by_s_genera>0).sum(axis=1)/s_by_s_genera.shape[1]

            print(environment)
            if environment == 'microbial mat metagenome':
                genera_to_keep = (occupancy_genera >= 0.75)
            else:
                genera_to_keep = (occupancy_genera >= min_occupancy)

            #genera_to_keep = (occupancy_genera >= min_occupancy)
            rel_s_by_s_genera = rel_s_by_s_genera[genera_to_keep,:]

            clade_log10_rescaled_all = []
            for clade in rel_s_by_s_genera:

                clade = clade[clade>0]
                clade_log10 = numpy.log10(clade)

                if len(clade_log10) < 20:
                    continue

                clade_log10_rescaled = (clade_log10 - numpy.mean(clade_log10))/numpy.std(clade_log10)
                clade_log10_rescaled_all.extend(clade_log10_rescaled)


            #clade_log10_rescaled_all_flat = numpy.concatenate(clade_log10_rescaled_all).ravel()
            #hist_, bin_edges_ = numpy.histogram(clade_log10_rescaled_all_flat, density=True, bins=10)
            #bins_mean_ = numpy.asarray([0.5 * (bin_edges_[i] + bin_edges_[i+1]) for i in range(0, len(bin_edges_)-1 )])
            #hist_to_plot = hist_[hist_>0]
            #bins_mean_to_plot = bins_mean_[hist_>0]

            asv_log10_rescaled_all_to_fit.extend(clade_log10_rescaled_all)
            hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(clade_log10_rescaled_all)

            ax.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=plot_utils.rgb_blue_taxon(rank_idx+1), alpha=0.9, lw=2, label=rank)

        asv_log10_rescaled_all_to_fit = numpy.asarray(asv_log10_rescaled_all_to_fit)
        shape_gamma, loc_gamma, scale_gamma = stats.loggamma.fit(asv_log10_rescaled_all_to_fit)
        x = numpy.linspace(stats.loggamma.ppf(0.0001, shape_gamma, loc=loc_gamma, scale=scale_gamma), stats.loggamma.ppf(0.9999, shape_gamma, loc=loc_gamma, scale=scale_gamma), 100)
        pdf_loggamma_to_plot = stats.loggamma.pdf(x, shape_gamma, loc=loc_gamma, scale=scale_gamma)
        ax.plot(x, pdf_loggamma_to_plot, 'k', ls='--', lw=3, label='Gamma fit')
        #ax.legend(loc="upper left", fontsize=9, frameon=False)
        ax.set_ylim([0.0004, 1.7])


        ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)
        ax.set_yscale('log', base=10)

        if environment_idx == 0:
            ax.set_ylabel("Probability density", fontsize = 10)


        if environment_chunk_idx == len(environment_chunk_all)-1:
            ax.set_xlabel("Rescaled log relative abundance", fontsize = 10)

        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="upper left", fontsize=7, frameon=False)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sfigure-2-figure-supplement-1-coarse_afd_taxon%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
