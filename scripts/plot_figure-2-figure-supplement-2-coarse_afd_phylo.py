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


fig = plt.figure(figsize = (11, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

min_occupancy = 1.0

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

# keep same color scale for all plots
# get all the distances in single plot


coarse_grained_tree_dict_all = {}
distances_all = []
for environment in environments_to_keep:
    sys.stderr.write("Loading tree dict for %s...\n" % environment)
    coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)
    coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict
    distances = list(coarse_grained_tree_dict.keys())
    distances_all.extend(distances)

distances_all = list(set(distances_all))
distances_all.sort()
color_all = plot_utils.make_blue_cmap(len(distances_all))




environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
ax_all = []
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]

        distances = list(coarse_grained_tree_dict.keys())
        distances.sort()

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #samples = diversity_utils.subset_observations(environment=environment)
        pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)

        samples = pres_abs_dict['samples']
        taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
        sys.stderr.write("Getting site-by-species matrix...\n")
        #s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)

        s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))
        ax_all.append(ax)
        sys.stderr.write("Running phylogenetic coarse-graining.....\n")

        asv_log10_rescaled_all_to_fit = []
        for distance_idx, distance in enumerate(distances):

            sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
            coarse_grained_list = coarse_grained_tree_dict[distance]
            coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

            #if len(coarse_grained_n) == len(taxonomy_names):
            #    continue

            # get indexes for each clade
            coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

            # coarse grain s-by-s for all clades
            s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
            rel_s_by_s_all_clades = (s_by_s_all_clades/s_by_s_all_clades.sum(axis=0))

            occupancy_clades = (s_by_s_all_clades>0).sum(axis=1)/s_by_s_all_clades.shape[1]


            clades_to_keep = (occupancy_clades >= min_occupancy)
            
            rel_s_by_s_all_clades = rel_s_by_s_all_clades[clades_to_keep,:]


            clade_log10_rescaled_all = []
            for clade in rel_s_by_s_all_clades:

                clade = clade[clade>0]
                clade_log10 = numpy.log10(clade)

                if len(clade_log10) < 4:
                    continue

                clade_log10_rescaled = (clade_log10 - numpy.mean(clade_log10))/numpy.std(clade_log10)
                clade_log10_rescaled_all.extend(clade_log10_rescaled)

            asv_log10_rescaled_all_to_fit.extend(clade_log10_rescaled_all)
            hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(clade_log10_rescaled_all)
            ax.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=color_all(distances_all.index(distance)), alpha=0.9, lw=2)


        asv_log10_rescaled_all_to_fit = numpy.asarray(asv_log10_rescaled_all_to_fit)
        shape_gamma, loc_gamma, scale_gamma = stats.loggamma.fit(asv_log10_rescaled_all_to_fit)
        x = numpy.linspace(stats.loggamma.ppf(0.0001, shape_gamma, loc=loc_gamma, scale=scale_gamma), stats.loggamma.ppf(0.9999, shape_gamma, loc=loc_gamma, scale=scale_gamma), 100)
        pdf_loggamma_to_plot = stats.loggamma.pdf(x, shape_gamma, loc=loc_gamma, scale=scale_gamma)
        ax.plot(x, pdf_loggamma_to_plot, 'k', ls='--', lw=3, label='Gamma fit')
        #ax.legend(loc="upper left", fontsize=9, frameon=False)
        ax.set_ylim([0.0004, 1.5])


        ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)
        ax.set_yscale('log', base=10)

        if environment_idx == 0:
            ax.set_ylabel("Probability density", fontsize = 10)


        if environment_chunk_idx == len(environment_chunk_all)-1:
            ax.set_xlabel("Rescaled log relative abundance", fontsize = 10)

        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="upper left", fontsize=9, frameon=False)



norm = mpl.colors.LogNorm(vmin=min(distances_all), vmax=max(distances_all))
pcm = plt.cm.ScalarMappable(cmap=color_all, norm=norm)

ax_all = numpy.array(ax_all)
#ax_all = numpy.reshape(ax_all, (3, 3))  # or

fig.subplots_adjust(hspace=0.2,wspace=0.25)
clb = fig.colorbar(pcm, ax=ax_all, shrink=0.9, pad = 0.04)
clb.set_label(label='Phylogenetic distance', fontsize=12)



fig_name = "%sfigure-2-figure-supplement-2-coarse_afd_phylo%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
