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

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label



taxa_ranks.insert(0, 'ASV')

environment = 'human gut metagenome'


fig = plt.figure(figsize = (8, 12.5)) #
fig.subplots_adjust(bottom= 0.15)

transfers = [12,18]

#ax_occupancy = plt.subplot2grid((1, 3), (0,0), colspan=1)
ax_asv = plt.subplot2grid((3,2), (0,0), colspan=1)
ax_genus = plt.subplot2grid((3,2), (0,1), colspan=1)

ax_family = plt.subplot2grid((3,2), (1,0), colspan=1)
ax_order = plt.subplot2grid((3,2), (1,1), colspan=1)

ax_class = plt.subplot2grid((3,2), (2,0), colspan=1)
ax_phylum = plt.subplot2grid((3,2), (2,1), colspan=1)






ax_all = [ax_asv, ax_genus, ax_family, ax_order, ax_class, ax_phylum]



sys.stderr.write("Subsetting samples for %s...\n" % environment)
#samples = diversity_utils.subset_observations(environment=environment)
sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment)
samples = sad_annotated_dict['samples']

taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])




for rank_idx, rank in enumerate(taxa_ranks):

    sys.stderr.write("Starting %s level analysis...\n" % rank)

    if rank == 'ASV':

        s_by_s_genera = s_by_s

    else:

        all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
        all_genera_idx = numpy.arange(len(all_genera))

        genus_to_taxa_dict = {}
        sad_genera_all = []
        for genus_idx, genus in enumerate(all_genera):
            genus_to_taxa_dict[genus] = []
            for t in taxa:
                if sad_annotated_dict['taxa'][t][rank] == genus:
                    genus_to_taxa_dict[genus].append(t)


            g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])

            s_by_s_genus = s_by_s[g_taxa_idx,:]
            sad_genera_all.append(s_by_s_genus.sum(axis=0))

        s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
        # remove sites where there are no observations
        s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]

    rel_s_by_s_genera = s_by_s_genera/numpy.sum(s_by_s_genera,axis=0)

    mean_genera = numpy.mean(rel_s_by_s_genera, axis=1)
    var_genera = numpy.var(rel_s_by_s_genera, axis=1)

    #slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(mean_genera), numpy.log10(var_genera))

    occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s_genera, range(s_by_s_genera.shape[0]))

    idx_to_keep = (occupancies>0) & (predicted_occupancies > 0) & (mad > 0)
    mad = mad[idx_to_keep]
    occupancies = occupancies[idx_to_keep]
    predicted_occupancies = predicted_occupancies[idx_to_keep]


    ax = ax_all[rank_idx]

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
        bins_mean_all_to_keep.append(bins_mean_all[i])
        bins_occupancies.append(numpy.mean(predicted_occupancies_log10_i))


    bins_mean_all_to_keep = numpy.asarray(bins_mean_all_to_keep)
    bins_occupancies = numpy.asarray(bins_occupancies)

    bins_mean_all_to_keep_no_nan = bins_mean_all_to_keep[(~numpy.isnan(bins_mean_all_to_keep)) & (~numpy.isnan(bins_occupancies))]
    bins_occupancies_no_nan = bins_occupancies[(~numpy.isnan(bins_mean_all_to_keep)) & (~numpy.isnan(bins_occupancies))]

    ax.plot(10**bins_mean_all_to_keep_no_nan, 10**bins_occupancies_no_nan, lw=3, ls='--',c='k', zorder=2, label='Prediction')





    #ax.set_xlim([min_*0.5, 1.08])
    #ax.set_ylim([min_*0.5, 1.08])

    ax.set_xscale('log', base=10)
    ax.set_yscale('log', base=10)
    ax.set_xlabel('Mean relative abundance', fontsize=11)
    ax.set_ylabel('Occupancy', fontsize=11)
    ax.tick_params(axis='both', which='minor', labelsize=9)
    ax.tick_params(axis='both', which='major', labelsize=9)

    ax.set_title(diversity_utils.taxa_ranks_label_with_asv[rank_idx], fontsize=12)

    if rank_idx == 0:
        ax.legend(loc="upper left", fontsize=7)


    #ax_occupancy.legend(loc="upper left", fontsize=8)





fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sabundance_occupancy.png" % config.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
