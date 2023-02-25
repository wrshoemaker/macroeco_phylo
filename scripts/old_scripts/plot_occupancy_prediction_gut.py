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

    slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(mean_genera), numpy.log10(var_genera))

    print(slope)



    occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s_genera, range(s_by_s_genera.shape[0]))

    idx_to_keep = (occupancies>0) & (predicted_occupancies > 0)
    occupancies = occupancies[idx_to_keep]
    predicted_occupancies = predicted_occupancies[idx_to_keep]


    ax = ax_all[rank_idx]

    predicted_and_observed_occupancies = numpy.concatenate((occupancies,predicted_occupancies),axis=0)
    min_ = min(predicted_and_observed_occupancies)
    max_ = max(predicted_and_observed_occupancies)

    sorted_plot_data = plot_utils.plot_color_by_pt_dens(occupancies, predicted_occupancies, radius=plot_utils.color_radius, loglog=1)
    x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]


    ax.scatter(x, y, c=numpy.sqrt(z), cmap='Blues', s=70, alpha=0.9, edgecolors='none', zorder=1)
    all_ = numpy.concatenate([x, y])



    # occupancy
    ax.plot([min_*0.5,1.08],[min_*0.5,1.08], lw=2,ls='--',c='k',zorder=2, label='1:1')

    #ax.scatter(occupancies, predicted_occupancies, alpha=0.08, c='k', s=12, zorder=1)#, linewidth=0.8, edgecolors='k')


    ax.set_xlim([min_*0.5, 1.08])
    ax.set_ylim([min_*0.5, 1.08])

    ax.set_xscale('log', base=10)
    ax.set_yscale('log', base=10)
    ax.set_xlabel('Observed occupancy', fontsize=11)
    ax.set_ylabel('Predicted occupancy, gamma', fontsize=11)
    ax.tick_params(axis='both', which='minor', labelsize=9)
    ax.tick_params(axis='both', which='major', labelsize=9)

    ax.set_title(diversity_utils.taxa_ranks_label_with_asv[rank_idx], fontsize=12)

    if rank_idx == 0:
        ax.legend(loc="upper left", fontsize=7)


    #ax_occupancy.legend(loc="upper left", fontsize=8)





fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%soccupancy_no_subsampling_gut.png" % config.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
