import os
#from Bio import Phylo
import config
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
import itertools
import matplotlib as mpl

import dbd_utils
import diversity_utils
import mle_utils
import plot_utils


environment = 'human gut metagenome'
rarefied = False


#fig = plt.figure(figsize = (9, 12)) #
#fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

fig = plt.figure(figsize = (14, 4)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


#ax_afd_taxon = plt.subplot2grid((3, 2), (0, 0))
ax_afd_phylo = plt.subplot2grid((1,3), (0, 0))

#ax_occupancy_taxon = plt.subplot2grid((3, 2), (1, 0))
ax_occupancy_phylo = plt.subplot2grid((1, 3), (0, 1))

#ax_mad_vs_occupancy_taxon = plt.subplot2grid((3, 2), (2, 0))
ax_mad_vs_occupancy_phylo = plt.subplot2grid((1, 3), (0, 2))


#ax_afd_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_afd_taxon.transAxes)
ax_afd_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_afd_phylo.transAxes)
#ax_occupancy_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[2], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_occupancy_taxon.transAxes)
ax_occupancy_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_occupancy_phylo.transAxes)
#ax_mad_vs_occupancy_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[4], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_mad_vs_occupancy_taxon.transAxes)
ax_mad_vs_occupancy_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[2], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_mad_vs_occupancy_phylo.transAxes)


#ax_afd_taxon.set_title('Taxonomic coarse-graining', fontsize=12, fontweight='bold')
#ax_afd_phylo.set_title('Phylogenetic coarse-graining', fontsize=12, fontweight='bold')

#ax_diversity_dist_taxon.text(1.2, 1.2, "Intra vs. inter-group diversity ", fontsize=18, fontweight='bold', ha='center', va='center', transform=ax_diversity_dist_taxon.transAxes)



occupancy_dict = {}
occupancy_dict['taxon'] = {}
occupancy_dict['phylo'] = {}

sys.stderr.write("Subsetting samples for %s...\n" % environment)
pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)
samples = pres_abs_dict['samples']

taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
s_by_s_taxon = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])
rel_s_by_s_taxon = (s_by_s_taxon/s_by_s_taxon.sum(axis=0))

asv_log10_rescaled_all = []
for clade in rel_s_by_s_taxon:

    clade = clade[clade>0]
    clade_log10 = numpy.log10(clade)

    if len(clade_log10) < 5:
        continue

    clade_log10_rescaled = (clade_log10 - numpy.mean(clade_log10))/numpy.std(clade_log10)
    asv_log10_rescaled_all.extend(clade_log10_rescaled.tolist())


hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(asv_log10_rescaled_all)
#ax_afd_taxon.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=plot_utils.rgb_blue_taxon(0), alpha=0.9, lw=2, label='ASV')

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

    clade_log10_rescaled_all = []
    for clade in rel_s_by_s_genera:

        clade = clade[clade>0]
        clade_log10 = numpy.log10(clade)

        if len(clade_log10) < 50:
            continue

        clade_log10_rescaled = (clade_log10 - numpy.mean(clade_log10))/numpy.std(clade_log10)
        clade_log10_rescaled_all.extend(clade_log10_rescaled)


    occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s_genera, s_by_s_genera.shape[0])
    idx_to_keep = (occupancies>0) & (predicted_occupancies > 0)
    occupancies = occupancies[idx_to_keep]
    predicted_occupancies = predicted_occupancies[idx_to_keep]
    mad = mad[idx_to_keep]
    occupancy_dict['taxon'][rank] = {}
    occupancy_dict['taxon'][rank]['occupancy'] = occupancies
    occupancy_dict['taxon'][rank]['predicted_occupancies'] = predicted_occupancies
    occupancy_dict['taxon'][rank]['mad'] = mad


    hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(clade_log10_rescaled_all)
    #ax_afd_taxon.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=plot_utils.rgb_blue_taxon(rank_idx+1), alpha=0.9, lw=2, label= diversity_utils.taxa_ranks_label[rank_idx])

    #for clade in rel_s_by_s_taxon:

    #    clade = clade[clade>0]
    #    clade_log10 = numpy.log10(clade)

    #    if len(clade_log10) < 80:
    #        continue

    #    clade_log10_rescaled = (clade_log10 - numpy.mean(clade_log10))/numpy.std(clade_log10)
    #    #asv_log10_rescaled_all.extend(clade_log10_rescaled.tolist())
    #    hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(clade_log10_rescaled, bins=15)
    #    ax_afd_taxon.plot(bins_mean_to_plot, hist_to_plot, color=plot_utils.rgb_blue_taxon(rank_idx+1), alpha=0.9, lw=0.4)



#ax_afd_taxon.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)
#ax_afd_taxon.set_yscale('log', base=10)
#ax_afd_taxon.set_ylabel("Probability density", fontsize = 12)
#ax_afd_taxon.set_xlabel("Rescaled log relative abundance", fontsize = 12)
#ax_afd_taxon.legend(loc="upper left", fontsize=7)




# AFD phylo

#coarse_grained_tree_dict_all = {}
#distances_all = []
sys.stderr.write("Loading tree dict for %s...\n" % environment)
coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)
#coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict
#coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]
distances = list(coarse_grained_tree_dict.keys())
distances.sort()
#distances_all.extend(distances)

#distances_all = list(set(distances_all))
#distances_all.sort()

#distances_all = distances_all[::2]
#color_all = plot_utils.make_blue_cmap(len(distances_all))
color_all = plot_utils.make_blue_cmap(len(distances))

sys.stderr.write("Subsetting samples for %s...\n" % environment)
pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)

samples = pres_abs_dict['samples']
taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))

sys.stderr.write("Getting site-by-species matrix...\n")
s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

sys.stderr.write("Running phylogenetic coarse-graining.....\n")
for distance_idx, distance in enumerate(distances):

    sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
    coarse_grained_list = coarse_grained_tree_dict[distance]
    coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

    # get indexes for each clade
    coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

    # coarse grain s-by-s for all clades
    s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
    rel_s_by_s_all_clades = (s_by_s_all_clades/s_by_s_all_clades.sum(axis=0))

    clade_log10_rescaled_all = []
    for clade in rel_s_by_s_all_clades:

        clade = clade[clade>0]
        clade_log10 = numpy.log10(clade)

        if len(clade_log10) < 4:
            continue

        clade_log10_rescaled = (clade_log10 - numpy.mean(clade_log10))/numpy.std(clade_log10)
        clade_log10_rescaled_all.extend(clade_log10_rescaled)


    occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s_all_clades, s_by_s_all_clades.shape[0])
    
    idx_to_keep = (occupancies>0) & (predicted_occupancies > 0)
    occupancies = occupancies[idx_to_keep]
    predicted_occupancies = predicted_occupancies[idx_to_keep]
    mad = mad[idx_to_keep]

    occupancy_dict['phylo'][distance] = {}
    occupancy_dict['phylo'][distance]['occupancy'] = occupancies
    occupancy_dict['phylo'][distance]['predicted_occupancies'] = predicted_occupancies
    occupancy_dict['phylo'][distance]['mad'] = mad


    hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(clade_log10_rescaled_all)
    ax_afd_phylo.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=color_all(distances.index(distance)), alpha=0.9, lw=2)


#ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)
ax_afd_phylo.set_yscale('log', base=10)
ax_afd_phylo.set_ylabel("Probability density", fontsize = 12)
ax_afd_phylo.set_xlabel("Rescaled log relative abundance", fontsize = 12)



norm = mpl.colors.LogNorm(vmin=min(distances), vmax=max(distances))
pcm = plt.cm.ScalarMappable(cmap=color_all, norm=norm)


cbar_ax = fig.add_axes([0.82, 0.15, 0.03, 0.7])
#fig.colorbar(im, cax=cbar_ax)


fig.subplots_adjust(hspace=0.2,wspace=0.25)
fig.subplots_adjust(right=0.8)

#clb = fig.colorbar(pcm, ax=ax_afd_phylo, shrink=0.9, pad = 0.04)
clb = fig.colorbar(pcm, cax=cbar_ax, shrink=0.9, pad = 0.04)
clb.set_label(label='Phylogenetic distance', fontsize=16)


# ax_occupancy_taxon
occupancy_all = []
predicted_occupancies_all = []
for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

    occupancy = occupancy_dict['taxon'][rank]['occupancy']
    predicted_occupancies = occupancy_dict['taxon'][rank]['predicted_occupancies']

    log10_occupancy = numpy.log10(occupancy)
    log10_predicted_occupancies = numpy.log10(predicted_occupancies)
    hist_, bin_edges_ = numpy.histogram(log10_occupancy, density=True, bins=20)

    hist_to_plot = []
    log10_predicted_occupancies_to_plot = []
    for i in range(0, len(bin_edges_)-1):

        log10_predicted_occupancies_i = log10_predicted_occupancies[(log10_occupancy >= bin_edges_[i]) & (log10_occupancy < bin_edges_[i+1])]
        if len(log10_predicted_occupancies_i) < 5:
            continue

        hist_to_plot.append(hist_[i])
        log10_predicted_occupancies_to_plot.append(numpy.mean(log10_predicted_occupancies_i))


    #bins_mean_ = numpy.asarray([numpy.mean(log10_predicted_occupancies[bin_edges_[i]: bin_edges_[i+1]]) for i in range(0, len(bin_edges_)-1 )])
    hist_to_plot = numpy.asarray(hist_to_plot)
    log10_predicted_occupancies_to_plot = numpy.asarray(log10_predicted_occupancies_to_plot)

    #ax_occupancy_taxon.scatter(10**hist_to_plot, 10**log10_predicted_occupancies_to_plot, c='k')

    bins_x_to_keep_no_nan, bins_y_no_nan = plot_utils.get_bin_mean_x_y(occupancy, predicted_occupancies, bins=40)

    #ax_occupancy_taxon.scatter(occupancy, predicted_occupancies, color=plot_utils.rgb_blue_taxon(rank_idx), s=40, linewidth=0.8, edgecolors='k')
    #ax_occupancy_taxon.scatter(bins_x_to_keep_no_nan, bins_y_no_nan, color=plot_utils.rgb_blue_taxon(rank_idx), s=40, alpha=0.9, linewidth=0.8, edgecolors='k')

    occupancy_all.extend(occupancy.tolist())
    predicted_occupancies_all.extend(predicted_occupancies.tolist())


occupancy_all = numpy.asarray(occupancy_all)
predicted_occupancies_all = numpy.asarray(predicted_occupancies_all)

predicted_and_observed_occupancies = numpy.concatenate((occupancy_all,predicted_occupancies_all),axis=0)
min_ = min(predicted_and_observed_occupancies)
max_ = max(predicted_and_observed_occupancies)

#ax_occupancy_taxon.set_xscale('log', base=10)
#ax_occupancy_taxon.set_yscale('log', base=10)
#ax_occupancy_taxon.set_xlabel("Observed occupancy", fontsize = 12)
#ax_occupancy_taxon.set_ylabel("Predicted occupancy", fontsize = 12)
#ax_occupancy_taxon.plot([min_,max_],[min_,max_], lw=4, ls=':',c='k', zorder=2, label='1:1')
#ax_occupancy_taxon.plot([min_*0.5,1.1],[min_*0.5,1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
#ax_occupancy_taxon.set_xlim([min_*0.5, 1.1])
#ax_occupancy_taxon.set_ylim([min_*0.5, 1.1])
#ax_occupancy_taxon.set_xlim([0.004, 1.1])
#ax_occupancy_taxon.set_ylim([0.004, 1.1])
#ax_occupancy_taxon.legend(loc="upper left", fontsize=7)




# occupancy phylo
occupancy_all = []
predicted_occupancies_all = []
for distance in distances:

    occupancy = occupancy_dict['phylo'][distance]['occupancy']
    predicted_occupancies = occupancy_dict['phylo'][distance]['predicted_occupancies']

    log10_occupancy = numpy.log10(occupancy)
    log10_predicted_occupancies = numpy.log10(predicted_occupancies)
    hist_, bin_edges_ = numpy.histogram(log10_occupancy, density=True, bins=20)

    hist_to_plot = []
    log10_predicted_occupancies_to_plot = []
    for i in range(0, len(bin_edges_)-1):

        log10_predicted_occupancies_i = log10_predicted_occupancies[(log10_occupancy >= bin_edges_[i]) & (log10_occupancy < bin_edges_[i+1])]
        if len(log10_predicted_occupancies_i) < 5:
            continue

        hist_to_plot.append(hist_[i])
        log10_predicted_occupancies_to_plot.append(numpy.mean(log10_predicted_occupancies_i))


    #bins_mean_ = numpy.asarray([numpy.mean(log10_predicted_occupancies[bin_edges_[i]: bin_edges_[i+1]]) for i in range(0, len(bin_edges_)-1 )])
    hist_to_plot = numpy.asarray(hist_to_plot)
    log10_predicted_occupancies_to_plot = numpy.asarray(log10_predicted_occupancies_to_plot)

    #ax_occupancy_taxon.scatter(10**hist_to_plot, 10**log10_predicted_occupancies_to_plot, c='k')
    #ax_occupancy_phylo.scatter(occupancy, predicted_occupancies, color=color_all(distances_all.index(distance)), s=35, alpha=0.2, linewidth=0.8, edgecolors='k')

    bins_x_to_keep_no_nan, bins_y_no_nan = plot_utils.get_bin_mean_x_y(occupancy, predicted_occupancies, bins=20)

    #ax_occupancy_phylo.scatter(occupancy, predicted_occupancies, color=color_all(distances.index(distance)), s=40, linewidth=0.8, edgecolors='k')
    ax_occupancy_phylo.scatter(bins_x_to_keep_no_nan, bins_y_no_nan, color=color_all(distances.index(distance)), s=35, alpha=0.9, linewidth=0.8, edgecolors='k')


    occupancy_all.extend(occupancy.tolist())
    predicted_occupancies_all.extend(predicted_occupancies.tolist())




occupancy_all = numpy.asarray(occupancy_all)
predicted_occupancies_all = numpy.asarray(predicted_occupancies_all)

predicted_and_observed_occupancies = numpy.concatenate((occupancy_all,predicted_occupancies_all),axis=0)
min_ = min(predicted_and_observed_occupancies)
max_ = max(predicted_and_observed_occupancies)

ax_occupancy_phylo.set_xscale('log', base=10)
ax_occupancy_phylo.set_yscale('log', base=10)
ax_occupancy_phylo.set_xlabel("Observed occupancy", fontsize = 12)
ax_occupancy_phylo.set_ylabel("Predicted occupancy", fontsize = 12)
#ax_occupancy_taxon.plot([min_,max_],[min_,max_], lw=4, ls=':',c='k', zorder=2, label='1:1')
ax_occupancy_phylo.plot([min_*0.5,1.1],[min_*0.5,1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax_occupancy_phylo.set_xlim([0.004, 1.1])
ax_occupancy_phylo.set_ylim([0.004, 1.1])
ax_occupancy_phylo.legend(loc="upper left", fontsize=7)



#ax_occupancy_phylo.set_xscale('log', base=10)
#ax_occupancy_phylo.set_yscale('log', base=10)
#ax_occupancy_phylo.set_xlabel("Observed occupancy", fontsize = 12)
#ax_occupancy_phylo.set_ylabel("Predicted occupancy", fontsize = 12)
#ax_occupancy_phylo.plot([0.009,1.1],[0.009,1.1], lw=5, ls=':',c='k', zorder=2, label='1:1')
#ax_occupancy_phylo.legend(loc="upper left", fontsize=7)



# MAD vs occupancy taxon
mad_all = []
predict_occupancy_all = []
for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

    mad = occupancy_dict['taxon'][rank]['mad']
    occupancy = occupancy_dict['taxon'][rank]['occupancy']
    predicted_occupancies = occupancy_dict['taxon'][rank]['predicted_occupancies']

    mad_all.extend(mad.tolist())
    predict_occupancy_all.extend(predicted_occupancies.tolist())

    bins_x_to_keep_no_nan, bins_y_no_nan = plot_utils.get_bin_mean_x_y(mad, occupancy)

    bins_x_to_keep_no_nan, bins_y_no_nan = plot_utils.get_bin_mean_x_y(mad, occupancy)

    #ax_mad_vs_occupancy_taxon.scatter(bins_x_to_keep_no_nan, bins_y_no_nan, color=plot_utils.rgb_blue_taxon(rank_idx), s=40, linewidth=0.8, alpha=0.9, edgecolors='k')
    

    #ax_mad_vs_occupancy_phylo.scatter(mad, occupancy, color=color_all(distances_all.index(distance)), s=35, alpha=0.2, linewidth=0.8, edgecolors='k')
    #ax_mad_vs_occupancy_phylo.scatter(bins_x_to_keep_no_nan, bins_y_no_nan, color=color_all(distances_all.index(distance)), s=35, alpha=0.9, linewidth=0.8, edgecolors='k')
 




mad_all = numpy.asarray(mad_all)
predict_occupancy_all = numpy.asarray(predict_occupancy_all)

#ax_mad_vs_occupancy_taxon.set_xscale('log', base=10)
#ax_mad_vs_occupancy_taxon.set_yscale('log', base=10)
#ax_mad_vs_occupancy_taxon.set_xlabel("Mean relative abundance", fontsize = 12)
#ax_mad_vs_occupancy_taxon.set_ylabel("Occupancy", fontsize = 12)


mad_log10_all = numpy.log10(mad_all)
predicted_occupancies_log10_all = numpy.log10(predict_occupancy_all)
hist_all, bin_edges_all = numpy.histogram(mad_log10_all, density=True, bins=25)
bins_mean_all = [0.5 * (bin_edges_all[i] + bin_edges_all[i+1]) for i in range(0, len(bin_edges_all)-1 )]
bins_mean_all_to_keep = []
bins_occupancies = []
for i in range(0, len(bin_edges_all)-1 ):
    predicted_occupancies_log10_i = predicted_occupancies_log10_all[(mad_log10_all>=bin_edges_all[i]) & (mad_log10_all<bin_edges_all[i+1])]

    if len(predicted_occupancies_log10_i) >= 3:
        #bins_mean_all_to_keep.append(bins_mean_all[i])
        bins_mean_all_to_keep.append(bin_edges_all[i])
        bins_occupancies.append(numpy.mean(predicted_occupancies_log10_i))


bins_mean_all_to_keep = numpy.asarray(bins_mean_all_to_keep)
bins_occupancies = numpy.asarray(bins_occupancies)

bins_mean_all_to_keep_no_nan = bins_mean_all_to_keep[(~numpy.isnan(bins_mean_all_to_keep)) & (~numpy.isnan(bins_occupancies))]
bins_occupancies_no_nan = bins_occupancies[(~numpy.isnan(bins_mean_all_to_keep)) & (~numpy.isnan(bins_occupancies))]

#ax_mad_vs_occupancy_taxon.plot(10**bins_mean_all_to_keep_no_nan, 10**bins_occupancies_no_nan, lw=4, ls='--',c='k', zorder=2, label='SLM prediction')
#ax_mad_vs_occupancy_taxon.legend(loc="upper left", fontsize=7)
#ax_mad_vs_occupancy_taxon.set_ylim([0.008, 1.08])



# mad vs occupancy, phylo
mad_all = []
predict_occupancy_all = []
for distance_idx, distance in enumerate(distances):

    if distance_idx %2 == 0:
        continue

    mad = occupancy_dict['phylo'][distance]['mad']
    occupancy = occupancy_dict['phylo'][distance]['occupancy']
    predicted_occupancies = occupancy_dict['phylo'][distance]['predicted_occupancies']

    mad_all.extend(mad.tolist())
    predict_occupancy_all.extend(predicted_occupancies.tolist())
    bins_x_to_keep_no_nan, bins_y_no_nan = plot_utils.get_bin_mean_x_y(mad, occupancy)

    #ax_mad_vs_occupancy_phylo.scatter(mad, occupancy, color=color_all(distances_all.index(distance)), s=35, alpha=0.2, linewidth=0.8, edgecolors='k')
    ax_mad_vs_occupancy_phylo.scatter(bins_x_to_keep_no_nan, bins_y_no_nan, color=color_all(distances.index(distance)), s=35, alpha=0.9, linewidth=0.8, edgecolors='k')
 


mad_all = numpy.asarray(mad_all)
predict_occupancy_all = numpy.asarray(predict_occupancy_all)

ax_mad_vs_occupancy_phylo.set_xscale('log', base=10)
ax_mad_vs_occupancy_phylo.set_yscale('log', base=10)
ax_mad_vs_occupancy_phylo.set_xlabel("Mean relative abundance", fontsize = 12)
ax_mad_vs_occupancy_phylo.set_ylabel("Occupancy", fontsize = 12)
ax_mad_vs_occupancy_phylo.set_ylim([0.008, 1.08])


mad_log10_all = numpy.log10(mad_all)
predicted_occupancies_log10_all = numpy.log10(predict_occupancy_all)
hist_all, bin_edges_all = numpy.histogram(mad_log10_all, density=True, bins=25)
#bins_mean_all = [0.5 * (bin_edges_all[i] + bin_edges_all[i+1]) for i in range(0, len(bin_edges_all)-1 )]
bins_mean_all_to_keep = []
bins_occupancies = []
for i in range(0, len(bin_edges_all)-1 ):
    predicted_occupancies_log10_i = predicted_occupancies_log10_all[(mad_log10_all>=bin_edges_all[i]) & (mad_log10_all<bin_edges_all[i+1])]

    if len(predicted_occupancies_log10_i) >= 3:
        #bins_mean_all_to_keep.append(bins_mean_all[i])
        bins_mean_all_to_keep.append(bin_edges_all[i])
        bins_occupancies.append(numpy.mean(predicted_occupancies_log10_i))


bins_mean_all_to_keep = numpy.asarray(bins_mean_all_to_keep)
bins_occupancies = numpy.asarray(bins_occupancies)

bins_mean_all_to_keep_no_nan = bins_mean_all_to_keep[(~numpy.isnan(bins_mean_all_to_keep)) & (~numpy.isnan(bins_occupancies))]
bins_occupancies_no_nan = bins_occupancies[(~numpy.isnan(bins_mean_all_to_keep)) & (~numpy.isnan(bins_occupancies))]

ax_mad_vs_occupancy_phylo.plot(10**bins_mean_all_to_keep_no_nan, 10**bins_occupancies_no_nan, lw=3, ls='--',c='k', zorder=2, label='Gamma prediction')
ax_mad_vs_occupancy_phylo.legend(loc="upper left", fontsize=7)




fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sgamma_summary.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
