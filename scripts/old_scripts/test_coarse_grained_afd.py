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


environment = 'human gut metagenome'
rarefied = True


sys.stderr.write("Loading tree dict for %s...\n" % environment)
coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_dict(environment=environment, rarefied=rarefied)

distances = list(coarse_grained_tree_dict.keys())
distances.sort()

#distances = [0.01, 0.1274274985703134]


sys.stderr.write("Subsetting samples...\n")
#samples = diversity_utils.subset_observations(environment=environment)
sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
samples = sad_annotated_dict['samples']

sys.stderr.write("Getting site-by-species matrix...\n")
s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)


#print(numpy.sum(s_by_s > 0, axis=1))

#print(sum(numpy.all(s_by_s > 0, axis=1)))


fig = plt.figure(figsize = (8, 8)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

ax_taxon = plt.subplot2grid((2, 2), (0, 0))
ax_phylo = plt.subplot2grid((2, 2), (0, 1))

ax_taylor_taxon = plt.subplot2grid((2, 2), (1, 0))
ax_taylor_phylo = plt.subplot2grid((2, 2), (1, 1))


mean_genera_all = []
var_genera_all = []
sys.stderr.write("Running taxonomic coarse-graining.....\n")
for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

    sys.stderr.write("Starting %s level analysis...\n" % rank)

    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
    all_genera_idx = numpy.arange(len(all_genera))

    s_by_s_taxon = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
    rel_s_by_s_taxon = (s_by_s_taxon/s_by_s_taxon.sum(axis=0))

    genus_to_taxa_dict = {}
    sad_genera_all = []

    mean_ratio_all = []
    beta_div_mean_ratio_all = []
    for genus in all_genera:
        genus_to_taxa_dict[genus] = []
        for t in taxa:
            if sad_annotated_dict['taxa'][t][rank] == genus:
                genus_to_taxa_dict[genus].append(t)

        g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
        sad_genera_all.append(s_by_s_taxon[g_taxa_idx,:].sum(axis=0))

        rel_s_by_s_taxon_genus = rel_s_by_s_taxon[g_taxa_idx,:]

        mean_coarse_rel_s_by_s_taxon_genus = numpy.mean(numpy.sum(rel_s_by_s_taxon_genus, axis=0))
        std_coarse_rel_s_by_s_taxon_genus = numpy.std(numpy.sum(rel_s_by_s_taxon_genus, axis=0))

        beta_div_mean_coarse = (mean_coarse_rel_s_by_s_taxon_genus)/(std_coarse_rel_s_by_s_taxon_genus**2)

        mean_genus_all = numpy.mean(rel_s_by_s_taxon_genus, axis=1)
        std_genus_all = numpy.std(rel_s_by_s_taxon_genus, axis=1)
        #print(rel_s_by_s_taxon[g_taxa_idx,:])
        #print(len(mean_genus_all))
        beta_div_mean = (mean_genus_all)/(std_genus_all**2)

        #if len(mean_genus_all) < 2:
        #    continue


        #mode_beta_div_mean = stats.mode(beta_div_mean)[0][0]
        #median_beta_div_mean = numpy.median(beta_div_mean)


        mean_ratio_all.append(max(mean_genus_all)/mean_coarse_rel_s_by_s_taxon_genus)
        #beta_div_mean_ratio_all.append(median_beta_div_mean/beta_div_mean_coarse)

        #print(beta)
        #print(len(std_genus_all))
        #if len(beta) > 4:
        #    cv_beta = numpy.std(beta)/numpy.mean(beta)

        #    print(len(beta), cv_beta)

    s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
    # remove sites where there are no observations
    s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]
    rel_s_by_s_genera = (s_by_s_genera/s_by_s_genera.sum(axis=0))

    clade_log10_rescaled_all = []
    for clade in rel_s_by_s_genera:

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

    ax_taxon.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=diversity_utils.rgb_blue_taxon(rank_idx), alpha=0.9, lw=2, label=rank)

    #mean_genera = numpy.mean(rel_s_by_s_genera, axis=1)
    #var_genera = numpy.var(rel_s_by_s_genera, axis=1)

    #mean_genera_all.append(mean_genera)
    #var_genera_all.append(var_genera)

    #print(mean_ratio_all)

    #print(len(mean_ratio_all))

    mean_ratio_log10_all = numpy.log10(mean_ratio_all)

    hist_, bin_edges_ = numpy.histogram(mean_ratio_log10_all, density=True, bins=10)
    bins_mean_ = numpy.asarray([0.5 * (bin_edges_[i] + bin_edges_[i+1]) for i in range(0, len(bin_edges_)-1 )])
    hist_to_plot = hist_[hist_>0]
    bins_mean_to_plot = bins_mean_[hist_>0]

    #ax_taylor_taxon.scatter(10**bins_mean_to_plot, hist_to_plot, s=10, color=diversity_utils.rgb_blue_taxon(rank_idx), alpha=0.9, lw=2, label=rank)

    ax_taylor_taxon.plot(10**bins_mean_to_plot, hist_to_plot, lw=2, linestyle='-', c=diversity_utils.rgb_blue_taxon(rank_idx), alpha=0.9)

    #ax_taylor_taxon.hist(numpy.log10(mean_ratio_all), 8, lw=2, histtype='step', density=True, color=diversity_utils.rgb_blue_taxon(rank_idx), alpha=0.8)


    #ax_taylor_taxon.scatter(mean_genera, var_genera,  s=6, color=diversity_utils.rgb_blue_taxon(rank_idx), alpha=0.6)
    #slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(mean_genera), numpy.log10(var_genera))




#mean_genera_all = numpy.concatenate(mean_genera_all).ravel()
#var_genera_all = numpy.concatenate(var_genera_all).ravel()

#slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(mean_genera_all), numpy.log10(var_genera_all))

#x_log10_range = numpy.linspace(min(numpy.log10(mean_genera_all)) , max(numpy.log10(mean_genera_all)) , 10000)
#y_log10_null_range = 10 ** (2*x_log10_range + intercept)
#ax_taylor_taxon.plot(10**x_log10_range, y_log10_null_range, c='k', lw=2.5, linestyle='-', zorder=2, label= r'$y \sim x^{2}$')





sys.stderr.write("Running phylogenetic coarse-graining.....\n")
#fig, ax = plt.subplots(figsize=(4,4))
color_all = diversity_utils.make_blue_cmap(len(distances))
mean_genera_all = []
var_genera_all = []
for distance_idx, distance in enumerate(distances):

    sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
    coarse_grained_list = coarse_grained_tree_dict[distance]
    coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

    #if len(coarse_grained_n) == len(taxonomy_names):
    #    continue

    # get indexes for each clade
    coarse_grained_idx_all = [numpy.asarray([numpy.where(taxonomy_names==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
    # coarse grain s-by-s for all clades
    s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
    rel_s_by_s_all_clades = (s_by_s_all_clades/s_by_s_all_clades.sum(axis=0))

    mean_genus_all = []
    max_mean_fine_all = [max(numpy.mean(s_by_s[coarse_grained_idx,:], axis=1)) for coarse_grained_idx in coarse_grained_idx_all]
    mean_s_by_s_all_clades = numpy.mean(s_by_s_all_clades, axis=1)

    #print(len(mean_mean_fine_all), mean_s_by_s_all_clades.shape)

    #if len(mean_genus_all) < 5:
    #    continue

    #mean_ratio_all.append(max(mean_genus_all)/mean_coarse_rel_s_by_s_taxon_genus)


    clade_log10_rescaled_all = []
    for clade in rel_s_by_s_all_clades:

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
    #print(len(hist_to_plot), len(hist_))
    bins_mean_to_plot = bins_mean_[hist_>0]

    #ax.plot(bins_mean_to_plot, hist_to_plot, ls='-', c=color_all(distance_idx), alpha=0.3, lw=2)
    ax_phylo.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=color_all(distance_idx), alpha=0.9, lw=2)


    mean_genera = numpy.mean(rel_s_by_s_all_clades, axis=1)
    var_genera = numpy.var(rel_s_by_s_all_clades, axis=1)

    mean_genera_all.append(mean_genera)
    var_genera_all.append(var_genera)

    slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(mean_genera), numpy.log10(var_genera))


    # plot mean ratio

    mean_ratio_all = max_mean_fine_all/mean_s_by_s_all_clades

    mean_ratio_log10_all = numpy.log10(mean_ratio_all)

    hist_, bin_edges_ = numpy.histogram(mean_ratio_log10_all, density=True, bins=10)
    bins_mean_ = numpy.asarray([0.5 * (bin_edges_[i] + bin_edges_[i+1]) for i in range(0, len(bin_edges_)-1 )])
    hist_to_plot = hist_[hist_>0]
    bins_mean_to_plot = bins_mean_[hist_>0]

    #ax_taylor_taxon.scatter(10**bins_mean_to_plot, hist_to_plot, s=10, color=diversity_utils.rgb_blue_taxon(rank_idx), alpha=0.9, lw=2, label=rank)

    #ax_taylor_phylo.scatter(10**bins_mean_to_plot, hist_to_plot, s=10, color=diversity_utils.rgb_blue_taxon(distance_idx), alpha=0.9)

    ax_taylor_phylo.plot(10**bins_mean_to_plot, hist_to_plot, lw=2, linestyle='-', c=color_all(distance_idx), alpha=0.9)



ax_taxon.set_yscale('log', base=10)
ax_phylo.set_yscale('log', base=10)

ax_taxon.legend(loc="upper left", fontsize=7)
ax_taxon.set_xlabel("Rescaled log relative abundance", fontsize = 12)
ax_taxon.set_ylabel("Probability density", fontsize = 11)


ax_phylo.set_xlabel("Rescaled log relative abundance", fontsize = 12)
ax_phylo.set_ylabel("Probability density", fontsize = 11)


#norm = mpl.colors.LogNorm(vmin=0,vmax=2)

norm = mpl.colors.LogNorm(vmin=min(distances), vmax=max(distances))

pcm = plt.cm.ScalarMappable(cmap=color_all, norm=norm)
clb = plt.colorbar(pcm, ax=ax_phylo)
clb.set_label(label='Phylogenetic distance')


ax_taylor_taxon.set_xscale('log', base=10)
ax_taylor_taxon.set_yscale('log', base=10)


ax_taylor_phylo.set_xscale('log', base=10)
ax_taylor_phylo.set_yscale('log', base=10)


ax_taylor_taxon.set_xlabel("Ratio of max. mean fine-grained and\nmean coarse-grained rel. abund." , fontsize = 12)
ax_taylor_taxon.set_ylabel("Probability density", fontsize = 12)

ax_taylor_phylo.set_xlabel("Ratio of max. mean fine-grained and\nmean coarse-grained rel. abund." , fontsize = 12)
ax_taylor_phylo.set_ylabel("Probability density", fontsize = 12)


#ax_taylor_phylo.set_xlabel("Mean relative abundance", fontsize = 12)
#ax_taylor_phylo.set_ylabel("Variance of relative abundance", fontsize = 12)


ax_taylor_taxon.set_xlim([0.08,1])
ax_taylor_phylo.set_xlim([0.05,1])


ax_taylor_taxon.xaxis.set_tick_params(labelsize=6)
ax_taylor_phylo.xaxis.set_tick_params(labelsize=6)





fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%stest_coarse_grained_afd%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
