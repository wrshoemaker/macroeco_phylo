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

import diversity_utils
import config



rarefied = True
iter = 1000


#environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

environment = 'marine metagenome'

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = list(range(len(taxa_ranks)))
idx_taxa_ranks.append(len(taxa_ranks))
taxa_ranks_label = diversity_utils.taxa_ranks_label
taxa_ranks_label.insert(0,'OTU')


sys.stderr.write("Subsetting samples for %s...\n" % environment)
#samples = diversity_utils.subset_observations(environment=environment)
sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
samples = sad_annotated_dict['samples']

taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

sys.stderr.write("Running taxonomic coarse-graining.....\n")
afd_dict = {}
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

                if t not in afd_dict:
                    afd_dict[t] = {}

        g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])

        s_by_s_genus = s_by_s[g_taxa_idx,:]
        sad_genera_all.append(s_by_s_genus.sum(axis=0))
        rel_s_by_s_genus = rel_s_by_s[g_taxa_idx,:]

        for t_idx, t in enumerate(genus_to_taxa_dict[genus]):
            if 'otu' not in afd_dict[t]:
                afd_dict[t]['otu'] = rel_s_by_s_genus[t_idx,:]

                #if sum(rel_s_by_s_genus[t_idx,:]>0) > 180:
                #    print(t)
    s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
    # remove sites where there are no observations
    s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]
    rel_s_by_s_genera = (s_by_s_genera/s_by_s_genera.sum(axis=0))


    for genus_idx, genus in enumerate(all_genera):
        for t_idx, t in enumerate(genus_to_taxa_dict[genus]):
            afd_dict[t][rank] = {}
            afd_dict[t][rank]['name'] = genus
            afd_dict[t][rank]['afd'] = rel_s_by_s_genera[genus_idx,:]


focal_asv = 'DQ808483.1.1374'


fig = plt.figure(figsize = (8, 8)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

ax = plt.subplot2grid((2, 2), (0, 0))
ax_log = plt.subplot2grid((2, 2), (0, 1))

ax_genus = plt.subplot2grid((2, 2), (1, 0))
ax_log_genus = plt.subplot2grid((2, 2), (1, 1))


asv_all = list(afd_dict.keys())


finished_rank_dict = {}
for rank_idx, rank in enumerate(taxa_ranks):
    finished_rank_dict[rank] = []


for asv in asv_all:

    afd_asv = afd_dict[asv]['otu']
    afd_asv_no_zeros = afd_asv[afd_asv>0]

    if len(afd_asv_no_zeros) < 50:
        continue

    afd_asv_no_zeros_log = numpy.log10(afd_asv_no_zeros)

    rescaled_afd_asv_no_zeros = (afd_asv_no_zeros - numpy.mean(afd_asv_no_zeros))/numpy.std(afd_asv_no_zeros)
    rescaled_afd_asv_no_zeros_log = (afd_asv_no_zeros_log - numpy.mean(afd_asv_no_zeros_log))/numpy.std(afd_asv_no_zeros_log)

    hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(rescaled_afd_asv_no_zeros)
    log_hist_to_plot, log_bins_mean_to_plot = diversity_utils.get_hist_and_bins(rescaled_afd_asv_no_zeros_log)


    #ax.scatter(bins_mean_to_plot, hist_to_plot, s=10, color='k', alpha=0.9, lw=2, label=rank)
    #ax_log.scatter(log_bins_mean_to_plot, log_hist_to_plot, s=10, color='k', alpha=0.9, lw=2, label=rank)

    ax.plot(bins_mean_to_plot, hist_to_plot, color='k', alpha=0.05, lw=1)
    ax_log.plot(log_bins_mean_to_plot, log_hist_to_plot, color='k', alpha=0.05, lw=1)

    # run ranks
    for rank_idx, rank in enumerate(taxa_ranks):


        if rank != 'family':
            continue

        name_rank = afd_dict[asv][rank]['name']

        if name_rank in finished_rank_dict[rank]:
            continue


        afd_rank = afd_dict[asv][rank]['afd']

        afd_rank_no_zeros = afd_rank[afd_asv>0]
        afd_rank_no_zeros_log = numpy.log10(afd_rank_no_zeros)

        rescaled_afd_rank_no_zeros = (afd_rank_no_zeros - numpy.mean(afd_rank_no_zeros))/numpy.std(afd_rank_no_zeros)
        rescaled_afd_rank_no_zeros_log = (afd_rank_no_zeros_log - numpy.mean(afd_rank_no_zeros_log))/numpy.std(afd_rank_no_zeros_log)

        hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(rescaled_afd_rank_no_zeros)
        log_hist_to_plot, log_bins_mean_to_plot = diversity_utils.get_hist_and_bins(rescaled_afd_rank_no_zeros_log)

        ax_genus.plot(bins_mean_to_plot, hist_to_plot, color='k', alpha=0.05, lw=1)
        ax_log_genus.plot(log_bins_mean_to_plot, log_hist_to_plot, color='k', alpha=0.05, lw=1)


        finished_rank_dict[rank].append(afd_dict[asv][rank]['name'])






ax.set_yscale('log', base=10)
ax_log.set_yscale('log', base=10)
ax_genus.set_yscale('log', base=10)
ax_log_genus.set_yscale('log', base=10)



ax.set_xlabel("Rescaled AFD", fontsize = 10)
ax_genus.set_xlabel("Rescaled AFD", fontsize = 10)


ax_log.set_xlabel("Rescaled log AFD", fontsize = 10)
ax_log_genus.set_xlabel("Rescaled log AFD", fontsize = 10)


fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%ssingle_afd_coarse_taxon%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()

#for rank_idx, rank in enumerate(taxa_ranks):
