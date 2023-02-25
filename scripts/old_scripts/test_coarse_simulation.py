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
import simulation_utils



environment = 'human gut metagenome'
rarefied=True

sys.stderr.write("Subsetting samples for %s...\n" % environment)
#samples = diversity_utils.subset_observations(environment=environment)
sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
samples = sad_annotated_dict['samples']

taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))


coarse_grain_idx_dict = {}
for rank in diversity_utils.taxa_ranks:

    sys.stderr.write("Starting %s level analysis...\n" % rank)
    all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
    all_genera_idx = numpy.arange(len(all_genera))

    genus_to_taxa_dict = {}
    g_taxa_idx_all = []
    for genus in all_genera:
        genus_to_taxa_dict[genus] = []
        for t in taxa:
            if sad_annotated_dict['taxa'][t][rank] == genus:
                genus_to_taxa_dict[genus].append(t)

        g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
        g_taxa_idx_all.append(g_taxa_idx)
    counts_genera = [len(n) for n in g_taxa_idx_all]
    coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_genera))[:-1]
    coarse_grain_idx_dict[rank] = {}
    coarse_grain_idx_dict[rank]['counts'] = numpy.asarray(counts_genera)
    coarse_grain_idx_dict[rank]['coarse_grain_idx'] = coarse_grain_by_genus_idx



S, n_sites = s_by_s.shape
gm=0.4 #mean of sigma;
mu=-6 #mean of log(K)
s=3 #variance of log(K)
N=3*10**4 # number of reads

sim_s_by_s = simulation_utils.generate_community(mu, s, S, N, 'unif', gm, n_sites)
#rel_sim_s_by_s = sim_s_by_s/sim_s_by_s.sum(axis=0)
richness = numpy.sum(sim_s_by_s>0, axis=0)


diversity_fine = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, sim_s_by_s)

fig = plt.figure(figsize = (20, 3.5)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

    sim_s_by_s_coarse = numpy.add.reduceat(sim_s_by_s, coarse_grain_idx_dict[rank]['coarse_grain_idx'], axis=0)
    coarse_grain_counts = coarse_grain_idx_dict[rank]['counts']

    #richness_coarse = numpy.sum(sim_s_by_s_coarse>0, axis=0)

    predict_coarse_all = []
    for sad in sim_s_by_s_coarse.T:

        coarse_grain_counts_i = coarse_grain_counts[sad>0]
        predict_coarse = simulation_utils.predict_coarse_shannon_diversity(mu, s, gm, N, coarse_grain_counts_i)
        predict_coarse_all.append(predict_coarse)



    diversity_coarse = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, sim_s_by_s_coarse)

    idx_to_keep = (diversity_fine>0) & (diversity_coarse>0) & (numpy.isfinite(diversity_fine)) & (numpy.isfinite(diversity_coarse))
    diversity_fine_to_plot = diversity_fine[idx_to_keep]
    diversity_coarse_to_plot = diversity_coarse[idx_to_keep]

    ax = plt.subplot2grid((1, 5), (0, rank_idx))

    ax.scatter(diversity_fine, diversity_coarse, s=10, color='dodgerblue', alpha=0.3, lw=2, label=rank)

    predict_fine = simulation_utils.predict_shannon_diversity(mu, s, gm, richness, N)
    #predict_coarse = predict_shannon_diversity(mu, s, sigma, richness, N)

    ax.scatter(predict_fine, predict_coarse_all, s=7, color='k', alpha=1, lw=3)

    ax.set_title(rank, fontsize=12)
    #ax.set_xscale('log', base=10)
    #ax.set_yscale('log', base=10)

    min_ = min(diversity_fine)
    max_ = max(diversity_fine)
    ax.plot([min_,max_],[min_,max_], lw=3,ls=':',c='k',zorder=1)

    ax.set_xlabel("Fine-grained diversity", fontsize = 12)
    ax.set_ylabel("Coarse-grained diversity", fontsize = 12)

    ax.set_title(rank, fontsize=12)





fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%stest_sim_diversity_vs_diversity.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
