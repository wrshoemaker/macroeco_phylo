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

#import process_dalbello_data




dbd_utils.make_diversity_vs_diversity_taxon_dict()
dbd_slope_dict = dbd_utils.load_diversity_slope_taxon_dalballo_dict()

n_resournces = list(dbd_slope_dict.keys())
n_resournces.sort()


fig = plt.figure(figsize = (5, 8)) #
fig.subplots_adjust(bottom= 0.15,  wspace=0.25)

ax_slope = plt.subplot2grid((2, 1), (0, 0), colspan=1)
ax_percentile = plt.subplot2grid((2, 1), (1, 0), colspan=1)

for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

    slope_median = []
    percentile_median = []

    for n in n_resournces:

        taxa = list(dbd_slope_dict[n]['ranks'][rank].keys())

        slope_all = [dbd_slope_dict[n]['ranks'][rank][t]['slope'] for t in taxa]
        percentile_all = [dbd_slope_dict[n]['ranks'][rank][t]['percentile'] for t in taxa]

        slope_median.append(numpy.median(slope_all))
        percentile_median.append(numpy.median(percentile_all))



    ax_slope.plot(n_resournces, slope_median, ls='-',  lw=1.5, alpha=0.9, c=diversity_utils.rgb_blue_taxon(rank_idx), label=diversity_utils.taxa_ranks_label[rank_idx])
    ax_percentile.plot(n_resournces, percentile_median, ls='-',  lw=1.5, alpha=0.9, c=diversity_utils.rgb_blue_taxon(rank_idx), label=diversity_utils.taxa_ranks_label[rank_idx])




ax_slope.set_xscale('log', base=2)
ax_percentile.set_xscale('log', base=2)

#ax_percentile.set_yscale('log', base=10)


ax_slope.set_xlabel('Number of resources', fontsize=12)
ax_slope.set_ylabel('Median DBD slope', fontsize=12)

ax_percentile.set_xlabel('Number of resources', fontsize=12)
ax_percentile.set_ylabel('Median position of observed DBD\nslope relative to the null distribution', fontsize=12)

ax_slope.legend(loc="upper left", fontsize=9)
ax_percentile.legend(loc="upper left", fontsize=9)

fig.subplots_adjust(wspace=0.3, hspace=0.3)
fig.savefig("%sdbd_slope_resource.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()








dbd_slope_phylo_dict = dbd_utils.load_diversity_slope_phylo_dalballo_dict()


n_resournces = list(dbd_slope_dict.keys())
n_resournces.sort()


fig = plt.figure(figsize = (5, 8)) #
fig.subplots_adjust(bottom= 0.15,  wspace=0.25)

ax_slope = plt.subplot2grid((2, 1), (0, 0), colspan=1)
ax_percentile = plt.subplot2grid((2, 1), (1, 0), colspan=1)


for n_idx, n in enumerate(n_resournces):

    distances = list(dbd_slope_phylo_dict[n].keys())
    distances.sort()

    slope_median = []
    percentile_median = []

    for d in distances:

        slope_median.append(numpy.mean(dbd_slope_phylo_dict[n][d]['slope']))
        percentile_median.append(numpy.mean(dbd_slope_phylo_dict[n][d]['percentile']))

    ax_slope.plot(distances, slope_median, ls='-',  lw=1.5, alpha=0.9, c=diversity_utils.rgb_blue_taxon(n_idx), label='%d resources' % n)
    ax_percentile.plot(distances, percentile_median, ls='-',  lw=1.5, alpha=0.9, c=diversity_utils.rgb_blue_taxon(n_idx), label='%d resources' % n)


ax_slope.set_xscale('log', base=10)
ax_percentile.set_xscale('log', base=10)

#ax_percentile.set_yscale('log', base=10)


ax_slope.set_xlabel('Phylogenetic distance', fontsize=12)
ax_slope.set_ylabel('Median DBD slope', fontsize=12)

ax_percentile.set_xlabel('Phylogenetic distance', fontsize=12)
ax_percentile.set_ylabel('Median position of observed DBD\nslope relative to the null distribution', fontsize=12)

ax_slope.legend(loc="upper left", fontsize=7)
ax_percentile.legend(loc="upper left", fontsize=7)




fig.subplots_adjust(wspace=0.3, hspace=0.3)
fig.savefig("%sdbd_slope_resource_phylo.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
