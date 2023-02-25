from __future__ import division
import config
import os
import sys
import subprocess
import random
import re
import itertools
import pickle
from collections import Counter

import numpy
import diversity_utils
import scipy.stats as stats

import ete3
import tree_utils
import random
import dbd_utils
import plot_utils

import matplotlib.pyplot as plt

rarefied = False








dbd_dict = dbd_utils.load_richness_dbd_dict()




fig = plt.figure(figsize = (8, 4)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

ax_taxon = plt.subplot2grid((1, 2), (0,0))
ax_phylo = plt.subplot2grid((1, 2), (0,1))

ax_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_taxon.transAxes)
ax_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_phylo.transAxes)

for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

distances = list(dbd_dict[environment]['phylo'].keys())
distances.sort()
taxa_ranks = list(dbd_dict[environment]['taxon'].keys())

color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)
color_environment_phylo = plot_utils.get_custom_cmap_phylo(environment, len(distances))

for coarse_rank_idx, coarse_rank in enumerate(taxa_ranks):
    
    observed = dbd_dict[environment]['taxon'][coarse_rank]['mean_slope']
    predicted = dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm']

    color_taxon = color_environment_taxon[coarse_rank_idx+1]

    ax_taxon.scatter(observed, predicted, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)

for d_idx, d in enumerate(distances):

    observed = dbd_dict[environment]['phylo'][d]['mean_slope']
    predicted = dbd_dict[environment]['phylo'][d]['mean_slope_slm']

    color_phylo = color_environment_phylo[d_idx]
    ax_phylo.scatter(observed, predicted, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)


ax_taxon.set_xscale('log', base=10)
ax_taxon.set_yscale('log', base=10)
#ax.plot([min_,max_], [min_,max_], lw=2,ls='--',c='k',zorder=1, label='1:1')
#ax.set_xlim([min_,max_])
#ax.set_ylim([min_,max_])
ax_taxon.plot([0.004,0.7], [0.004,0.7], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax_taxon.set_xlabel("Observed mean DBD slope", fontsize = 12)
ax_taxon.set_ylabel("Predicted mean DBD slope", fontsize = 12)
ax_taxon.legend(loc="upper left", fontsize=7)
ax_taxon.set_title("Taxonomic coarse-graining", fontsize=12)



ax_phylo.set_xscale('log', base=10)
ax_phylo.set_yscale('log', base=10)
ax_phylo.plot([0.004,0.7], [0.004,0.7], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax_phylo.set_xlabel("Observed mean DBD slope", fontsize = 12)
ax_phylo.set_ylabel("Predicted mean DBD slope", fontsize = 12)
ax_phylo.set_title("Phylogenetic coarse-graining", fontsize=12)
#ax_phylo.legend(loc="lower left", fontsize=7)




fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sdbd_slope_slm.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()




