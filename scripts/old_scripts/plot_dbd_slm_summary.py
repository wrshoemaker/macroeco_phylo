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


#dbd_utils.make_dbd_slm_dict()
dbd_dict = dbd_utils.load_richness_dbd_dict()


fig = plt.figure(figsize = (8, 8)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

ax_taxon_dist = plt.subplot2grid((2, 2), (0,0))
ax_phylo_dist = plt.subplot2grid((2, 2), (0,1))

ax_taxon_predict = plt.subplot2grid((2, 2), (1,0))
ax_phylo_predict = plt.subplot2grid((2, 2), (1,1))

ax_taxon_dist.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_taxon_dist.transAxes)
ax_phylo_dist.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_phylo_dist.transAxes)
ax_taxon_predict.text(-0.1, 1.04, plot_utils.sub_plot_labels[2], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_taxon_predict.transAxes)
ax_phylo_predict.text(-0.1, 1.04, plot_utils.sub_plot_labels[3], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_phylo_predict.transAxes)


for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    distances = list(dbd_dict[environment]['phylo'].keys())
    distances.sort()
    taxa_ranks = list(dbd_dict[environment]['taxon'].keys())

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)
    color_environment_phylo = plot_utils.get_custom_cmap_phylo(environment, len(distances))

    slope_all_taxa = []
    for coarse_rank_idx, coarse_rank in enumerate(taxa_ranks):
        
        observed = dbd_dict[environment]['taxon'][coarse_rank]['mean_slope']
        #predicted = dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm_fix_beta']
        #predicted = dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm']
        predicted = dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm_fix_beta_across']


        

        color_taxon = color_environment_taxon[coarse_rank_idx+1]

        ax_taxon_predict.scatter(observed, predicted, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)

        # for distribution
        #slope_all_taxa.extend(dbd_dict[environment]['taxon'][coarse_rank]['slope_slm_all'])
        slope_all_taxa.extend(dbd_dict[environment]['taxon'][coarse_rank]['slope_slm_fix_beta_across_all'])

        
    hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(slope_all_taxa)
    ax_taxon_dist.scatter(bins_mean_to_plot, hist_to_plot, color=color_environment_taxon[-2], alpha=0.9, s=30, linewidth=0.8, edgecolors='k', label=diversity_utils.format_environment_label(environment))

   
    slope_all_phylo = []
    for d_idx, d in enumerate(distances):

        observed = dbd_dict[environment]['phylo'][d]['mean_slope']
        #predicted = dbd_dict[environment]['phylo'][d]['mean_slope_slm_fix_beta']
        #predicted = dbd_dict[environment]['phylo'][d]['mean_slope_slm']
        predicted = dbd_dict[environment]['phylo'][d]['mean_slope_slm_fix_beta_across']


        color_phylo = color_environment_phylo[d_idx]
        ax_phylo_predict.scatter(observed, predicted, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)

        # for distributions 
        #slope_all_phylo.extend(dbd_dict[environment]['phylo'][d]['slope_slm_fix_beta_all'])
        slope_all_phylo.extend(dbd_dict[environment]['phylo'][d]['slope_slm_fix_beta_across_all'])


    hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(slope_all_phylo)
    ax_phylo_dist.scatter(bins_mean_to_plot, hist_to_plot, color=color_environment_phylo[-3], alpha=0.9, s=30, linewidth=0.8, edgecolors='k', label=diversity_utils.format_environment_label(environment))



# plot distributions


ax_taxon_dist.set_xscale('log', base=10)
ax_taxon_dist.set_yscale('log', base=10)
ax_taxon_dist.set_xlabel("Predicted DBD slope", fontsize = 12)
ax_taxon_dist.set_ylabel("Probability density", fontsize = 12)
ax_taxon_dist.legend(loc="upper right", fontsize=6)

ax_phylo_dist.set_xscale('log', base=10)
ax_phylo_dist.set_yscale('log', base=10)
ax_phylo_dist.set_xlabel("Predicted DBD slope", fontsize = 12)
ax_phylo_dist.set_ylabel("Probability density", fontsize = 12)


ax_taxon_predict.set_xscale('log', base=10)
ax_taxon_predict.set_yscale('log', base=10)
#ax.plot([min_,max_], [min_,max_], lw=2,ls='--',c='k',zorder=1, label='1:1')
#ax.set_xlim([min_,max_])
#ax.set_ylim([min_,max_])
ax_taxon_predict.plot([0.004,0.7], [0.004,0.7], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax_taxon_predict.set_xlim([0.004,0.7])
ax_taxon_predict.set_ylim([0.004,0.7])
ax_taxon_predict.set_xlabel("Observed mean DBD slope", fontsize = 12)
ax_taxon_predict.set_ylabel("Predicted mean DBD slope", fontsize = 12)
ax_taxon_predict.legend(loc="upper left", fontsize=7)
ax_taxon_predict.set_title("Taxonomic coarse-graining", fontsize=12)



ax_phylo_predict.set_xscale('log', base=10)
ax_phylo_predict.set_yscale('log', base=10)
ax_phylo_predict.plot([0.0002,1.05], [0.0002,1.05], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax_phylo_predict.set_xlim([0.0002,1.05])
ax_phylo_predict.set_ylim([0.0002,1.05])
ax_phylo_predict.set_xlabel("Observed mean DBD slope", fontsize = 12)
ax_phylo_predict.set_ylabel("Predicted mean DBD slope", fontsize = 12)
ax_phylo_predict.set_title("Phylogenetic coarse-graining", fontsize=12)
#ax_phylo.legend(loc="lower left", fontsize=7)




fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sdbd_slope_slm.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()





