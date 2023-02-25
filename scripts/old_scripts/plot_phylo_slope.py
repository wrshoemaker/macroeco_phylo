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


environment='human gut metagenome'

diversity_slope_phylo_dict = dbd_utils.load_diversity_slope_phylo_dict(environment)


distacnces = list(diversity_slope_phylo_dict.keys())
distacnces.sort()

mean_slope_to_plot = []
for d in distacnces:
    print(d, len(diversity_slope_phylo_dict[d]['slope']))
    mean_slope_to_plot.append(numpy.mean(diversity_slope_phylo_dict[d]['slope']))



fig, ax = plt.subplots(figsize=(4,4))

ax.plot(distacnces, mean_slope_to_plot)

ax.set_xscale('log', base=10)
ax.set_xlabel("Phylogenetic distance", fontsize = 12)
ax.set_ylabel("Mean slope", fontsize = 11)


fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%stest_phylo_slope.png" % config.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
