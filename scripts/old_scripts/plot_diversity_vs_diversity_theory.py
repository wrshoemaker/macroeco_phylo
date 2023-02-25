from __future__ import division
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

from scipy.special import polygamma




N = 10**6
s_tot = 10**4

mu = numpy.linspace(-6, -3, num=100)

#mu = -4
s = 0.5
sigma = 0.5
g = 1000

def predict_entropy(s):

    return -1 * (s_tot/N) * (1 - (sigma/2))* ((numpy.exp(mu + (s**2)/2) * (mu + s**2 - numpy.log(N))) +  numpy.log(sigma/2) + polygamma(0, 2/sigma) + mu)


def predict_coarse_entropy(s):
    return -1 * (s_tot/N) * (1 - (sigma/2))* ((numpy.exp(mu + (s**2)/2) * (mu + s**2 - numpy.log(N) +  numpy.log(g))) +  numpy.log(sigma/2) + polygamma(0, 2/sigma) + mu)



color_all = ['lightblue', 'dodgerblue', 'darkblue']
s_all = [0.5, 1, 1.5]


fig, ax = plt.subplots(figsize=(4,4))

H_all = []

for s_idx, s in enumerate(s_all):

    H = predict_entropy(s)
    H_coarse = predict_coarse_entropy(s)

    ax.plot(H, H_coarse, lw=1.5, ls='-', c=color_all[s_idx], label = 's = %s' % str(s))

    H_all.extend(H.tolist())

ax.plot([min(H_all),max(H_all)],[min(H_all),max(H_all)], lw=2,ls=':',c='k',zorder=1, label='1:1')

ax.legend(loc="upper left", fontsize=8)


ax.set_xlabel('Shannons diversity', fontsize=12)
ax.set_ylabel('Coarse-grained diversity', fontsize=12)

fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sfine_vs_coarse_theory.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()

print(H, H_coarse)
