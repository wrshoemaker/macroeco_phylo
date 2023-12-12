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
import plot_utils

import diversity_utils
import config

from scipy.special import rel_entr
import plot_utils
from collections import Counter
import tree_utils


distance_correlation_dict = config.data_directory + "distance_correlation_dict/distance_correlation_dict_%s.pickle"


rarefied = False

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

def make_rho_dict():

    tree = tree_utils.get_emp_tree()

    for environment_idx, environment in enumerate(environments_to_keep):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]
        pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)

        samples = pres_abs_dict['samples']
        taxa = list(pres_abs_dict['taxa'].keys())
        taxa_array = numpy.asarray(taxa)

        sys.stderr.write("Getting site-by-species matrix...\n")
        s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

        rho = numpy.corrcoef(rel_s_by_s)

        idx_to_keep = numpy.random.choice(rho.shape[0], size=500, replace=False)
        #test = rho[idx_to_keep, idx_to_keep]
        # dont know why you cant index both axes of a matrix simultaneously...
        rho = rho[idx_to_keep, :]
        rho = rho[:, idx_to_keep]
        taxa_array = taxa_array[idx_to_keep]
        taxa = taxa_array.tolist()
        tree_environment = tree_utils.subset_tree(taxa, tree)

        distance_all = []
        rho_all = []
        #for t_i in range(len(idx_to_keep)):
        for t_i in range(len(taxa)):

            for t_j in range(t_i):

                distance_ij = tree_environment.get_distance(taxa[t_i], taxa[t_j])
                rho_ij = rho[t_i, t_j]

                distance_all.append(distance_ij)
                rho_all.append(rho_ij)


        rho_dict = {}
        #rho_dict[environment] = {}
        rho_dict['distance'] = distance_all
        rho_dict['rho'] = rho_all

        dict_path = distance_correlation_dict % diversity_utils.get_environment_label(environment)

        sys.stderr.write("Saving diversity dictionary...\n")
        with open(dict_path, 'wb') as handle:
            pickle.dump(rho_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


#make_rho_dict()


fig = plt.figure(figsize = (20, 20)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


count_ = 0
environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):
    for environment_idx, environment in enumerate(environment_chunk):

        print(environment)

        #if count_ > 0:
        #    continue

        dict_path = distance_correlation_dict % diversity_utils.get_environment_label(environment)

        with open(dict_path, 'rb') as handle:
            rho_dict = pickle.load(handle)

        distance_all = numpy.asarray(rho_dict['distance'])
        rho_all = numpy.asarray(rho_dict['rho'])

        #fig, ax = plt.subplots(figsize=(4,4))
        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        x, y, z = plot_utils.get_scatter_density_arrays_for_linearlog(distance_all, rho_all, color_radius=0.5)
        ax.scatter(10**x, y, c=numpy.sqrt(z), cmap='Blues', s=30, alpha=0.9, edgecolors='none', zorder=1)
        #ax.scatter(distance_all, rho_all, c='k', s=70, alpha=0.9, edgecolors='none', zorder=1)

        ax.set_xscale('log', base=10)

        ax.set_xlabel("Phylogenetic distance", fontsize = 12)
        ax.set_ylabel("Correlation", fontsize = 12)

        ax.set_title(' '.join(environments_to_keep[0].split(' ')[:-1]).capitalize(), fontsize=11)

        count_ += 1



fig.subplots_adjust(hspace=0.37, wspace=0.37)
fig_name = "%sdistance_vs_rho.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()

#environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
#ax_all = []
#for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

#    for environment_idx, environment in enumerate(environment_chunk):