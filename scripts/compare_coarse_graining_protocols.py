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
from matplotlib.lines import Line2D

import functools
import operator

import dbd_utils

import plot_utils

import diversity_utils
import config
import simulation_utils



rarefied = False

richness_diversity_prediction_phylo_dict_path = config.data_directory + "compare_coarse_graining_protocols_dict.pickle"

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]
environments_to_keep = ['human gut metagenome']



def make_dict():

    coarse_grained_tree_dict_all = {}
    distances_all = []
    for environment in environments_to_keep:
        sys.stderr.write("Loading tree dict for %s...\n" % environment)
        coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)
        coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict
        distances = list(coarse_grained_tree_dict.keys())
        distances_all.extend(distances)

    distances_all = list(set(distances_all))
    distances_all.sort()



    diversity_dict = {}

    #environment = 'human gut metagenome'

    for environment in environments_to_keep:

        coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]
        #coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)

        samples = pres_abs_dict['samples']
        taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
        sys.stderr.write("Getting site-by-species matrix...\n")

        s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])
        n_reads = numpy.sum(s_by_s, axis=0)
        rel_s_by_s = s_by_s/n_reads
        mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
        var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
        

        diversity_dict[environment] = {}

        distances_to_keep = []
        for distance_idx, distance in enumerate(distances_all):

            if distance in coarse_grained_tree_dict:
                distances_to_keep.append(distance)
                diversity_dict[environment][distance] = {}
                diversity_dict[environment][distance]['richness_simulate_then_coarse'] = []
                diversity_dict[environment][distance]['richness_coarse_then_simulate'] = []
                diversity_dict[environment][distance]['diversity_simulate_then_coarse'] = []
                diversity_dict[environment][distance]['diversity_coarse_then_simulate'] = []


        # simulate with observed matrix
        for distance_idx, distance in enumerate(distances_to_keep):
            sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
            coarse_grained_list = coarse_grained_tree_dict[distance]
            coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

            # get indexes for each clade
            coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

            # coarse grain s-by-s for all clades
            s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
            # remove zeros
            s_by_s_all_clades = s_by_s_all_clades[(numpy.sum(s_by_s_all_clades, axis=1)>0),:]

            mean_richness_observed, var_richness_observed, mean_richness_predicted, var_richness_predicted, mean_diversity_observed, var_diversity_observed, mean_diversity_predicted, var_diversity_predicted = diversity_utils.predict_mean_and_var_richness_and_diversity_using_gamma_rv(s_by_s_all_clades) 

            diversity_dict[environment][distance]['richness_observed'] = mean_richness_observed
            diversity_dict[environment][distance]['richness_simulated'] = mean_richness_predicted
            diversity_dict[environment][distance]['diversity_observed'] = mean_diversity_observed
            diversity_dict[environment][distance]['diversity_simulated'] = mean_diversity_predicted

            rel_error = numpy.absolute(mean_diversity_observed - mean_diversity_predicted)/mean_diversity_observed
            print("Observed", mean_diversity_observed, mean_diversity_predicted, rel_error)


        
        # simulate from simulated matrix.
        for i in range(100):

            # generate null matrix
            rel_abundances_gamma, read_counts_gamma = simulation_utils.genrate_community_from_mean_and_var_rvs(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads))

            for distance_idx, distance in enumerate(distances_to_keep):

                coarse_grained_list = coarse_grained_tree_dict[distance]
                coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

                # get indexes for each clade
                coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

                # coarse grain s-by-s for all clades
                s_by_s_all_clades = numpy.stack([numpy.sum(read_counts_gamma[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
                # remove zeros
                s_by_s_all_clades = s_by_s_all_clades[(numpy.sum(s_by_s_all_clades, axis=1)>0),:]

                mean_richness_observed_clades, var_richness_observed_clades, mean_richness_predicted_clades, var_richness_predicted_clades, mean_diversity_observed_clades, var_diversity_observed_clades, mean_diversity_predicted_clades, var_diversity_predicted_clades = diversity_utils.predict_mean_and_var_richness_and_diversity_using_gamma_rv(s_by_s_all_clades) 

                diversity_dict[environment][distance]['richness_simulate_then_coarse'].append(mean_richness_observed_clades)
                diversity_dict[environment][distance]['richness_coarse_then_simulate'].append(mean_richness_predicted_clades)
                diversity_dict[environment][distance]['diversity_simulate_then_coarse'].append(mean_diversity_observed_clades)
                diversity_dict[environment][distance]['diversity_coarse_then_simulate'].append(mean_diversity_predicted_clades)

                rel_error = numpy.absolute(mean_diversity_observed_clades - mean_diversity_predicted_clades)/mean_diversity_observed_clades
                print("Predicted", distance, mean_diversity_observed_clades, mean_diversity_predicted_clades, rel_error)



    sys.stderr.write("Saving prediction dictionary...\n")
    with open(richness_diversity_prediction_phylo_dict_path, 'wb') as handle:
        pickle.dump(diversity_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




def load_dict():

    with open(richness_diversity_prediction_phylo_dict_path, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_



diversity_dict = load_dict()



fig, ax = plt.subplots(figsize=(4,4))

environment = 'human gut metagenome'
distance_all = list(diversity_dict[environment].keys())
distance_all.sort()

for distance in distance_all:

    diversity_observed = diversity_dict[environment][distance]['diversity_observed']
    diversity_simulated = diversity_dict[environment][distance]['diversity_simulated']

    diversity_simulate_then_coarse = numpy.asarray(diversity_dict[environment][distance]['diversity_simulate_then_coarse'])
    diversity_coarse_then_simulate = numpy.asarray(diversity_dict[environment][distance]['diversity_coarse_then_simulate'])

    relative_error_observed = numpy.absolute(diversity_observed - diversity_simulated)/diversity_observed
    relative_error_simulation = numpy.mean(numpy.absolute(diversity_simulate_then_coarse - diversity_coarse_then_simulate)/diversity_simulate_then_coarse)

    ax.scatter(distance, relative_error_observed, c='blue')
    ax.scatter(distance, relative_error_simulation, c='k')

    print(relative_error_observed, relative_error_simulation)


ax.set_xscale('log', base=10)
ax.set_xlabel("Phylogenetic distance", fontsize = 12)
ax.set_ylabel("Diversity relative error", fontsize = 12)

legend_elements = [Line2D([0], [0], marker='o', color='w', label='Observed', markerfacecolor='blue', markersize=8),
                   Line2D([0], [0], marker='o', color='w', label='Gamma simulation', markerfacecolor='k', markersize=8)]

# Create the figure
ax.legend(handles=legend_elements, loc='upper left')



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%scompare_coarse_graining_protocols.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()