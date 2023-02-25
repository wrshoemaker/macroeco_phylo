from __future__ import division
import os
#from Bio import Phylo
import random
import copy
import sys
import numpy
import random
import itertools
import scipy.stats as stats

import ete3
import diversity_utils
import config
import tree_utils

#import treeswift
import matplotlib.pyplot as plt

random.seed(123456789)


def get_c():
    total_c = 0
    for nodes in t.get_descendants():
        total_c +=  len(list(nodes.traverse("preorder")))



def get_c_and_a(t):

    #t = ete3.Tree('((((H,K)D,(F,I)G)B,E)A,((L,(N,Q)O)J,(P,S)M)C);', format=1)

    subtrees_dict = {}
    for nodes in t.get_descendants():

        # counts internal and tips
        n_nodes = len(list(nodes.traverse("preorder")))

        subtrees_dict[nodes] = {}
        subtrees_dict[nodes]['A'] = n_nodes
        subtrees_dict[nodes]['C'] = n_nodes

        for subtree in list(subtrees_dict.keys()):
            if nodes in subtree:
                subtrees_dict[subtree]['C'] += n_nodes

    subtrees = list(subtrees_dict.keys())

    A_all = [subtrees_dict[s]['A'] for s in subtrees]
    C_all = [subtrees_dict[s]['C'] for s in subtrees]

    return A_all, C_all


#get_c_and_a()


def test():

    samples_all = diversity_utils.all_observations()
    SADs_all, taxonomy_names_all, samples_keep_all = diversity_utils.get_SADs(samples_all)
    taxonomy_names_all_set = list(set(itertools.chain.from_iterable(taxonomy_names_all)))

    sys.stderr.write("Subsetting samples...\n")
    samples = diversity_utils.subset_observations()

    sys.stderr.write("Getting SADs...\n")
    SADs, taxonomy_names, samples_keep = diversity_utils.get_SADs(samples)

    richness_all = [len(s) for s in SADs]
    print(max(richness_all))


    #SAD_idx_list = random.sample(range(len(SADs)), 1)

    SAD_idx_list = [richness_all.index(max(richness_all))]

    SADs_sample = [SADs[i] for i in SAD_idx_list]
    taxonomy_names_sample = [taxonomy_names[i] for i in SAD_idx_list]
    samples_keep = [samples_keep[i] for i in SAD_idx_list]


    sys.stderr.write("Loading tree...\n")
    tree = ete3.Tree('%semp/otu_info/silva_123/97_otus.tre' % config.data_directory, quoted_node_names=True, format=1)
    # only keep species in human guts
    #tree = tree_utils.subset_tree(taxonomy_names_sample_i_null, tree)


    sys.stderr.write("Calculating measures of tree topology...\n")
    #record_strs = [",".join(['SampleID', 'evenness', 'symmetry'])]

    fig, ax = plt.subplots(figsize=(4,4))

    obs_a_all = []
    obs_c_all = []

    obs_a_null_all = []
    obs_c_null_all = []

    for i in range(len(SADs_sample)):

        SADs_sample_i = SADs_sample[i]
        taxonomy_names_sample_i = taxonomy_names_sample[i]
        samples_keep_i = samples_keep[i]

        tree_subset = tree_utils.subset_tree(taxonomy_names_sample_i, tree)

        a, c = get_c_and_a(tree_subset)

        obs_a_all.extend(a)
        obs_c_all.extend(c)

        taxonomy_names_sample_i_null = random.sample(taxonomy_names_all_set, len(taxonomy_names_sample_i))

        tree_subset_null = tree_utils.subset_tree(taxonomy_names_sample_i_null, tree)

        a_null, c_null = get_c_and_a(tree_subset_null)

        obs_a_null_all.extend(a_null)
        obs_c_null_all.extend(c_null)


    slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(obs_a_all), numpy.log10(obs_c_all))
    slope_null, intercept_null, r_value_null, p_value_null, std_err_null = stats.linregress(numpy.log10(obs_a_null_all), numpy.log10(obs_c_null_all))

    print(slope, slope_null)


    ax.scatter(obs_a_all, obs_c_all, alpha=0.2, s=4, c='dodgerblue', zorder=2)
    ax.scatter(obs_a_null_all, obs_c_null_all, alpha=0.2, s=4, c='orangered', zorder=2)

    # plot null
    a_range = numpy.logspace(min(numpy.log10(obs_a_all)), max(numpy.log10(obs_a_all)), num=400)
    c_min = 1 + (a_range+1)*((numpy.log(a_range + 1) / numpy.log(2))-1)
    c_max = (a_range**2)/4 + a_range - (1/4)

    ax.plot(a_range, c_min, 'k', ls='--', lw=1, zorder=1)
    ax.plot(a_range, c_max, 'k', ls='--', lw=1, zorder=1)

    ax.set_xlabel('Branch size', fontsize=12)
    ax.set_ylabel('Cumulative branch size', fontsize=12)

    ax.set_xscale('log', base=10)
    ax.set_yscale('log', base=10)

    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%stest_topology_max.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()



test()
