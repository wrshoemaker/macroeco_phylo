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




environment =  'marine metagenome'

sys.stderr.write("Subsetting samples for %s...\n" % environment)
sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = False)

samples = sad_annotated_dict['samples']
taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
sys.stderr.write("Getting site-by-species matrix...\n")

s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

taxa_ranks = ['family']
#diversity_utils.taxa_ranks

fig, ax = plt.subplots(figsize=(4,4))

for rank_idx, rank in enumerate(taxa_ranks):

    if rank == 'OTU':

        s_by_s_genera = s_by_s

    else:

        sys.stderr.write("Running taxonomic coarse-graining.....\n")

        sys.stderr.write("Starting %s level analysis...\n" % rank)
        all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
        all_genera_idx = numpy.arange(len(all_genera))

        genus_to_taxa_dict = {}
        sad_genera_all = []
        for genus_idx, genus in enumerate(all_genera):
            genus_to_taxa_dict[genus] = []
            var_asv = 0
            for t in taxa:
                if sad_annotated_dict['taxa'][t][rank] == genus:
                    genus_to_taxa_dict[genus].append(t)


            g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
            s_by_s_genus = s_by_s[g_taxa_idx,:]
            sad_genera_all.append(s_by_s_genus.sum(axis=0))


        s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
        # remove sites where there are no observations
        s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]
    

    expected_diversity_second_moment_first_term, expected_diversity_second_moment_first_term_rvs, expected_diversity_second_moment_second_term, expected_diversity_second_moment_second_term_rvs = diversity_utils.predict_second_moment(s_by_s_genera)

    # test variance
    idx_to_remove = (expected_diversity_second_moment_second_term > 10**-10.5)

    sum_second_rvs = sum(expected_diversity_second_moment_second_term[idx_to_remove])
    sum_second = sum(expected_diversity_second_moment_second_term_rvs[idx_to_remove])

    #print(sum_second_rvs, sum_second)


    expected_diversity_second_moment = expected_diversity_second_moment_second_term
    expected_diversity_second_moment_rvs = expected_diversity_second_moment_second_term_rvs

    idx_to_keep = (expected_diversity_second_moment>0) & (expected_diversity_second_moment_rvs>0)
    expected_diversity_second_moment = expected_diversity_second_moment[idx_to_keep]
    expected_diversity_second_moment_rvs = expected_diversity_second_moment_rvs[idx_to_keep]

    min_ = min(numpy.union1d(expected_diversity_second_moment, expected_diversity_second_moment_rvs))
    max_ = max(numpy.union1d(expected_diversity_second_moment, expected_diversity_second_moment_rvs))

    ax.plot([min_*0.5,max_*1.1],[min_*0.5,max_*1.1], lw=2,ls='--',c='k',zorder=3, label='1:1')
    ax.set_xlim([min_*0.5,max_*1.1])
    ax.set_ylim([min_*0.5,max_*1.1])
    ax.set_xlabel("Predicted second moment " +r'$ \left< ( x_{i} \mathrm{ln}[x_{i}])^{2}  \right>$'  + ", analytic", fontsize = 10)
    ax.set_ylabel("Predicted second momen "  +r'$ \left< ( x_{i} \mathrm{ln}[x_{i}])^{2}  \right>$' ", simulation", fontsize = 10)

    ax.scatter(expected_diversity_second_moment, expected_diversity_second_moment_rvs, s=7, color='dodgerblue', alpha=0.08, zorder=2)

    ax.set_xscale('log', base=10)
    ax.set_yscale('log', base=10)

    ax.legend(loc='upper left', fontsize=8, frameon=False)


    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%ssecond_moment_diversity.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()



