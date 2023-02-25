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
import mle_utils



rarefied = True
iter = 1000

#fig = plt.figure(figsize = (10, 10)) #
#fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

#environments_to_keep = diversity_utils.environments_to_keep
environment = 'human gut metagenome'

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = range(len(taxa_ranks))
taxa_ranks_label = diversity_utils.taxa_ranks_label


sys.stderr.write("Subsetting samples for %s...\n" % environment)
#samples = diversity_utils.subset_observations(environment=environment)
sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
samples = sad_annotated_dict['samples']

taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))


#print(taxa_ranks)

s_all = []
sum_div_max_all = []

sys.stderr.write("Running taxonomic coarse-graining.....\n")
for rank_idx, rank in enumerate(taxa_ranks):

    sys.stderr.write("Starting %s level analysis...\n" % rank)

    all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
    all_genera_idx = numpy.arange(len(all_genera))

    genus_to_taxa_dict = {}
    sad_genera_all = []
    for genus_idx, genus in enumerate(all_genera):
        genus_to_taxa_dict[genus] = []
        for t in taxa:
            if sad_annotated_dict['taxa'][t][rank] == genus:
                genus_to_taxa_dict[genus].append(t)

        g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])

        rel_s_by_s_genus = rel_s_by_s[g_taxa_idx,:]
        sum_rel_s_by_s_genus = rel_s_by_s_genus.sum(axis=0)
        max_rel_s_by_s_genus = numpy.amax(rel_s_by_s_genus, axis=0)

        idx_to_keep = max_rel_s_by_s_genus>0

        sum_div_max = max_rel_s_by_s_genus[idx_to_keep]/sum_rel_s_by_s_genus[idx_to_keep]
        sum_div_max_all.append(numpy.mean(sum_div_max))
        #sum_div_max_all.extend(sum_div_max.tolist())
        s_all.append(rel_s_by_s_genus.shape[0])
        #s_all.extend([rel_s_by_s_genus.shape[0]] * len(sum_div_max))
        #print(len(sum_rel_s_by_s_genus))





fig, ax = plt.subplots(figsize=(4,4))


ax.scatter(s_all, sum_div_max_all, s=7, color='k', alpha=0.5)


ax.set_xscale('log', base=10)
#ax.set_yscale('log', base=10)


fig.subplots_adjust(wspace=0.3, hspace=0.3)
fig.savefig("%safd_sum_vs_max.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
plt.close()
