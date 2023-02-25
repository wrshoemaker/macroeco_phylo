
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



rarefied = False

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]


fig = plt.figure(figsize = (10, 20)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


for environment_idx, environment in enumerate(environments_to_keep):

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)
    samples = sad_annotated_dict['samples']

    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)

    sys.stderr.write("Running taxonomic coarse-graining.....\n")
    afd_dict = {}
    for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

        sys.stderr.write("Starting %s level analysis...\n" % rank)

        all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
        all_genera_idx = numpy.arange(len(all_genera))

        genus_to_taxa_dict = {}
        sad_genera_all = []

        var_asv_all = []
        var_genera_all = []
        for genus_idx, genus in enumerate(all_genera):
            genus_to_taxa_dict[genus] = []
            var_asv = 0
            for t in taxa:
                if sad_annotated_dict['taxa'][t][rank] == genus:
                    genus_to_taxa_dict[genus].append(t)
                    var_asv += var_rel_s_by_s[numpy.where(taxa==t)[0][0]]


            g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])

            s_by_s_genus = s_by_s[g_taxa_idx,:]
            sad_genera_all.append(s_by_s_genus.sum(axis=0))
            rel_s_by_s_genus = rel_s_by_s[g_taxa_idx,:]

            var_genus = numpy.var(rel_s_by_s_genus.sum(axis=0))

            var_asv_all.append(var_asv)
            var_genera_all.append(var_genus)


        ax = plt.subplot2grid((9, 5), (environment_idx, rank_idx))

        var_asv_all = numpy.asarray(var_asv_all)
        var_genera_all = numpy.asarray(var_genera_all)

        x, y, z = plot_utils.get_scatter_density_arrays_for_loglog(var_asv_all, var_genera_all)
        ax.scatter(x, y, c=numpy.sqrt(z), cmap='Blues', s=70, alpha=0.9, edgecolors='none', zorder=1)

        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)

        min_ = min(var_asv_all+var_genera_all)
        max_ = max(var_asv_all+var_genera_all)
        ax.plot([min_,max_],[min_,max_], lw=3,ls=':',c='k',zorder=1)

        if environment_idx == 0:
            ax.set_title(diversity_utils.taxa_ranks_label[rank_idx], fontsize=12)

        if rank_idx == 0:
            ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize = 12)


        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)


        #ax.set_xlabel("Sum of ASV variances", fontsize = 9)
        #ax.set_ylabel("Coarse-grained variance", fontsize = 9)


fig.text(0.3, 0.06, "Sum of OTU variances", va='center', fontsize=28)
fig.text(0.02, 0.5, "Coarse-grained variance", va='center', rotation='vertical', fontsize=28)


fig.subplots_adjust(hspace=0.37, wspace=0.37)
fig_name = "%sfine_vs_coarse_variance_taxon%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
