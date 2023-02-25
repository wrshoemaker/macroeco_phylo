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
import simulation_utils
import mle_utils

import config



rarefied = True
iter = 1000

fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = list(range(len(taxa_ranks)))
idx_taxa_ranks.append(len(taxa_ranks))
taxa_ranks_label = diversity_utils.taxa_ranks_label
taxa_ranks_label.insert(0,'ASV')


color_all = numpy.asarray([diversity_utils.rgb_blue_taxon(r) for r in idx_taxa_ranks])
environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):
    #environment = 'freshwater metagenome'

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s, taxa)
        error = numpy.absolute(occupancies - predicted_occupancies)/occupancies

        idx_to_keep =  (~numpy.isinf(error)) & (error>0)
        species = species[idx_to_keep]
        error = error[idx_to_keep]

        error_dict = {}
        for s_idx, s in enumerate(species):
            error_dict[s] = {}
            error_dict[s]['otu'] = error[s_idx]

        for rank in diversity_utils.taxa_ranks:

            sys.stderr.write("Starting %s level analysis...\n" % rank)
            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))

            genus_to_taxa_dict = {}
            sad_genera_all = []
            g_taxa_idx_all = []
            for genus in all_genera:
                genus_to_taxa_dict[genus] = []
                for t in taxa:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
                g_taxa_idx_all.append(g_taxa_idx)
                sad_genera_all.append(s_by_s[g_taxa_idx,:].sum(axis=0))

            n_clades = [len(g_taxa_idx) for g_taxa_idx in g_taxa_idx_all]

            s_by_s_all_clades = numpy.stack(sad_genera_all, axis=0)
            # remove sites where there are no observations
            s_by_s_all_clades = s_by_s_all_clades[:,~(numpy.all(s_by_s_all_clades == 0, axis=0))]

            rel_s_by_s_all_clades =  s_by_s_all_clades/numpy.sum(s_by_s_all_clades, axis=0)
            shape = (numpy.mean(rel_s_by_s_all_clades, axis=1)**2)/numpy.var(rel_s_by_s_all_clades, axis=1)
            scale = shape/numpy.mean(rel_s_by_s_all_clades, axis=1)

            occupancies_clades, predicted_occupancies_clades, mad_clades, beta_clades, species_clades = diversity_utils.predict_occupancy(s_by_s_all_clades, all_genera)
            error_clades = numpy.absolute(occupancies_clades - predicted_occupancies_clades)/occupancies_clades

            #print(n_clades)
            #print(rank, shape, scale)

            idx_to_keep =  (~numpy.isinf(error_clades)) & (error_clades>0)
            error_clades = error_clades[idx_to_keep]
            species_clades = species_clades[idx_to_keep]

            for s_idx, s in enumerate(species_clades):
                for t_idx, t in enumerate(genus_to_taxa_dict[s]):

                    if t in error_dict:
                        error_dict[t][rank] = error_clades[s_idx]


        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        species = list(error_dict.keys())
        ranks_in_dict = ['otu', 'genus', 'family', 'order',  'class',  'phylum']
        for s in species:
            #error_to_plot = [error_dict[s][r] for r in ranks_in_dict if r in error_dict[s]]
            idx_taxa_ranks_to_plot = []
            error_to_plot = []
            for r_idx, r in enumerate(ranks_in_dict):

                if r in error_dict[s]:
                    idx_taxa_ranks_to_plot.append(r_idx)
                    error_to_plot.append(error_dict[s][r])

            if len(error_to_plot) < 3:
                continue

            ax.plot(idx_taxa_ranks_to_plot, error_to_plot, ls='-', lw=0.5, c='dodgerblue', alpha=0.005,  zorder=2)


        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)
        ax.set_yscale('log', base=10)

        ax.set_xticks(idx_taxa_ranks)
        ax.set_xticklabels(taxa_ranks_label, fontsize=8)

        #ax.set_ylim([-1.03,0.05])

        if environment_idx == 0:
            ax.set_ylabel("Occupancy prediction error, gamma", fontsize = 10)

        #if environment_chunk_idx == len(environment_chunk_all)-1:
        #    ax.set_xlabel("Rescaled log relative abundance", fontsize = 10)

        #if (environment_chunk_idx ==0) and (environment_idx == 0):
        #    ax.legend(loc="upper left", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%serror_coarse_emp_taxon%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
