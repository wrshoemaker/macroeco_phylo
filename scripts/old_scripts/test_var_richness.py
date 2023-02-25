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
import simulation_utils



iter_ = 100

fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environments_to_keep = diversity_utils.environments_to_keep

environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = False)

        samples = sad_annotated_dict['samples']
        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        sys.stderr.write("Getting site-by-species matrix...\n")
        #s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)

        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

        prob_absence_all = []
        mean_prob_absence_gamma_all = []
        prob_absence_theory_all = []

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            print(rank)

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

            s_by_s_genera_pres_abs = (s_by_s_genera > 0)
            occupancy = numpy.sum(s_by_s_genera_pres_abs, axis=1)/s_by_s_genera_pres_abs.shape[1]
            prob_absence = 1-occupancy

            #prob_absence_theory = 
            #predi

            rel_s_by_s_genera = s_by_s_genera/numpy.sum(s_by_s_genera, axis=0)
            n_reads = numpy.sum(s_by_s_genera, axis=0)

            mean_rel_s_by_s = numpy.mean(rel_s_by_s_genera, axis=1)
            var_rel_s_by_s = numpy.var(rel_s_by_s_genera, axis=1)

            beta = (mean_rel_s_by_s**2)/var_rel_s_by_s
            theta = mean_rel_s_by_s/beta

            prob_absence_theory = []
            for i in range(len(theta)):
                prob_absence_theory.append(numpy.mean((1+theta[i]*n_reads)**(-1*beta[i])))
            prob_absence_theory = numpy.asarray(prob_absence_theory)


            prob_absence_gamma_all = []
            for i in range(iter_):
            
                rel_abundances, read_counts = simulation_utils.genrate_community_from_mean_and_var(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads))

                read_counts_pres_abs = (read_counts>0)
                occupancy_gamma = numpy.sum(read_counts_pres_abs, axis=1)/read_counts_pres_abs.shape[1]
                prob_absence_gamma = 1-occupancy_gamma

                prob_absence_gamma_all.append(prob_absence_gamma)

            mean_prob_absence_gamma = numpy.mean(prob_absence_gamma_all, axis=0)

            idx_to_keep = (prob_absence>0) & (mean_prob_absence_gamma>0) & (prob_absence_theory>0)
            prob_absence = prob_absence[idx_to_keep]
            mean_prob_absence_gamma = mean_prob_absence_gamma[idx_to_keep]
            prob_absence_theory = prob_absence_theory[idx_to_keep]

            prob_absence_all.extend(prob_absence.tolist())
            mean_prob_absence_gamma_all.extend(mean_prob_absence_gamma.tolist())
            prob_absence_theory_all.extend(prob_absence_theory.tolist())

            ax.scatter(prob_absence_theory, mean_prob_absence_gamma, c='k', edgecolors='k', alpha=0.3, s=10, zorder=3)  # , c='#87CEEB')


        min_ = min(mean_prob_absence_gamma_all + prob_absence_theory_all)
        max_ = max(mean_prob_absence_gamma_all + prob_absence_theory_all)

        ax.plot([min_*0.5,max_*1.1],[min_*0.5,max_*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
        ax.set_xlim([min_*0.5,max_*1.1])
        ax.set_ylim([min_*0.5,max_*1.1])

        #ax.set_xscale('log', base=10)
        #ax.set_yscale('log', base=10)

        ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)

        ax.set_xlabel("Predicted P(0|...), theory", fontsize = 10)
        ax.set_ylabel("Predicted P(0|...), simulation", fontsize = 10)






fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%stest_prob_absence_theory.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()

     
