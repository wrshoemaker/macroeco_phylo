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



rarefied = True
iter = 1000

fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = range(len(taxa_ranks))
taxa_ranks_label = diversity_utils.taxa_ranks_label

color_all = numpy.asarray([diversity_utils.rgb_blue_taxon(r) for r in idx_taxa_ranks])

environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #samples = diversity_utils.subset_observations(environment=environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))


        sys.stderr.write("Running taxonomic coarse-graining.....\n")
        evenness_vs_error_rho_all = []
        richness_vs_error_rho_all = []
        for rank_idx, rank in enumerate(taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)

            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))

            genus_to_taxa_dict = {}
            sad_genera_all = []
            evenness_fine_all = []
            evenness_fine_all_idx = []
            #sparsity_fine_all = []
            cv_beta_div_mean_all = []
            cv_beta_div_mean_all_idx = []
            richness_fine_all = []
            for genus_idx, genus in enumerate(all_genera):
                genus_to_taxa_dict[genus] = []
                for t in taxa:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])

                #rel_s_by_s_genus = rel_s_by_s[g_taxa_idx,:]
                s_by_s_genus = s_by_s[g_taxa_idx,:]
                evenness_fine = numpy.apply_along_axis(diversity_utils.calculate_pielou_evenness, 0, s_by_s_genus)
                sad_genera_all.append(s_by_s_genus.sum(axis=0))

                #sparsity_fine = numpy.apply_along_axis(diversity_utils.calculate_sparsity, 0, s_by_s_genus)
                #sparsity_fine_all.append(numpy.mean(sparsity_fine))

                richness_fine = numpy.apply_along_axis(diversity_utils.calculate_relative_richness, 0, s_by_s_genus)
                richness_fine_all.append(numpy.mean(richness_fine))
                #richness_fine_all_idx.append(True)

                #evenness_fine = evenness_fine[evenness_fine>0]
                rel_s_by_s_genus = rel_s_by_s[g_taxa_idx,:]

                mean_genus = numpy.mean(rel_s_by_s_genus, axis=1)
                std_genus = numpy.std(rel_s_by_s_genus, axis=1)
                beta_div_mean = mean_genus/(std_genus**2)

                if len(beta_div_mean)>5:

                    cv_beta_div_mean = numpy.std(beta_div_mean)/numpy.mean(beta_div_mean)
                    cv_beta_div_mean_all.append(cv_beta_div_mean)
                    cv_beta_div_mean_all_idx.append(True)

                else:
                    cv_beta_div_mean_all.append(0)
                    cv_beta_div_mean_all_idx.append(False)

                if sum(evenness_fine>0) > 5:
                    evenness_fine_all_idx.append(True)
                    evenness_fine_all.append(numpy.mean(evenness_fine))
                else:
                    evenness_fine_all_idx.append(False)
                    evenness_fine_all.append(0)


            s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
            # remove sites where there are no observations
            s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]
            rel_s_by_s_genera = (s_by_s_genera/s_by_s_genera.sum(axis=0))

            occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s_genera, list(range(s_by_s_genera.shape[0])))

            error = numpy.absolute(occupancies - predicted_occupancies)/occupancies
            error_to_keep_idx = (error>0)&(~numpy.isnan(error))

            #sparsity_fine_all = numpy.asarray(sparsity_fine_all)
            evenness_fine_all = numpy.asarray(evenness_fine_all)
            evenness_fine_all_idx = numpy.asarray(evenness_fine_all_idx)
            evenness_fine_all_to_keep = evenness_fine_all[(evenness_fine_all_idx) & (error_to_keep_idx)]
            error_to_keep = error[(evenness_fine_all_idx) & (error_to_keep_idx)]
            evenness_vs_error_rho = numpy.corrcoef(evenness_fine_all_to_keep, numpy.log10(error_to_keep))[1,0]
            evenness_vs_error_rho_all.append(evenness_vs_error_rho)


            # richness vs error
            richness_fine_all = numpy.asarray(richness_fine_all)
            richness_fine_all_to_keep = richness_fine_all[(error_to_keep_idx)]
            error_to_keep = error[(error_to_keep_idx)]
            richness_vs_error_rho = numpy.corrcoef(numpy.log10(richness_fine_all_to_keep), numpy.log10(error_to_keep))[1,0]
            richness_vs_error_rho_all.append(richness_vs_error_rho)


            #cv_beta_div_mean_all_idx = numpy.asarray(cv_beta_div_mean_all_idx)
            #error_to_keep = error[(cv_beta_div_mean_all_idx) & (error_to_keep_idx)]
            #cv_beta_div_mean_all = numpy.asarray(cv_beta_div_mean_all)
            #cv_beta_div_mean_all_to_keep = cv_beta_div_mean_all[(cv_beta_div_mean_all_idx) & (error_to_keep_idx)]


            #print(cv_beta_div_mean_all_to_keep)
            #cv_beta_div_mean_vs_error_rho = numpy.corrcoef(cv_beta_div_mean_all_to_keep, numpy.log10(error_to_keep))[1,0]
            #sparsity_vs_error_rho = numpy.corrcoef(sparsity_fine_all[error_to_keep_idx], numpy.log10(error[error_to_keep_idx]))[1,0]

            #ax.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=diversity_utils.rgb_blue_taxon(rank_idx), alpha=0.9, lw=2, label=rank)



        ax.plot(idx_taxa_ranks, evenness_vs_error_rho_all, ls='-', lw=1.5, c='k',  zorder=2)
        ax.scatter(idx_taxa_ranks, evenness_vs_error_rho_all, c=color_all, label='Observed',  zorder=3)


        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)
        #ax.set_yscale('log', base=10)

        ax.set_xticks(idx_taxa_ranks)
        ax.set_xticklabels(taxa_ranks_label, fontsize=9)

        ax.set_ylim([-1.03,0.05])

        if environment_idx == 0:
            ax.set_ylabel("Corr. b/w mean fine-grain evenness and\nmean log error of coarse-grained AFD", fontsize = 10)

        #if environment_chunk_idx == len(environment_chunk_all)-1:
        #    ax.set_xlabel("Rescaled log relative abundance", fontsize = 10)

        #if (environment_chunk_idx ==0) and (environment_idx == 0):
        #    ax.legend(loc="upper left", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sevenness_vs_error_taxon%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
