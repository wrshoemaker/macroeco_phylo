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



rarefied = False



#dbd_utils.make_richness_diversity_prediction_taxon_dict()

dbd_utils.make_richness_diversity_prediction_phylo_dict()

#dict_taxon = dbd_utils.load_richness_diversity_prediction_taxon_dict()

#print(dict_taxon)

fig, ax = plt.subplots(figsize=(4,4))

variance_richness_observed_all = []
variance_richness_predicted_all = []
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    continue


    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)

    samples = sad_annotated_dict['samples']
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    sys.stderr.write("Getting site-by-species matrix...\n")
    #s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)

    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

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


        #var_richness_observed, var_richness_predicted = diversity_utils.predict_var_richness(s_by_s_genera, s_by_s_genera.shape[0])

        #mean_richness_observed, mean_richness_predicted = diversity_utils.predict_mean_richness(s_by_s_genera, s_by_s_genera.shape[0])
        mean_richness_observed_gamma_rvs, var_richness_observed_gamma_rvs, mean_richness_predicted_gamma_rvs, var_richness_predicted_gamma_rvs, mean_richness_predicted_gamma_rvs_corr, var_richness_predicted_gamma_rvs_corr, mean_diversity_observed_gamma_rvs, var_diversity_observed_gamma_rvs, mean_diversity_predicted_gamma_rvs, var_diversity_predicted_gamma_rvs, mean_diversity_predicted_gamma_rvs_corr, var_diversity_predicted_gamma_rvs_corr = diversity_utils.predict_mean_and_var_richness_and_diversity_using_gamma_rv(s_by_s_genera)

        mean_diversity_observed, mean_diversity_predicted, var_diversity_observed, var_diversity_predicted = diversity_utils.predict_mean_and_var_diversity_analytic(s_by_s_genera, s_by_s_genera.shape[0])

        print(var_diversity_predicted, var_diversity_predicted_gamma_rvs)


        color = color_environment_taxon[rank_idx+1]

        continue

        #print(var_richness_predicted_gamma_rvs, variance_prediction_test)

        if numpy.isnan(variance_prediction_test) == True:
            continue

        #print(var_diversity_observed_gamma_rvs, var_diversity_predicted_gamma_rvs_corr)

        ax.scatter(var_richness_predicted_gamma_rvs, variance_prediction_test, color=color, s=40, linewidth=0.8, edgecolors='k', zorder=2)

        #variance_richness_observed_all.append(variance_richness_observed)
        #variance_richness_predicted_all.append(variance_richness_predicted)

        variance_richness_observed_all.append(var_richness_predicted_gamma_rvs)
        variance_richness_predicted_all.append(variance_prediction_test)

        



min_ = min(variance_richness_observed_all + variance_richness_predicted_all)
max_ = max(variance_richness_observed_all + variance_richness_predicted_all)

ax.plot([min_*0.5,max_*1.1],[min_*0.5,max_*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
ax.set_xlim([min_*0.5,max_*1.1])
ax.set_ylim([min_*0.5,max_*1.1])


ax.set_xscale('log', base=10)
ax.set_yscale('log', base=10)
#ax_occupancy_taxon.plot([min_,max_],[min_,max_], lw=4, ls=':',c='k', zorder=2, label='1:1')
#ax.plot([min_*0.5,1.1],[min_*0.5,1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')


ax.set_xlabel("Predicted variance of richness, simulation", fontsize = 9)
ax.set_ylabel("Predicted variance of diversity, analytic", fontsize = 9)

fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sgamma_diversity_mean_taxon.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
