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


rarefied = False



#dbd_utils.make_richness_diversity_prediction_taxon_dict()

#dbd_utils.make_richness_diversity_prediction_phylo_dict()

#dict_taxon = dbd_utils.load_richness_diversity_prediction_taxon_dict()

#print(dict_taxon)

#fig, ax = plt.subplots(figsize=(4,4))

fig = plt.figure(figsize = (11, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

variance_richness_observed_all = []
variance_richness_predicted_all = []

environments_to_keep = diversity_utils.environments_to_keep

environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):    

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)

        samples = sad_annotated_dict['samples']
        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        sys.stderr.write("Getting site-by-species matrix...\n")
        #s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)

        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

        for rank_idx, rank in enumerate([diversity_utils.taxa_ranks[0]]):

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


            rel_s_by_s_genera = s_by_s_genera/numpy.sum(s_by_s_genera, axis=0)
            mean_rel_s_by_s_genera = numpy.mean(rel_s_by_s_genera, axis=1)
            var_rel_s_by_s_genera = numpy.var(rel_s_by_s_genera, axis=1)
            corr_rel_s_by_s_genera = numpy.corrcoef(rel_s_by_s_genera)
            total_reads_genera = numpy.sum(s_by_s_genera, axis=0)

            
            corr_flat = (corr_rel_s_by_s_genera[numpy.triu_indices(corr_rel_s_by_s_genera.shape[0], k = 1)])**2
            corr_flat_null_all = []
            for i in range(100):

                rel_abundances_gamma_genera, read_counts_gamma_genera = simulation_utils.genrate_community_from_mean_and_var(mean_rel_s_by_s_genera, var_rel_s_by_s_genera, total_reads_genera, len(total_reads_genera), corr_matrix=corr_rel_s_by_s_genera)

                rel_read_counts_gamma_genera = read_counts_gamma_genera/numpy.sum(read_counts_gamma_genera, axis=0)
                corr_rel_read_counts_gamma_genera = numpy.corrcoef(rel_read_counts_gamma_genera)

                corr_flat_null_all.append((corr_rel_read_counts_gamma_genera[numpy.triu_indices(corr_rel_read_counts_gamma_genera.shape[0], k = 1)])**2)

            
            mean_corr_flat_null = numpy.mean(corr_flat_null_all, axis=0)
            
            idx_to_keep = (mean_corr_flat_null>0) & (corr_flat>0)&(~numpy.isnan(mean_corr_flat_null))&(~numpy.isnan(corr_flat))
 
            corr_flat = corr_flat[idx_to_keep]
            mean_corr_flat_null = mean_corr_flat_null[idx_to_keep]

            print(corr_flat)
            print(mean_corr_flat_null)

            ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

            print(min(corr_flat), min(mean_corr_flat_null))

            ax.scatter(corr_flat, mean_corr_flat_null, alpha=0.1, s=5)
            ax.set_xscale('log', base=10)
            ax.set_yscale('log', base=10)
        
            ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)
            ax.set_xlabel("Observed rho^2", fontsize = 10)
            ax.set_ylabel("Predicted rho^2", fontsize = 10)


            min_ = min(corr_flat.tolist() + mean_corr_flat_null.tolist())
            max_ = max(corr_flat.tolist() + mean_corr_flat_null.tolist())

            ax.plot([min_*0.5,max_*1.1],[min_*0.5,max_*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
            ax.set_xlim([min_*0.5,max_*1.1])
            ax.set_ylim([min_*0.5,max_*1.1])




fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%stest_rho_prediction.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()

