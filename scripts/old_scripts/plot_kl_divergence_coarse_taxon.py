
import os
#from Bio import Phylo
import random
import copy
import sys
import numpy
import random
import pickle
import scipy.stats as stats
import scipy.special as special
import matplotlib.pyplot as plt
import functools
import operator

import dbd_utils

import diversity_utils
import config

from scipy.special import rel_entr
import simulation_utils


if len(sys.argv) > 1:
    x_axis_n = True
else:
    x_axis_n = False

rarefied = False


kl_divergence_taxon_path = "%skl_divergence_taxon.pickle" % config.data_directory


environments_to_keep = diversity_utils.environments_to_keep





def gamma_kl_integral(shape_1, rate_1, shape_2, rate_2):

    scale_1 = 1/rate_1
    scale_2 = 1/rate_2

    return (-1*shape_2*scale_2/shape_1) - numpy.log((shape_1**(scale_1))  * special.gamma(scale_1)) + (scale_1-1)*(polygamma(0,rate_2) +  numpy.log(shape_2))



def calculate_kl_divergence(array_1, array_2):

    idx_to_keep = (array_1>0) & (array_2>0)

    array_1 = array_1[idx_to_keep]
    array_2 = array_2[idx_to_keep]

    array_1 = array_1/sum(array_1)
    array_2 = array_2/sum(array_2)

    rel_entr_sum = sum(rel_entr(array_1, array_2))

    return rel_entr_sum





def make_kl_prediction_dict():

    kl_prediction_dict = {}

    for environment_idx, environment in enumerate(environments_to_keep):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #samples = diversity_utils.subset_observations(environment=environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

        mean = numpy.mean(rel_s_by_s, axis=1)
        var = numpy.var(rel_s_by_s, axis=1)

        N = s_by_s.sum(axis=0)

        kl_prediction_dict[environment] = {}
        for taxa_i in taxa:
            kl_prediction_dict[environment][taxa_i] = {}
            for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):
                kl_prediction_dict[environment][taxa_i][rank] = {}
                kl_prediction_dict[environment][taxa_i][rank]['kl_divergence'] = []
                kl_prediction_dict[environment][taxa_i][rank]['kl_divergence_coarse'] = []


        rank_parameter_dict = {}
        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)
            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))

            genus_to_taxa_dict = {}
            sad_genera_all = []
            g_taxa_idx_all = []
            for genus_idx, genus in enumerate(all_genera):
                genus_to_taxa_dict[genus] = []
                for t in taxa:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])

                s_by_s_genus = s_by_s[g_taxa_idx,:]
                sad_genera_all.append(s_by_s_genus.sum(axis=0))
                g_taxa_idx_all.append(g_taxa_idx)


            s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
            rel_s_by_s_genera = (s_by_s_genera/s_by_s_genera.sum(axis=0))
            mean_genera = numpy.mean(rel_s_by_s_genera, axis=1)
            var_genera = numpy.var(rel_s_by_s_genera, axis=1)

            rank_parameter_dict[rank] = {}
            rank_parameter_dict[rank]['mean'] = mean_genera
            rank_parameter_dict[rank]['var'] = var_genera
            rank_parameter_dict[rank]['g_taxa_idx_all'] = g_taxa_idx_all


            for g_taxa_idx_idx, g_taxa_idx_i in enumerate(g_taxa_idx_all):

                for g_taxa_idx_i_j in g_taxa_idx_i:

                    kl_prediction_dict[environment][taxa[g_taxa_idx_i_j]][rank]['n'] = len(g_taxa_idx_i)

        sys.stderr.write("Starting simulation...\n")
        for i in range(100):

            s_by_s_sim, s_by_s_sim_counts = simulation_utils.genrate_community_from_mean_and_var(mean, var, N, len(N))
            rel_s_by_s_sim_counts = (s_by_s_sim_counts/s_by_s_sim_counts.sum(axis=0))

            for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

                s_by_s_genera_sim, s_by_s_genera_sim_counts = simulation_utils.genrate_community_from_mean_and_var(rank_parameter_dict[rank]['mean'], rank_parameter_dict[rank]['var'], N, len(N))

                rel_s_by_s_genera_sim_counts = (s_by_s_genera_sim_counts/s_by_s_genera_sim_counts.sum(axis=0))

                g_taxa_idx_all = rank_parameter_dict[rank]['g_taxa_idx_all']

                # here we are getting a prediction by coarse-graining the simulated ASVs
                s_by_s_sim_counts_coarse = [s_by_s_sim[g_taxa_idx_i,:].sum(axis=0) for g_taxa_idx_i in g_taxa_idx_all]
                s_by_s_sim_counts_coarse = numpy.asarray(s_by_s_sim_counts_coarse)
                rel_s_by_s_sim_counts_coarse = (s_by_s_sim_counts_coarse/s_by_s_sim_counts_coarse.sum(axis=0))


                for g_taxa_idx_idx, g_taxa_idx_i in enumerate(g_taxa_idx_all):

                    sad_genus_i = rel_s_by_s_genera_sim_counts[g_taxa_idx_idx,:]
                    sad_genus_coarse_i = rel_s_by_s_sim_counts_coarse[g_taxa_idx_idx,:]

                    for g_taxa_idx_i_j in g_taxa_idx_i:
                        sad_asv_i_j = rel_s_by_s_sim_counts[g_taxa_idx_i_j,:]

                        rel_entr_sum = calculate_kl_divergence(sad_asv_i_j, sad_genus_i)
                        rel_entr_sum_coarse = calculate_kl_divergence(sad_asv_i_j, sad_genus_coarse_i)

                        #if kl_prediction_dict[environment][taxa[g_taxa_idx_i_j]][rank]['n']>100:

                        #    print(sad_asv_i_j)
                        #    print(sad_genus_i)

                        #    print(kl_prediction_dict[environment][taxa[g_taxa_idx_i_j]][rank]['n'], rel_entr_sum)
                        taxa_i_j = taxa[g_taxa_idx_i_j]
                        kl_prediction_dict[environment][taxa_i_j][rank]['kl_divergence'].append(rel_entr_sum)
                        kl_prediction_dict[environment][taxa_i_j][rank]['kl_divergence_coarse'].append(rel_entr_sum_coarse)


        # get the mean
        for taxa_i in taxa:
            for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):
                kl_divergence = numpy.asarray(kl_prediction_dict[environment][taxa_i][rank]['kl_divergence'])
                kl_divergence = kl_divergence[kl_divergence>0]

                if len(kl_divergence) == 0:
                    mean_kl_divergence = 0
                else:
                    mean_kl_divergence = numpy.mean(kl_divergence)


                kl_divergence_coarse = numpy.asarray(kl_prediction_dict[environment][taxa_i][rank]['kl_divergence_coarse'])
                kl_divergence_coarse = kl_divergence_coarse[kl_divergence_coarse>0]

                if len(kl_divergence_coarse) == 0:
                    mean_kl_divergence_coarse = 0
                else:
                    mean_kl_divergence_coarse = numpy.mean(kl_divergence_coarse)

                kl_prediction_dict[environment][taxa_i][rank]['mean_kl_divergence'] = mean_kl_divergence
                del kl_prediction_dict[environment][taxa_i][rank]['kl_divergence']

                kl_prediction_dict[environment][taxa_i][rank]['mean_kl_divergence_coarse'] = mean_kl_divergence_coarse
                del kl_prediction_dict[environment][taxa_i][rank]['kl_divergence_coarse']



    sys.stderr.write("Saving diversity vs. diversity dictionary...\n")
    with open(kl_divergence_taxon_path, 'wb') as handle:
        pickle.dump(kl_prediction_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_kl_prediction_dict():

    with open(kl_divergence_taxon_path, 'rb') as handle:
        kl_prediction_dict = pickle.load(handle)
    return kl_prediction_dict




make_kl_prediction_dict()



kl_prediction_dict = load_kl_prediction_dict()


fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)



#environments_to_keep = [environments_to_keep[0]]

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

        sys.stderr.write("Running taxonomic coarse-graining.....\n")
        afd_dict = {}
        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

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

                        if t not in afd_dict:
                            afd_dict[t] = {}

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])

                s_by_s_genus = s_by_s[g_taxa_idx,:]
                sad_genera_all.append(s_by_s_genus.sum(axis=0))
                rel_s_by_s_genus = rel_s_by_s[g_taxa_idx,:]

                for t_idx, t in enumerate(genus_to_taxa_dict[genus]):
                    if 'otu' not in afd_dict[t]:
                        afd_dict[t]['otu'] = rel_s_by_s_genus[t_idx,:]


            s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
            # remove sites where there are no observations
            s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]
            rel_s_by_s_genera = (s_by_s_genera/s_by_s_genera.sum(axis=0))

            for genus_idx, genus in enumerate(all_genera):
                for t_idx, t in enumerate(genus_to_taxa_dict[genus]):
                    afd_dict[t][rank] = {}
                    afd_dict[t][rank]['n'] = len(genus_to_taxa_dict[genus])
                    afd_dict[t][rank]['name'] = genus
                    afd_dict[t][rank]['afd'] = rel_s_by_s_genera[genus_idx,:]

        n_all = []
        mean_kl_divergence_all = []
        for t in taxa:
            for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

                n_all.append(kl_prediction_dict[environment][t][rank]['n'])
                mean_kl_divergence_all.append(kl_prediction_dict[environment][t][rank]['mean_kl_divergence'])

        n_all, mean_kl_divergence_all = zip(*sorted(zip(n_all, mean_kl_divergence_all)))

        n_all = numpy.asarray(n_all)
        mean_kl_divergence_all = numpy.asarray(mean_kl_divergence_all)

        n_all_to_plot = n_all[mean_kl_divergence_all>0]
        mean_kl_divergence_all_to_plot = mean_kl_divergence_all[mean_kl_divergence_all>0]

        n_all_to_plot_log10 = numpy.log10(n_all_to_plot)
        mean_kl_divergence_all_to_plot_log10 = numpy.log10(mean_kl_divergence_all_to_plot)

        hist, bin_edges = numpy.histogram(n_all_to_plot_log10, density=True, bins=10)
        bins_mean = [0.5 * (bin_edges[i] + bin_edges[i+1]) for i in range(0, len(bin_edges)-1 )]

        bins_mean_to_plot = []
        mean_kl_divergence_bin_to_plot = []
        for i in range(0, len(bin_edges)-1):
            idx_i = (n_all_to_plot_log10 > bin_edges[i]) & (n_all_to_plot_log10 <= bin_edges[i+1])

            if len(idx_i) == 0:
                continue

            bins_mean_to_plot.append(bins_mean[i])

            mean_kl_divergence_bin_to_plot.append(numpy.mean(mean_kl_divergence_all_to_plot_log10[idx_i]))

        bins_mean_to_plot = numpy.asarray(bins_mean_to_plot)
        mean_kl_divergence_bin_to_plot = numpy.asarray(mean_kl_divergence_bin_to_plot)

        bins_mean_to_plot = 10**bins_mean_to_plot
        mean_kl_divergence_bin_to_plot = 10**mean_kl_divergence_bin_to_plot

        idx_to_keep = ((~numpy.isnan(bins_mean_to_plot)) & (~numpy.isnan(mean_kl_divergence_bin_to_plot)))

        bins_mean_to_plot = bins_mean_to_plot[idx_to_keep]
        mean_kl_divergence_bin_to_plot = mean_kl_divergence_bin_to_plot[idx_to_keep]


        # plot
        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        asv_all = list(afd_dict.keys())
        finished_rank_dict = {}
        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):
            finished_rank_dict[rank] = []

        for asv in asv_all:

            afd_asv = afd_dict[asv]['otu']

            if sum(afd_asv>0) < 20:
                continue

            # run ranks
            rank_idx_to_plot = []
            n_to_plot = []
            divergence_to_plot = []
            for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

                name_rank = afd_dict[asv][rank]['name']
                afd_rank = afd_dict[asv][rank]['afd']
                n = afd_dict[asv][rank]['n']

                rel_entr_sum = calculate_kl_divergence(afd_asv, afd_rank)

                if rel_entr_sum > 0:

                    rank_idx_to_plot.append(rank_idx)
                    divergence_to_plot.append(rel_entr_sum)
                    n_to_plot.append(n)


            if x_axis_n == True:
                ax.plot(n_to_plot, divergence_to_plot, color='k', alpha=0.05, lw=1)

            else:
                ax.plot(rank_idx_to_plot, divergence_to_plot, color='k', alpha=0.05, lw=1)


        ax.plot(bins_mean_to_plot, mean_kl_divergence_bin_to_plot)


        ax.set_ylim([10**-3,10**1])
        ax.set_yscale('log', base=10)
        ax.set_ylabel("KL divergence b/w ASV and coarse-grained AFD", fontsize=7.5)
        ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)

        if x_axis_n == True:
            ax.set_xscale('log', base=10)
            ax.set_xlabel("Number coarse-grained ASVs", fontsize=9)


        else:
            taxa_ranks = diversity_utils.taxa_ranks
            idx_taxa_ranks = list(range(len(taxa_ranks)))
            taxa_ranks_label = diversity_utils.taxa_ranks_label

            ax.set_xticks(idx_taxa_ranks)
            ax.set_xticklabels(taxa_ranks_label, fontsize=8)





if x_axis_n == True:
    x_axis_n_label = '_n'
else:
    x_axis_n_label = ''



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%skl_divergence_taxon%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied) + x_axis_n_label)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
