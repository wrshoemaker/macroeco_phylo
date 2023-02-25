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
import matplotlib as mpl
import dbd_utils

import diversity_utils
import config
import plot_utils



rarefied = False
iter = 1000

rank_all_idx = range(len(diversity_utils.taxa_ranks))
color_all = [plot_utils.rgb_blue_taxon(rank_idx) for rank_idx in rank_all_idx]
taxa_ranks = diversity_utils.taxa_ranks[::-1]
ranks_all_format = [x.capitalize() for x in taxa_ranks]

fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

environments_to_keep = diversity_utils.environments_to_keep

environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        print(environment)

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        #samples = diversity_utils.subset_observations(environment=environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_axon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        sys.stderr.write("Running taxonomic coarse-graining.....\n")
        slope_all = []
        lower_ci_all = []
        upper_ci_all = []
        for rank_idx, rank in enumerate(taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)

            taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))

            s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
            rel_s_by = (s_by_s/s_by_s.sum(axis=0))

            genus_to_taxa_dict = {}
            sad_genera_all = []

            beta_div_mean_ratio_all = []
            for genus in all_genera:
                genus_to_taxa_dict[genus] = []
                for t in taxa:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
                sad_genera_all.append(s_by_s[g_taxa_idx,:].sum(axis=0))


            s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
            # remove sites where there are no observations
            s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]
            rel_s_by_s_genera = (s_by_s_genera/s_by_s_genera.sum(axis=0))

            mean_rel_s_by_s_genera = numpy.mean(rel_s_by_s_genera, axis=1)
            var_rel_s_by_s_genera = numpy.var(rel_s_by_s_genera, axis=1)


            slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(mean_rel_s_by_s_genera), numpy.log10(var_rel_s_by_s_genera))
            slope_all.append(slope)


            # generate null CIs
            # get indices for null genera coarse graining
            unique_genera, counts_genera = numpy.unique([sad_annotated_dict['taxa'][t][rank] for t in taxa], return_counts=True)
            coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_genera))[:-1]

            s_by_s_copy = numpy.copy(s_by_s)
            slope_null_all = []
            for i in range(iter):

                numpy.random.shuffle(s_by_s_copy)
                s_by_s_coarse_null = numpy.add.reduceat(s_by_s_copy, coarse_grain_by_genus_idx, axis=0)
                s_by_s_coarse_null = s_by_s_coarse_null[:,~(numpy.all(s_by_s_coarse_null == 0, axis=0))]
                rel_s_by_s_coarse_null = (s_by_s_coarse_null/s_by_s_coarse_null.sum(axis=0))

                mean_rel_s_by_s_genera_null = numpy.mean(rel_s_by_s_coarse_null, axis=1)
                var_rel_s_by_s_genera_null = numpy.var(rel_s_by_s_coarse_null, axis=1)

                #diversity_genera_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_coarse_null)
                slope_null, intercept_null, r_value_null, p_value_null, std_err_null = stats.linregress(numpy.log10(mean_rel_s_by_s_genera_null), numpy.log10(var_rel_s_by_s_genera_null))
                slope_null_all.append(slope_null)

            slope_null_all = numpy.asarray(slope_null_all)
            slope_null_all = numpy.sort(slope_null_all)

            percentile = sum(slope_null_all < slope)/iter

            lower_ci = slope_null_all[int(iter*0.025)]
            upper_ci = slope_null_all[int(iter*0.975)]

            lower_ci_all.append(lower_ci)
            upper_ci_all.append(upper_ci)

            print(slope, lower_ci, upper_ci)



        ax.fill_between(rank_all_idx, lower_ci_all, upper_ci_all, color='darkgrey', alpha=0.7, label='95% CI', zorder=1)


        #ax.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=diversity_utils.rgb_blue_taxon(rank_idx), alpha=0.9, lw=2, label=rank)
        ax.plot(rank_all_idx, slope_all, c='k', lw=2, ls='-', alpha=1, zorder=2)  # , c='#87CEEB')
        ax.scatter(rank_all_idx, slope_all, c=color_all, edgecolors='k', alpha=1, zorder=3)  # , c='#87CEEB')



        ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)
        #ax.set_yscale('log', base=10)

        if environment_idx == 0:
            ax.set_ylabel("Slope of Taylor's Law", fontsize = 10)


        #if environment_chunk_idx == len(environment_chunk_all)-1:
        ax.set_xticks(rank_all_idx)
        ax.set_xticklabels(ranks_all_format, fontsize=9)
        ax.set_xlim([-0.3,4.3])



        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="upper left", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%scoarse_taylor_taxon%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
