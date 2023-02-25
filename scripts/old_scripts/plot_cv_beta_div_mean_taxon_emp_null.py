import os
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
import plot_utils


rarefied = False
iter = 1000


environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label

cv_beta_div_mean_dict_path = config.data_directory + "cv_beta_div_mean_dict_taxon.pickle"

color_all = numpy.asarray([plot_utils.rgb_blue_taxon(r) for r in idx_taxa_ranks])
environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]

def make_cv_beta_div_mean_dict():

    cv_beta_div_mean_dict = {}

    for environment in environments_to_keep:

        cv_beta_div_mean_dict[environment] = {}

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
        rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

        mean_rel_abundance = numpy.mean(rel_s_by_s, axis=1)
        std_dev_rel_abundance = numpy.std(rel_s_by_s, axis=1)
        beta_div_mean_rel_abundance = mean_rel_abundance/(std_dev_rel_abundance**2)

        # dictionary to map taxon names to higher order rank
        sys.stderr.write("Running taxonomic coarse-graining.....\n")
        for rank_idx, rank in enumerate(taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)
            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))

            genus_to_taxa_dict = {}
            cv_beta_div_mean_g_taxa_all = []
            g_taxa_idx_all = []
            for genus_idx, genus in enumerate(all_genera):
                genus_to_taxa_dict[genus] = []
                for t in taxa:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])

                if len(g_taxa_idx) < 5:
                    continue

                beta_div_mean_g_taxa = beta_div_mean_rel_abundance[g_taxa_idx]
                cv_beta_div_mean_g_taxa = numpy.std(beta_div_mean_g_taxa)/numpy.mean(beta_div_mean_g_taxa)
                cv_beta_div_mean_g_taxa_all.append(cv_beta_div_mean_g_taxa)
                g_taxa_idx_all.append(g_taxa_idx)

            # calculate mean of CV
            mean_cv_beta_div_mean_g_taxa_all = numpy.mean(cv_beta_div_mean_g_taxa_all)
            # copy
            beta_div_mean_rel_abundance_copy = numpy.copy(beta_div_mean_rel_abundance)
            cv_beta_div_mean_g_taxa_null_all = []
            for i in range(iter):

                numpy.random.shuffle(beta_div_mean_rel_abundance_copy)

                cv_beta_div_mean_g_taxa_null_i = numpy.mean([numpy.std(beta_div_mean_rel_abundance_copy[g_taxa_idx])/numpy.mean(beta_div_mean_rel_abundance_copy[g_taxa_idx]) for g_taxa_idx in g_taxa_idx_all])
                cv_beta_div_mean_g_taxa_null_all.append(cv_beta_div_mean_g_taxa_null_i)

            cv_beta_div_mean_g_taxa_null_all = numpy.asarray(cv_beta_div_mean_g_taxa_null_all)
            cv_beta_div_mean_g_taxa_null_all = numpy.sort(cv_beta_div_mean_g_taxa_null_all)

            lower_ci = cv_beta_div_mean_g_taxa_null_all[int(iter*0.025)]
            upper_ci = cv_beta_div_mean_g_taxa_null_all[int(iter*0.975)]

            cv_beta_div_mean_dict[environment][rank] = {}
            cv_beta_div_mean_dict[environment][rank]['mean_cv'] = mean_cv_beta_div_mean_g_taxa_all
            cv_beta_div_mean_dict[environment][rank]['lower_ci_cv'] = lower_ci
            cv_beta_div_mean_dict[environment][rank]['upper_ci_cv'] = upper_ci



    sys.stderr.write("Saving slope dictionary...\n")
    with open(cv_beta_div_mean_dict_path, 'wb') as handle:
        pickle.dump(cv_beta_div_mean_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




def load_cv_beta_div_mean_dict():

    with open(cv_beta_div_mean_dict_path, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_



#make_cv_beta_div_mean_dict()

dict_ = load_cv_beta_div_mean_dict()



# make the plot
fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)



for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        mean_cv = numpy.asarray([dict_[environment][r]['mean_cv'] for r in taxa_ranks])
        lower_ci_cv = numpy.asarray([dict_[environment][r]['lower_ci_cv'] for r in taxa_ranks])
        upper_ci_cv= numpy.asarray([dict_[environment][r]['upper_ci_cv'] for r in taxa_ranks])


        ax.plot(idx_taxa_ranks, mean_cv, ls='-', lw=1.5, c='k',  zorder=2)
        ax.scatter(idx_taxa_ranks, mean_cv, c=color_all, s=60, label='Observed', linewidth=0.8, edgecolors='k', zorder=3)
        ax.fill_between(idx_taxa_ranks, lower_ci_cv, upper_ci_cv, color='darkgrey', alpha=0.7, label='95% Prediction null interval', zorder=1)

        ax.set_title(diversity_utils.format_environment_label(environment), fontsize=11)
        ax.set_xticks(idx_taxa_ranks)
        ax.set_xticklabels(taxa_ranks_label, fontsize=8)

        ax.axhline(y=1, color='k', linestyle=':', lw = 2, zorder=1, label=r'$\mathrm{CV} = 1$')

        ax.set_ylim([0.6,2])

        if environment_idx == 0:
            ax.set_ylabel("Mean CV of gamma rate parameters\nin a coarse-grained group, " + r'$\mathrm{CV}_{\beta/\bar{x}}$', fontsize = 10)


        if (environment_chunk_idx ==0) and (environment_idx == 0):
            ax.legend(loc="upper left", fontsize=7)



fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%scv_beta_div_mean_taxon_emp_null%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
