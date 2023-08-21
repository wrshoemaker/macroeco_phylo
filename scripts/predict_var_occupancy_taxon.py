from __future__ import division
import config
import biom
import numpy
import os
os.environ["SYMPY_USE_CACHE"]="no"
import sympy
#import math
import simulation_utils
import scipy.integrate as integrate
import scipy.special as special
#import dbd_utils


import numpy
import diversity_utils
import simulation_utils

import scipy.stats as stats
import scipy.integrate as integrate

#import ete3
import tree_utils
import random
import dbd_utils
import sys

import matplotlib.pyplot as plt
import plot_utils



def predict_var_occupancy(s_by_s, species, totreads=numpy.asarray([])):

    # get squared inverse cv
    # assume that entries are read counts.
    rel_s_by_s_np = (s_by_s/s_by_s.sum(axis=0))

    beta_all = []
    mean_all = []

    for s in rel_s_by_s_np:

        var = numpy.var(s)
        mean = numpy.mean(s)

        beta = (mean**2)/var

        mean_all.append(mean)
        beta_all.append(beta)

    beta_all = numpy.asarray(beta_all)
    mean_all = numpy.asarray(mean_all)


    s_by_s_presence_absence = numpy.where(s_by_s > 0, 1, 0)

    occupancies = s_by_s_presence_absence.sum(axis=1) / s_by_s_presence_absence.shape[1]

    # calcualte total reads if no argument is passed
    # sloppy quick fix
    if len(totreads) == 0:
        totreads = s_by_s.sum(axis=0)

    # calculate mean and variance excluding zeros
    # tf = mean relative abundances
    tf = []
    for afd in s_by_s:
        afd_no_zeros = afd[afd>0]
        tf.append(numpy.mean(afd_no_zeros/ totreads[afd>0]))

    tf = numpy.asarray(tf)
    # go through and calculate the variance for each species

    tvpf_list = []
    for afd in s_by_s:
        afd_no_zeros = afd[afd>0]

        N_reads = s_by_s.sum(axis=0)[numpy.nonzero(afd)[0]]
        tvpf_list.append(numpy.mean(  (afd_no_zeros**2 - afd_no_zeros) / (totreads[afd>0]**2) ))

    tvpf = numpy.asarray(tvpf_list)

    f = occupancies*tf
    vf= occupancies*tvpf

    # there's this command in Jacopo's code %>% mutate(vf = vf - f^2 )%>%
    # It's applied after f and vf are calculated, so I think I can use it
    # This should be equivalent to the mean and variance including zero
    vf = vf - (f**2)

    beta = (f**2)/vf
    theta = f/beta

    predicted_occupancies = []
    predicted_second_moment_occupancies = []
    # each species has it's own beta and theta, which is used to calculate predicted occupancy
    for beta_i, theta_i in zip(beta,theta):
        predicted_occupancies.append(1 - numpy.mean( ((1+theta_i*totreads)**(-1*beta_i ))   ))
        predicted_second_moment_occupancies.append(numpy.mean((1 -  ((1+theta_i*totreads)**(-1*beta_i ))) ** 2) )


    predicted_occupancies = numpy.asarray(predicted_occupancies)
    predicted_second_moment_occupancies = numpy.asarray(predicted_second_moment_occupancies)

    species = numpy.asarray(species)
    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
    mad = numpy.mean(rel_s_by_s, axis=1)

    predicted_variance_occupancy = predicted_second_moment_occupancies - predicted_occupancies**2

    # calculate observed variance of occupancy
    # observed second moment occupancy
    # <o^2> = 1/M \sum( \delta^2_{i} ) = <o>?
    # delta is kronecker delta
    #  

    observed_variance_occupancy = occupancies - (occupancies**2)

    # make sure both values greater than zero
    idx_to_keep = (predicted_variance_occupancy > 0) & (observed_variance_occupancy>0)

    observed_variance_occupancy = observed_variance_occupancy[idx_to_keep]
    predicted_variance_occupancy = predicted_variance_occupancy[idx_to_keep]


    return observed_variance_occupancy, predicted_variance_occupancy




rarefied = False

environments_to_keep = diversity_utils.environments_to_keep
#environments_to_keep = [environments_to_keep[0]]

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label
taxa_ranks.insert(0, 'OTU')



fig = plt.figure(figsize = (12, 20)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


for environment_idx, environment in enumerate(environments_to_keep):

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)
    samples = sad_annotated_dict['samples']

    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    for rank_idx, rank in enumerate(taxa_ranks):

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

        
        observed, predicted = predict_var_occupancy(s_by_s_genera, range(s_by_s_genera.shape[0]))


        idx_to_keep = (observed>0) & (predicted > 0)
        observed = observed[idx_to_keep]
        predicted = predicted[idx_to_keep]

        predicted_and_observed_occupancies = numpy.concatenate((observed,predicted), axis=0)
        min_ = min(predicted_and_observed_occupancies)
        max_ = max(predicted_and_observed_occupancies)


        ax = plt.subplot2grid((9, 6), (environment_idx, rank_idx))



        x, y, z = plot_utils.get_scatter_density_arrays_for_loglog(observed, predicted)
        ax.scatter(x, y, c=numpy.sqrt(z), cmap='Blues', s=70, alpha=0.9, edgecolors='none', zorder=1)

        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)


        ax.plot([min_,max_],[min_,max_], lw=3,ls=':',c='k',zorder=1)

        if environment_idx == 0:
            ax.set_title(diversity_utils.taxa_ranks_label_with_asv[rank_idx], fontsize=12)

        if rank_idx == 0:
            ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize = 12)


        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)


        #ax.set_xlabel("Sum of ASV variances", fontsize = 9)
        #ax.set_ylabel("Coarse-grained variance", fontsize = 9)


fig.text(0.3, 0.06, "Observed variance of occupancy", va='center', fontsize=28)
fig.text(0.02, 0.5, "Predicted variance of occupancy", va='center', rotation='vertical', fontsize=28)


fig.subplots_adjust(hspace=0.37, wspace=0.37)
fig_name = "%sgamma_occupancy_var_all_taxon%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()








