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


coarse_grained_tree_dict_all = {}
distances_all = []
for environment in environments_to_keep:
    sys.stderr.write("Loading tree dict for %s...\n" % environment)
    coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)
    coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict
    distances = list(coarse_grained_tree_dict.keys())
    distances_all.extend(distances)

distances_all = list(set(distances_all))
distances_all.sort()
color_all = plot_utils.make_blue_cmap(len(distances_all))



distances = numpy.asarray([0.01663614249384222, 0.02767612370754228, 0.05938601867590266, 0.16435747993726982, 0.3526699214174661])
#distances = numpy.asarray([0.3526699214174661])

#color_all = plot_utils.make_blue_cmap(len(distances_all))



fig = plt.figure(figsize = (12, 20)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


for environment_idx, environment in enumerate(environments_to_keep):

    coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    #samples = diversity_utils.subset_observations(environment=environment)
    pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = rarefied)

    samples = pres_abs_dict['samples']
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
    sys.stderr.write("Getting site-by-species matrix...\n")
    #s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)

    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

    for distance_idx, distance in enumerate(distances):

        sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
        coarse_grained_list = coarse_grained_tree_dict[distance]
        coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

        # get indexes for each clade
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)


        occupancies, predicted_occupancies = predict_var_occupancy(s_by_s_all_clades, s_by_s_all_clades.shape[0])

        idx_to_keep = (occupancies>0) & (predicted_occupancies > 0)
        occupancies = occupancies[idx_to_keep]
        predicted_occupancies = predicted_occupancies[idx_to_keep]

        predicted_and_observed_occupancies = numpy.concatenate((occupancies, predicted_occupancies),axis=0)
        min_ = min(predicted_and_observed_occupancies)
        max_ = max(predicted_and_observed_occupancies)


        ax = plt.subplot2grid((9, 6), (environment_idx, distance_idx))

        x, y, z = plot_utils.get_scatter_density_arrays_for_loglog(occupancies, predicted_occupancies)
        ax.scatter(x, y, c=numpy.sqrt(z), cmap='Blues', s=70, alpha=0.9, edgecolors='none', zorder=1)

        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)

        ax.plot([min_,max_],[min_,max_], lw=3,ls=':',c='k',zorder=1)

        if environment_idx == 0:
            ax.set_title('Phylo. dist. = %s' % str(round(distance, 3)), fontsize=11)

        if distance_idx == 0:
            ax.set_ylabel(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize = 12)


        ax.tick_params(axis='x', labelsize=8)
        ax.tick_params(axis='y', labelsize=8)


        #ax.set_xlabel("Sum of ASV variances", fontsize = 9)
        #ax.set_ylabel("Coarse-grained variance", fontsize = 9)


fig.text(0.3, 0.06, "Observed variance of occupancy", va='center', fontsize=28)
fig.text(0.02, 0.5, "Predicted variance of occupancy", va='center', rotation='vertical', fontsize=28)


fig.subplots_adjust(hspace=0.37, wspace=0.37)
fig_name = "%sgamma_occupancy_var_all_phylo%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()


