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
from collections import Counter

from scipy.special import rel_entr

import diversity_utils

from scipy.stats import hypergeom


#n_subsampled = 10
#N_subsampled = 100
#n =100
#N = 10000

#print(hypergeom.pmf(n_subsampled, N, n, N_subsampled))


sample_to_environment_dict = diversity_utils.get_sample_to_environment_dict()

n_sites = 2000
environment = 'human gut metagenome'

file_path = '%semp/otu_distributions/otu_summary.emp_deblur_90bp.subset_2k.rare_5000.tsv' % config.data_directory
file = open(file_path, 'r')
header = file.readline()
header_split = header.strip().split('\t')


environments = []
sites_all = []

environment_dict = {}

count_dict = {}
for line in file:

    line_split = line.strip().split('\t')

    sites = line_split[-1].split(',')

    for s in sites:
        environment_s = sample_to_environment_dict[s]

        if environment_s not in environment_dict:
            environment_dict[environment_s] = []

        environment_dict[environment_s].append(s)

        if environment == environment_s:

            sites_all.append(s)


file.close()


sites_all = numpy.asarray(sites_all)
numpy.random.shuffle(sites_all)

#for environment in diversity_utils.environments_to_keep:

#    sites_to_keep = list(set(environment_dict[environment]))
#    s_by_s, taxa, samples = diversity_utils.get_s_by_s_deblur(sites_to_keep)

#    print(environment, s_by_s.shape[1])



sites_to_keep = sites_all[:n_sites]

s_by_s_subsample, taxa_subsample, sites_subsample = diversity_utils.get_s_by_s_deblur_subsample(sites_to_keep)
s_by_s_sample, taxa_sample, sites_sample  = diversity_utils.get_s_by_s_deblur(sites_subsample)


n_reads_sample = s_by_s_sample.sum(axis=0)
n_reads_subsample = s_by_s_subsample.sum(axis=0)

rel_s_by_s = s_by_s_sample/s_by_s_sample.sum(axis=0)
rel_s_by_s_subsample = s_by_s_subsample/s_by_s_subsample.sum(axis=0)

mean_sample = numpy.mean(rel_s_by_s, axis=1)
var_sample = numpy.var(rel_s_by_s, axis=1)
beta_sample = (mean_sample**2)/var_sample


mean_subsample = numpy.mean(rel_s_by_s_subsample, axis=1)
var_subsample = numpy.var(rel_s_by_s_subsample, axis=1)
beta_subsample = (mean_subsample**2)/var_subsample



fig = plt.figure(figsize = (8, 4)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


ax_gamma = plt.subplot2grid((1, 2), (0, 0))
ax_gamma_hyper = plt.subplot2grid((1, 2), (0, 1))


observed_occupany_all = []
observed_occupany_subsample_all = []
predicted_occupancy_subsample_all = []
predicted_occupancy_hyper_all = []
for afd_subsample_idx, afd_subsample in enumerate(s_by_s_subsample):

    afd_subsample_taxon = taxa_subsample[afd_subsample_idx]

    afd_idx = numpy.where(taxa_sample==afd_subsample_taxon)[0][0]
    afd_sample = s_by_s_sample[afd_idx,:]

    mean_sample_afd = mean_sample[afd_idx]
    var_sample_afd = var_sample[afd_idx]
    beta_sample_afd = beta_sample[afd_idx]

    mean_subsample_afd = mean_subsample[afd_subsample_idx]
    var_subsample_afd = var_subsample[afd_subsample_idx]
    beta_subsample_afd = beta_subsample[afd_subsample_idx]

    prob_absence_subsample_hyper = hypergeom.pmf(numpy.asarray([0]*len(n_reads_sample)), n_reads_sample, afd_sample, n_reads_subsample)
    prob_absence_subsample_gamma = (beta_subsample_afd/(beta_subsample_afd + n_reads_subsample*mean_subsample_afd))**beta_subsample_afd
    prob_absence_sample_gamma = (beta_sample_afd/(beta_sample_afd + n_reads_sample*mean_sample_afd))**beta_sample_afd

    predicted_occupancy_hyper = 1 - numpy.mean(prob_absence_sample_gamma*prob_absence_subsample_hyper)
    predicted_occupancy = 1 - numpy.mean(prob_absence_sample_gamma)
    predicted_occupancy_subsample = 1 - numpy.mean(prob_absence_subsample_gamma)

    observed_occupany_subsample = sum(afd_subsample>0)/len(afd_subsample)
    observed_occupany = sum(afd_sample>0)/len(afd_sample)

    #error = numpy.absolute(predicted_occupancy-observed_occupany_subsampling)/observed_occupany_subsampling
    #error_subsample = numpy.absolute(predicted_occupancy_subsample-observed_occupany_subsampling)/observed_occupany_subsampling

    #print(observed_occupany, predicted_occupancy, predicted_occupancy_subsample, error, error_subsample)

    observed_occupany_all.append(observed_occupany)
    observed_occupany_subsample_all.append(observed_occupany_subsample)

    predicted_occupancy_subsample_all.append(predicted_occupancy_subsample)
    predicted_occupancy_hyper_all.append(predicted_occupancy_hyper)





observed_occupany_all = numpy.asarray(observed_occupany_all)
observed_occupany_subsample_all = numpy.asarray(observed_occupany_subsample_all)
predicted_occupancy_subsample_all = numpy.asarray(predicted_occupancy_subsample_all)
predicted_occupancy_hyper_all = numpy.asarray(predicted_occupancy_hyper_all)


idx_to_keep = (observed_occupany_subsample_all>0) & (predicted_occupancy_subsample_all > 0)
occupancies = observed_occupany_subsample_all[idx_to_keep]
predicted_occupancies = predicted_occupancy_subsample_all[idx_to_keep]

predicted_and_observed_occupancies = numpy.concatenate((occupancies,predicted_occupancies), axis=0)
min_ = min(predicted_and_observed_occupancies)
max_ = max(predicted_and_observed_occupancies)


x, y, z = plot_utils.get_scatter_density_arrays_for_loglog(occupancies, predicted_occupancies)
ax_gamma.scatter(x, y, c=numpy.sqrt(z), cmap='Blues', s=70, alpha=0.9, edgecolors='none', zorder=1)

ax_gamma.set_xscale('log', base=10)
ax_gamma.set_yscale('log', base=10)

ax_gamma.plot([min_,max_],[min_,max_], lw=3,ls=':',c='k',zorder=1)
ax_gamma.set_xlabel("Observed occupancy", fontsize = 9)
ax_gamma.set_ylabel("Predicted occupancy", fontsize = 9)





idx_to_keep = (observed_occupany_subsample_all>0) & (predicted_occupancy_hyper_all > 0)
occupancies = observed_occupany_subsample_all[idx_to_keep]
predicted_occupancies = predicted_occupancy_hyper_all[idx_to_keep]

predicted_and_observed_occupancies = numpy.concatenate((occupancies,predicted_occupancies), axis=0)
min_ = min(predicted_and_observed_occupancies)
max_ = max(predicted_and_observed_occupancies)


x, y, z = plot_utils.get_scatter_density_arrays_for_loglog(occupancies, predicted_occupancies)
ax_gamma_hyper.scatter(x, y, c=numpy.sqrt(z), cmap='Blues', s=70, alpha=0.9, edgecolors='none', zorder=1)
ax_gamma_hyper.set_xscale('log', base=10)
ax_gamma_hyper.set_yscale('log', base=10)
ax_gamma_hyper.plot([min_,max_],[min_,max_], lw=3,ls=':',c='k',zorder=1)


ax_gamma_hyper.set_xlabel("Observed occupancy", fontsize = 9)
ax_gamma_hyper.set_ylabel("Predicted occupancy, corrected for subsampling", fontsize = 9)



fig.subplots_adjust(hspace=0.37, wspace=0.37)
fig_name = "%stest_occupancy.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
