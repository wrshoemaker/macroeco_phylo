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
import scipy.integrate as integrate



# remove means *greater* than this value
max_mad = [1e-0, 1e-1, 1e-2, 1e-3]



def predict_mean_and_var_richness_and_diversity_using_gamma_rv(s_by_s, iter_=100):

    rel_s_by_s_np = (s_by_s/s_by_s.sum(axis=0))

    n_reads = s_by_s.sum(axis=0)

    mean_rel_s_by_s = numpy.mean(rel_s_by_s_np, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s_np, axis=1)

    mean_richness_all = []
    var_richness_all = []
    mean_diversity_all = []
    var_diversity_all = []

    # make SVD of identity matrix

    u = numpy.identity(len(mean_rel_s_by_s))
    s = numpy.asarray([1]*len(mean_rel_s_by_s))
    v = numpy.identity(len(mean_rel_s_by_s))
    svd = (u, s, v)

    var_dict = {}

    for max_mad_j in max_mad:
        var_dict[max_mad_j] = []

    for i in range(iter_):

        rel_abundances_gamma, read_counts_gamma = simulation_utils.genrate_community_from_mean_and_var_svd(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads), svd)

        read_counts_gamma_pres_abs = (read_counts_gamma>0)

        richness = read_counts_gamma_pres_abs.sum(axis=0)

        mean_richness_all.append(numpy.mean(richness))
        var_richness_all.append(numpy.var(richness))

        # diversity
        diversity_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, read_counts_gamma)
       
        mean_diversity_all.append(numpy.mean(diversity_null))
        var_diversity_all.append(numpy.var(diversity_null))

        for max_mad_j in max_mad:

            read_counts_gamma_pres_abs_max_mad_j = read_counts_gamma_pres_abs[(mean_rel_s_by_s<=max_mad_j),:]

            diversity_max_mad_j = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, read_counts_gamma_pres_abs_max_mad_j)
            var_diversity_max_mad_j = numpy.var(diversity_max_mad_j)

            var_dict[max_mad_j].append(var_diversity_max_mad_j)


    return var_dict



def predict_mean_and_var_diversity_analytic(s_by_s, species):

    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

    n_reads = s_by_s.sum(axis=0)

    mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
    beta_rel_s_by_s = (mean_rel_s_by_s**2)/var_rel_s_by_s

    diversity_first_moment_all = []
    diversity_second_moment_all = []


    var_dict = {}
    for max_mad_j in max_mad:
        var_dict[max_mad_j] = []


    # dict with integral for each species for each sample
    for m in range(len(n_reads)):

        N_m = int(n_reads[m])

        diversity_first_moment_m = 0
        diversity_second_moment_m = 0


        integrand_first_moment_all = []
        integrand_second_moment_first_term_all = []
        for i in range(len(mean_rel_s_by_s)):
            
            mean_i = mean_rel_s_by_s[i]
            beta_i = beta_rel_s_by_s[i]

            integral_first_moment_result = integrate.quad(diversity_utils.integrand_first_moment, 0, N_m, args=(N_m, mean_i, beta_i), epsabs=1e-20)
            #diversity_first_moment_m += integral_first_moment_result[0]
            integrand_first_moment_all.append(integral_first_moment_result[0])

            integrand_second_moment_result = integrate.quad(diversity_utils.integrand_second_moment, 0, N_m, args=(N_m, mean_i, beta_i),  epsabs=1e-25)
            integrand_second_moment_first_term_all.append(integrand_second_moment_result[0])


        integrand_first_moment_all = numpy.asarray(integrand_first_moment_all)
        integrand_second_moment_first_term_all = numpy.asarray(integrand_second_moment_first_term_all)

        for max_mad_j in max_mad:

            idx_to_keep = (mean_rel_s_by_s<=max_mad_j)

            integrand_first_moment_all_max_mad_j = integrand_first_moment_all[idx_to_keep]
            integrand_second_moment_first_term_all_max_mad_j = integrand_second_moment_first_term_all[idx_to_keep]

            diversity_second_moment_second_term_m = 0
            for i in range(len(integrand_first_moment_all_max_mad_j)):
                for j in range(i):
                    diversity_second_moment_second_term_m += integrand_first_moment_all_max_mad_j[i]*integrand_first_moment_all_max_mad_j[j]


            var_diversity_max_mad_j = sum(integrand_second_moment_first_term_all_max_mad_j)+(2*diversity_second_moment_second_term_m) - (sum(integrand_first_moment_all_max_mad_j)**2)
            var_dict[max_mad_j].append(var_diversity_max_mad_j)

    return var_dict











environment = 'marine metagenome'

sys.stderr.write("Subsetting samples for %s...\n" % environment)
sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = False)

samples = sad_annotated_dict['samples']
taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
sys.stderr.write("Getting site-by-species matrix...\n")

s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])


rank = 'genus'
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





var_dict_pred = predict_mean_and_var_diversity_analytic(s_by_s_genera, s_by_s_genera.shape[0])
var_dict_sim = predict_mean_and_var_richness_and_diversity_using_gamma_rv(s_by_s_genera)

for m in max_mad:

    mean_var_pred = numpy.mean(var_dict_pred[m])
    mean_var_sim = numpy.mean(var_dict_sim[m])

    rel_error = numpy.absolute(mean_var_pred - mean_var_sim)/mean_var_sim

    print(mean_var_pred, mean_var_sim, rel_error)



#mean_richness_observed_gamma_rvs, var_richness_observed_gamma_rvs, mean_richness_predicted_gamma_rvs, var_richness_predicted_gamma_rvs, mean_diversity_observed_gamma_rvs, var_diversity_observed_gamma_rvs, mean_diversity_predicted_gamma_rvs, var_diversity_predicted_gamma_rvs = diversity_utils.predict_mean_and_var_richness_and_diversity_using_gamma_rv(s_by_s_genera, iter_=iter_)


#mean_diversity_observed, mean_diversity_predicted, var_diversity_observed, var_diversity_predicted, var_diversity_predicted_plus_covariance = diversity_utils.predict_mean_and_var_diversity_analytic(s_by_s_genera, s_by_s_genera.shape[0])

#mean_richness_observed, mean_richness_predicted = diversity_utils.predict_mean_richness(s_by_s_genera, s_by_s_genera.shape[0])
#var_richness_observed, var_richness_predicted, var_richness_predicted_plus_covariance = diversity_utils.predict_var_richness(s_by_s_genera, s_by_s_genera.shape[0])


