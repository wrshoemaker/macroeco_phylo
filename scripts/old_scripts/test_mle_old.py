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
import mle_utils



rarefied = True
iter = 1000

#fig = plt.figure(figsize = (10, 10)) #
#fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

#environments_to_keep = diversity_utils.environments_to_keep
environment = 'human gut metagenome'

taxa_ranks = diversity_utils.taxa_ranks
idx_taxa_ranks = range(len(taxa_ranks))
taxa_ranks_label = diversity_utils.taxa_ranks_label


sys.stderr.write("Subsetting samples for %s...\n" % environment)
#samples = diversity_utils.subset_observations(environment=environment)
sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
samples = sad_annotated_dict['samples']

taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))


#print(taxa_ranks)

sys.stderr.write("Running taxonomic coarse-graining.....\n")
for rank_idx, rank in enumerate(['phylum']):

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

        g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])

        #rel_s_by_s_genus = rel_s_by_s[g_taxa_idx,:]
        s_by_s_genus = s_by_s[g_taxa_idx,:]
        sad_genera_all.append(s_by_s_genus.sum(axis=0))



    s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
    # remove sites where there are no observations
    s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]
    #rel_s_by_s_genera = (s_by_s_genera/s_by_s_genera.sum(axis=0))

    N_all = s_by_s_genera.sum(axis=0)

    for afd in s_by_s_genera:

        #if sum(afd>0)/len(afd) < 0.2:
        #    continue


        x = afd/N_all
        mu_start = numpy.mean(x)
        sigma_start = numpy.std(x)
        beta_start = (mu_start**2)/(sigma_start**2)

        #max_mu, max_sigma, min_mu, min_sigma = mle_utils.get_hyp1f1_bounds(afd, N_all)


        #trunc_gaussian_sampling_model = mle_utils.mle_truncated_gaussian_sampling(N_all, afd)
        #trunc_gaussian_sampling_result = trunc_gaussian_sampling_model.fit(method="lbfgs", disp = False, bounds=[(min_mu, max_mu), (min_sigma, max_sigma)])
        #trunc_gaussian_sampling_model_ll = trunc_gaussian_sampling_model.loglike(trunc_gaussian_sampling_result.params)

        gamma_sampling_model = mle_utils.mle_gamma_sampling(N_all, afd)
        gamma_sampling_result = gamma_sampling_model.fit(method="lbfgs", disp = False, bounds=[ (0.0001, 10000), (min(afd/N_all), max(afd/N_all))], start_params=numpy.array([beta_start, mu_start]))
        gamma_sampling_model_ll = gamma_sampling_model.loglike(gamma_sampling_result.params)


        N_all_no_0 = N_all[afd>0]
        afd_no_0 = afd[afd>0]
        #trunc_gaussian_model = mle_utils.mle_truncated_gaussian(N_all_no_0, afd_no_0)
        #trunc_gaussian_result = trunc_gaussian_model.fit(method="lbfgs", disp = False, bounds=[(min(afd_no_0/N_all_no_0), max(afd_no_0/N_all_no_0)), (0.0001, 10000)])
        #trunc_gaussian_model_ll = trunc_gaussian_model.loglike(trunc_gaussian_result.params)

        #gamma_model = mle_utils.mle_gamma(N_all_no_0, afd_no_0)
        #gamma_result = gamma_model.fit(method="lbfgs", disp = False, bounds= [(0.0000001,10000000), (0.000001,1)])
        #gamma_model_ll = gamma_model.loglike(gamma_result.params)

        trunc_gamma_model = mle_utils.mle_truncated_gamma(N_all_no_0, afd_no_0)
        trunc_gamma_result = trunc_gamma_model.fit(method="lbfgs", disp = False, bounds= [(0.0001,1000), (0.0001,1)], start_params=numpy.array([beta_start, mu_start]))
        trunc_gamma_model_ll = trunc_gamma_model.loglike(trunc_gamma_result.params)



        #print(trunc_gaussian_model_ll, trunc_gamma_model_ll)

        #print(sum(afd>0)/len(afd), gamma_sampling_result.params[0], trunc_gamma_result.params[0])

        print(gamma_sampling_result.params[0])

        #ll_ratio = -2*(trunc_gaussian_model_ll - gamma_model_ll)

        #print(ll_ratio)



        #gamma_sampling_model = mle_utils.mle_gamma_sampling(N_all, afd)
        #gamma_sampling_result = gamma_sampling_model.fit(method="lbfgs", disp = False, bounds= [(0.0000001,10000000), (0.000001,1)])
        #gamma_sampling_model_ll = gamma_sampling_model.loglike(gamma_sampling_result.params)

        #print(trunc_gaussian_sampling_model_ll, gamma_sampling_model_ll)
