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

from scipy.special import polygamma
from scipy.stats import norm
from scipy.stats import gamma
#import dbd_utils
#import plot_utils

#import diversity_utils
import config
import scipy.integrate as integrate
import scipy.special as special


# replace path with your own
s_by_s_path = config.data_directory + "s_by_s_test.pickle"






def prob_n_reads(n, N, mean_, beta_):

    # exp( gammaln(beta+n) - gammaln(n+1) - gammaln(beta) )
    # gamma of factorial results in numerical overflow, do logamma trick instead
    # gamma(beta+n) and gamma(n+1) are large, but their ratio is not, so gammaln(beta+n) - gammaln(n+1) is ok and can be exponentiated

    return numpy.exp( special.gammaln(beta_+n) - special.gammaln(n+1) - special.gammaln(beta_) )   * (((mean_*N)/(beta_ + mean_*N))**n) * ((beta_/(beta_ + mean_*N))**beta_)


def integrand_first_moment(n, N, mean_, beta_):
    return (n/N)*numpy.log(n/N) * prob_n_reads(n, N, mean_, beta_)


def integrand_second_moment(n, N, mean_, beta_):
    return (((n/N)*numpy.log(n/N))**2) * prob_n_reads(n, N, mean_, beta_)






# build and save genus-level site by family read count matrix for marine environment as pickle
def make_test_s_by_s():

    import dbd_utils

    environment =  'marine metagenome'

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = False)

    samples = sad_annotated_dict['samples']
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    sys.stderr.write("Getting site-by-species matrix...\n")

    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    taxa_ranks = ['family']
    #diversity_utils.taxa_ranks

    fig, ax = plt.subplots(figsize=(4,4))

    for rank_idx, rank in enumerate(taxa_ranks):

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
        




        s_by_s_genera_list = s_by_s_genera.tolist()


        dict_ = {}
        dict_['s_by_s'] = s_by_s_genera_list
        
        # protocol 2
        sys.stderr.write("Saving prediction dictionary...\n")
        with open(s_by_s_path, 'wb') as handle:
            pickle.dump(dict_, handle, protocol=2)






def load_test_s_by_s():


    with open(s_by_s_path, 'rb') as handle:
        distance_collapsed_dict = pickle.load(handle)
    return distance_collapsed_dict






def numpy_multivariate_normal(mean, svd, size=None, tol=1e-6):

    if size is None:
        shape = []
    elif isinstance(size, (int, numpy.integer)):
        shape = [size]
    else:
        shape = size

    # Compute shape of output and create a matrix of independent
    # standard normally distributed random numbers. The matrix has rows
    # with the same length as mean and as many rows are necessary to
    # form a matrix of shape final_shape.
    final_shape = list(shape[:])
    final_shape.append(mean.shape[0])
    x = numpy.random.standard_normal(final_shape).reshape(-1, mean.shape[0])

    # Transform matrix of standard normals into matrix where each row
    # contains multivariate normals with the desired covariance.
    # Compute A such that dot(transpose(A),A) == cov.
    # Then the matrix products of the rows of x and A has the desired
    # covariance. Note that sqrt(s)*v where (u,s,v) is the singular value
    # decomposition of cov is such an A.
    #
    # Also check that cov is positive-semidefinite. If so, the u.T and v
    # matrices should be equal up to roundoff error if cov is
    # symmetric and the singular value of the corresponding row is
    # not zero. We continue to use the SVD rather than Cholesky in
    # order to preserve current outputs. Note that symmetry has not
    # been checked.

    # GH10839, ensure double to make tol meaningful
    #cov = cov.astype(numpy.double)
    #(u, s, v) = svd(cov)

    u, s, v = svd

    x = numpy.dot(x, numpy.sqrt(s)[:, None] * v)
    x += mean
    x.shape = tuple(final_shape)
    return x



def genrate_community_from_mean_and_var_svd(mean, var, N, n_sites, svd):

    Z = numpy_multivariate_normal(numpy.asarray([0]*len(mean)), svd, size=n_sites)

    U = norm.cdf(Z)

    beta = (mean**2)/var

    abundances_all = []
    for idx in range(n_sites):
        # scale = 1/rate
        G = gamma.ppf(U[idx,:], beta, scale=mean/beta)
        abundances_all.append(G)

    abundances_all = numpy.asarray(abundances_all).T

    # Normalise, to have relative abundances sum to one
    rel_abundances_all =  abundances_all/numpy.sum(abundances_all, axis=0)

    read_counts_all = []
    for sad_idx, sad in enumerate(rel_abundances_all.T):
        read_counts_all.append(numpy.random.multinomial(int(N[sad_idx]), sad))

    read_counts_all = numpy.asarray(read_counts_all).T

    return rel_abundances_all, read_counts_all





def predict_second_moment(s_by_s, iter_=100):

    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

    n_reads = s_by_s.sum(axis=0)

    mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
    beta_rel_s_by_s = (mean_rel_s_by_s**2)/var_rel_s_by_s

    diversity_second_moment_all = []
    # dict with integral for each species for each sample
    product_pairs = []
    for m in range(len(n_reads)):

        N_m = n_reads[m]

        diversity_second_moment_m = []
        diversity_second_moment_second_term_m = []

        product_pairs_m = []
        
        for i in range(len(mean_rel_s_by_s)):

            mean_i = mean_rel_s_by_s[i]
            beta_i = beta_rel_s_by_s[i]

            integrand_second_moment_result = integrate.quad(integrand_second_moment, 0, N_m, args=(N_m, mean_i, beta_i),  epsabs=1e-25)
            diversity_second_moment_m.append(numpy.absolute(integrand_second_moment_result[0]))

            integrand_first_moment_result = integrate.quad(integrand_first_moment, 0, N_m, args=(N_m, mean_i, beta_i),  epsabs=1e-25)
            diversity_second_moment_second_term_m.append(integrand_first_moment_result[0])

        diversity_second_moment_all.append(diversity_second_moment_m)

        for i in range(len(mean_rel_s_by_s)):
            for j in range(i):

                product_pairs_m.append(diversity_second_moment_second_term_m[i]*diversity_second_moment_second_term_m[j])

        product_pairs.append(product_pairs_m)


    product_pairs = numpy.asarray(product_pairs)
    expected_product_pairs = numpy.mean(product_pairs, axis=0)

    diversity_second_moment_all = numpy.asarray(diversity_second_moment_all)
    expected_diversity_second_moment = numpy.mean(diversity_second_moment_all, axis=0)

    u = numpy.identity(len(mean_rel_s_by_s))
    s = numpy.asarray([1]*len(mean_rel_s_by_s))
    v = numpy.identity(len(mean_rel_s_by_s))
    svd = (u, s, v)

    expected_diversity_second_moment_rvs_all = []
    product_pairs_null = []
    for i in range(iter_):
        
        rel_abundances_gamma, read_counts_gamma = genrate_community_from_mean_and_var_svd(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads), svd)
        rel_read_counts_gamma = read_counts_gamma/numpy.sum(read_counts_gamma, axis=0)

        measure_rel_read_counts_gamma = (rel_read_counts_gamma*numpy.log(rel_read_counts_gamma))
        measure_rel_read_counts_gamma[numpy.isnan(measure_rel_read_counts_gamma)] = 0

        expected_diversity_first_moment_rvs = numpy.mean(measure_rel_read_counts_gamma, axis=1)
        
        rel_read_counts_gamma_second_moment = measure_rel_read_counts_gamma**2
        rel_read_counts_gamma_second_moment[numpy.isnan(rel_read_counts_gamma_second_moment)] = 0

        expected_diversity_second_moment_rvs = numpy.mean(rel_read_counts_gamma_second_moment, axis=1)
        expected_diversity_second_moment_rvs_all.append(expected_diversity_second_moment_rvs)

        product_pairs_null_i = []
        # calculates all product pairs of the first moment < x_i log x_i > * < x_j log x_j >
        for i in range(len(expected_diversity_first_moment_rvs)):
            for j in range(i):

                product_pairs_null_i.append(expected_diversity_first_moment_rvs[i] * expected_diversity_first_moment_rvs[j])

        product_pairs_null.append(product_pairs_null_i)


    # take the mean of < x_i log x_i > * < x_j log x_j > over the simulations
    product_pairs_null = numpy.asarray(product_pairs_null)
    expected_product_pairs_null_rvs = numpy.mean(product_pairs_null, axis=0)

    
    expected_diversity_second_moment_rvs_all = numpy.asarray(expected_diversity_second_moment_rvs_all)
    expected_diversity_second_moment_rvs = numpy.mean(expected_diversity_second_moment_rvs_all, axis=0)
    

    return expected_diversity_second_moment, expected_diversity_second_moment_rvs, expected_product_pairs, expected_product_pairs_null_rvs





def predict_product_first_moments():

    s_by_s = load_test_s_by_s()
    s_by_s = numpy.asarray(s_by_s['s_by_s'])

    expected_diversity_second_moment_first_term, expected_diversity_second_moment_first_term_rvs, expected_diversity_second_moment_second_term, expected_diversity_second_moment_second_term_rvs = predict_second_moment(s_by_s)

    # test variance
    expected_diversity_second_moment = expected_diversity_second_moment_second_term
    expected_diversity_second_moment_rvs = expected_diversity_second_moment_second_term_rvs

    idx_to_keep = (expected_diversity_second_moment>0) & (expected_diversity_second_moment_rvs>0)
    expected_diversity_second_moment = expected_diversity_second_moment[idx_to_keep]
    expected_diversity_second_moment_rvs = expected_diversity_second_moment_rvs[idx_to_keep]

    min_ = min(numpy.union1d(expected_diversity_second_moment, expected_diversity_second_moment_rvs))
    max_ = max(numpy.union1d(expected_diversity_second_moment, expected_diversity_second_moment_rvs))

    fig, ax = plt.subplots(figsize=(4,4))

    ax.plot([min_*0.5,max_*1.1],[min_*0.5,max_*1.1], lw=2,ls='--',c='k',zorder=2, label='1:1')
    ax.set_xlim([min_*0.5,max_*1.1])
    ax.set_ylim([min_*0.5,max_*1.1])
    ax.set_xlabel("Predicted second moment, analytic", fontsize = 10)
    ax.set_ylabel("Predicted second momen, simulation", fontsize = 10)

    ax.scatter(expected_diversity_second_moment, expected_diversity_second_moment_rvs, s=7, color='k', alpha=0.1, zorder=2)

    ax.set_xscale('log', base=10)
    ax.set_yscale('log', base=10)

    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%ssecond_moment_test.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()



#make_test_s_by_s()

predict_product_first_moments()


