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


#Set parameters of the community

N=3*10**4 # number of reads
mu=-19 #mean of log(K)
#mu=-10 #mean of log(K)
#s=5**2 #variance of log(K)
#gm=0.4 #mean of sigma;
gm=0.4 #mean of sigma;

n_sites = 200




def generate_comm(mu, s, S, dist, gm, rhogamma, n_sites):
    #This function  generates the parameters K and sigma for two communities with given parameters (mu, s),
    # S species, sigma^2 distributed as 'dist' (if dist is exponential, with average gm), correlation rho between the values of K
    # and correlation rhogamma between the Gamma-distributed fluctuations of abundance
    from scipy.stats import norm
    from scipy.stats import gamma

    K = numpy.exp(numpy.random.normal(mu, s, S)) #correlated K for the two communities extracted from lognormal dist
    sigmarnd=[]  #Exponentially distributed sigma, common for the two communities
    if dist=='exp':
        for k in range(S):
            tr=100
            while tr>1.95: # Values too close to 2 give numerical problems when extracting from the Gamma distribution
                tr=numpy.sqrt(numpy.random.exponential(gm))

            sigmarnd.append(tr)
    else:
        if dist=='unif':
            sigmarnd=numpy.random.uniform(0,1.95,size=S)


    # Extraction of the who vectors of abundances, distributed according to Gamma distributions with the correlation rhogamma
    cov = numpy.ones((n_sites, n_sites))
    I = numpy.identity(n_sites)
    cov = ((cov-I)*rhogamma) + I

    Z = numpy.random.multivariate_normal(numpy.asarray([0]*n_sites), cov, S, tol=1e-6)
    U = norm.cdf(Z)

    abundances_all = []

    for idx in range(n_sites):
        G = gamma.ppf(U[:,idx], numpy.divide(2,sigmarnd)-1, scale=sigmarnd*K/2)
        abundances_all.append(G)

    abundances_all = numpy.asarray(abundances_all).T
    # Normalise, to have relative abundances
    #rel_abundances_all = abundances_all
    rel_abundances_all =  abundances_all/numpy.sum(abundances_all, axis=0)



    return rel_abundances_all




environment = 'human gut metagenome'
rarefied=True

sys.stderr.write("Subsetting samples for %s...\n" % environment)
#samples = diversity_utils.subset_observations(environment=environment)
sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
samples = sad_annotated_dict['samples']

taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))


coarse_grain_idx_dict = {}
for rank in diversity_utils.taxa_ranks:

    sys.stderr.write("Starting %s level analysis...\n" % rank)
    all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
    all_genera_idx = numpy.arange(len(all_genera))

    genus_to_taxa_dict = {}
    g_taxa_idx_all = []
    for genus in all_genera:
        genus_to_taxa_dict[genus] = []
        for t in taxa:
            if sad_annotated_dict['taxa'][t][rank] == genus:
                genus_to_taxa_dict[genus].append(t)

        g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
        g_taxa_idx_all.append(g_taxa_idx)
    counts_genera = [len(n) for n in g_taxa_idx_all]
    coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_genera))[:-1]
    coarse_grain_idx_dict[rank] = coarse_grain_by_genus_idx



#S=10**4 # number of species
S = s_by_s.shape[0]
sigma_range = numpy.logspace(0.3, numpy.log10(50), num=9, base=10.0)

fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)
# plot afd
color_all = numpy.asarray([diversity_utils.rgb_blue_taxon(r) for r in range(len(diversity_utils.taxa_ranks))])
sigma_chunk_all = [sigma_range[x:x+3] for x in range(0, len(sigma_range), 3)]
for sigma_chunk_idx, sigma_chunk in enumerate(sigma_chunk_all):
    for sigma_idx, sigma in enumerate(sigma_chunk):

        ax = plt.subplot2grid((3, 3), (sigma_chunk_idx, sigma_idx))


        #rel_s_by_s_1 = generate_comm(-13, 1.9, 3000, 'unif', gm, -0.3, n_sites)
        #rel_s_by_s_2 = generate_comm(mu, sigma, 2325, 'unif', gm, -0.3, n_sites)
        #rel_s_by_s_2 = generate_comm(mu, sigma, S, 'unif', gm, 0, n_sites)

        #rel_s_by_s = numpy.concatenate((rel_s_by_s_1, rel_s_by_s_2), axis=0)
        rel_s_by_s = generate_comm(mu, sigma, S, 'exp', gm, 0, n_sites)

        beta_div_mean = numpy.mean(rel_s_by_s, axis=1)/numpy.var(rel_s_by_s, axis=1)

        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            rel_s_by_s_coarse = numpy.add.reduceat(rel_s_by_s, coarse_grain_idx_dict[rank], axis=0)

            cv_beta_div_mean = numpy.std(beta_div_mean)/numpy.mean(beta_div_mean)

            clade_log10_rescaled_all = []
            for clade in rel_s_by_s_coarse:

                clade = clade[clade>0]
                clade_log10 = numpy.log10(clade)

                if len(clade_log10) < 4:
                    continue

                clade_log10_rescaled = (clade_log10 - numpy.mean(clade_log10))/numpy.std(clade_log10)
                clade_log10_rescaled_all.append(clade_log10_rescaled)

            clade_log10_rescaled_all_flat = numpy.concatenate(clade_log10_rescaled_all).ravel()
            hist_, bin_edges_ = numpy.histogram(clade_log10_rescaled_all_flat, density=True, bins=20)
            bins_mean_ = numpy.asarray([0.5 * (bin_edges_[i] + bin_edges_[i+1]) for i in range(0, len(bin_edges_)-1 )])
            hist_to_plot = hist_[hist_>0]
            bins_mean_to_plot = bins_mean_[hist_>0]

            ax.scatter(bins_mean_to_plot, hist_to_plot, s=10, color=diversity_utils.rgb_blue_taxon(rank_idx), alpha=0.9, lw=2, label=rank)


        ax.set_title('sigma = %s, CV_{beta/mean} = %s' % (str(round(sigma, 3)), str(round(cv_beta_div_mean, 3))), fontsize=9)

        ax.set_yscale('log', base=10)




fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%ssimulate_afd_taxon_emp%s.png" % (config.analysis_directory, diversity_utils.get_rarefied_label(rarefied))
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
