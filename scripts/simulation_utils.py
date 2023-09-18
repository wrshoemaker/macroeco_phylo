import os
#from Bio import Phylo

import numpy


from scipy.special import polygamma
from scipy.stats import norm
from scipy.stats import gamma
#from scipy.linalg import cholesky


numpy.random.seed(123456789)

#Set parameters of the community
#N=3*10**4 # number of reads
#mu=-19 #mean of log(K)
#mu=-10 #mean of log(K)
#s=5**2 #variance of log(K)
##gm=0.4 #mean of sigma;
#gm=0.4 #mean of sigma;
#S=10**4
#n_sites = 200


def vectorized(prob_matrix, items):
    s = prob_matrix.cumsum(axis=1)
    r = numpy.random.rand(prob_matrix.shape[0])
    print(r.shape, s.shape)
    print(s < r)
    k = (s < r).sum(axis=1)
    return items[k]



def sample_constrained_on_richness(N, rel_abundances):

    asv_idx = []

    for sad_idx, sad in enumerate(rel_abundances_all.T):
        read_counts_all.append(numpy.random.multinomial(int(N[sad_idx]), sad))




def generate_community(mu, s, S, N, dist, gm, n_sites, rhogamma=0):
    #This function  generates the parameters K and sigma for two communities with given parameters (mu, s),
    # S species, sigma^2 distributed as 'dist' (if dist is exponential, with average gm), correlation rho between the values of K
    # and correlation rhogamma between the Gamma-distributed fluctuations of abundance

    if len(N) == 1:
        N = numpy.asarray([N]*n_sites)

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

        if dist == 'constant':
            sigmarnd=numpy.repeat(gm, S)


    # Extraction of the who vectors of abundances, distributed according to Gamma distributions with the correlation rhogamma
    cov = numpy.ones((n_sites, n_sites))
    I = numpy.identity(n_sites)
    cov = ((cov-I)*rhogamma) + I

    Z = numpy.random.multivariate_normal(numpy.asarray([0]*n_sites), cov, S, tol=1e-5)
    U = norm.cdf(Z)

    abundances_all = []
    for idx in range(n_sites):
        G = gamma.ppf(U[:,idx], numpy.divide(2,sigmarnd)-1, scale=sigmarnd*K/2)
        abundances_all.append(G)

    abundances_all = numpy.asarray(abundances_all).T
    # Normalise, to have relative abundances
    #rel_abundances_all = abundances_all
    rel_abundances_all =  abundances_all/numpy.sum(abundances_all, axis=0)

    # run multinomial
    #vals = numpy.asarray([N]*n_sites,dtype=int)
    #read_counts_all = vectorized(rel_abundances_all, vals)

    #read_counts_all = (numpy.random.multinomial(1,rel_abundances_all,size=N)==1).argmax(1)
    read_counts_multinomial_all = []
    read_counts_poisson_all = []
    for sad_idx, sad in enumerate(rel_abundances_all.T):
        read_counts_multinomial_all.append(numpy.random.multinomial(int(N[sad_idx]), sad))
        read_counts_poisson_all.append(numpy.random.poisson(lam=int(N[sad_idx])*rel_abundances_all))

    read_counts_multinomial_all = numpy.asarray(read_counts_multinomial_all).T
    read_counts_poisson_all = numpy.asarray(read_counts_poisson_all).T


    return rel_abundances_all, read_counts_multinomial_all







def genrate_community_from_mean_and_var_random_sites(mean, var, N, n_sites, rhogamma=0, corr_matrix=None):

    if corr_matrix is not None:
        cov = corr_matrix

    else:
        # Extraction of the who vectors of abundances, distributed according to Gamma distributions with the correlation rhogamma
        cov = numpy.ones((n_sites, n_sites))
        I = numpy.identity(n_sites)
        cov = ((cov-I)*rhogamma) + I

    Z = numpy.random.multivariate_normal(numpy.asarray([0]*n_sites), cov, len(mean), tol=1e-5)
    U = norm.cdf(Z)

    beta = (mean**2)/var

    abundances_all = []
    for idx in range(n_sites):
        # scale = 1/rate
        G = gamma.ppf(U[:,idx], beta, scale=mean/beta)
        abundances_all.append(G)

    abundances_all = numpy.asarray(abundances_all).T

    # Normalise, to have relative abundances
    rel_abundances_all =  abundances_all/numpy.sum(abundances_all, axis=0)

    # run multinomial
    #vals = numpy.asarray([N]*n_sites,dtype=int)
    #read_counts_all = vectorized(rel_abundances_all, vals)

    #read_counts_all = (numpy.random.multinomial(1,rel_abundances_all,size=N)==1).argmax(1)

    read_counts_all = []
    for sad_idx, sad in enumerate(rel_abundances_all.T):
        read_counts_all.append(numpy.random.multinomial(int(N[sad_idx]), sad))

    read_counts_all = numpy.asarray(read_counts_all).T

    return rel_abundances_all, read_counts_all




#import numba as nb

#@nb.njit('float64[:, :](float64[:, :])')
#def cholesky_numba(A):
#    n = A.shape[0]
#    L = numpy.zeros_like(A)
#    for i in range(n):
#        for j in range(i+1):
#            s = 0
#            for k in range(j):
#                s += L[i][k] * L[j][k]

#            if (i == j):
#                L[i][j] = (A[i][i] - s) ** 0.5
#            else:
#                L[i][j] = (1.0 / L[j][j] * (A[i][j] - s))
#    return L


def genrate_community_from_mean_and_var_rvs(mean, var, N, n_sites):

    beta = (mean**2)/var

    s_by_s = []

    for i in range(len(mean)):
        
        mean_i = mean[i]
        beta_i = beta[i]

        afd_i = gamma.rvs(beta_i, scale=mean_i/beta_i, size=n_sites)
        s_by_s.append(afd_i)

    rel_abundances_all = numpy.asarray(s_by_s)
    rel_abundances_all = rel_abundances_all/numpy.sum(rel_abundances_all, axis=0)

    read_counts_all = []
    for sad_idx, sad in enumerate(rel_abundances_all.T):
        read_counts_all.append(numpy.random.multinomial(int(N[sad_idx]), sad))

    read_counts_all = numpy.asarray(read_counts_all).T

    return rel_abundances_all, read_counts_all




def genrate_community_from_mean_and_var(mean, var, N, n_sites, rhogamma=0, corr_matrix=None):

    if corr_matrix is not None:
        cov = corr_matrix

    else:
        # Extraction of the who vectors of abundances, distributed according to Gamma distributions with the correlation rhogamma
        #cov = numpy.ones((len(mean), len(mean)))
        I = numpy.identity(len(mean))
        #cov = ((cov-I)*rhogamma) + I
        cov = I

    Z = numpy.random.multivariate_normal(numpy.asarray([0]*len(mean)), cov, n_sites, tol=1e-6)
    
    #l = cholesky(cov, check_finite=False, overwrite_a=True)
    #Z = numpy.asarray([0]*len(mean)) + l.dot(numpy.random.standard_normal(len(means)))

    #Z = numpy.asarray([0]*len(mean)) + numpy.linalg.cholesky(cov) @ numpy.random.standard_normal(len(mean))
    #Z = numpy.linalg.cholesky(cov) @ numpy.random.standard_normal(len(mean))

    U = norm.cdf(Z)

    beta = (mean**2)/var

    abundances_all = []
    for idx in range(n_sites):
        # scale = 1/rate
        #G = gamma.ppf(U[:,idx], beta, scale=mean/beta)
        G = gamma.ppf(U[idx,:], beta, scale=mean/beta)
        abundances_all.append(G)

    abundances_all = numpy.asarray(abundances_all).T

    # Normalise, to have relative abundances
    #rel_abundances_all = abundances_all
    rel_abundances_all =  abundances_all/numpy.sum(abundances_all, axis=0)

    # run multinomial
    #vals = numpy.asarray([N]*n_sites,dtype=int)
    #read_counts_all = vectorized(rel_abundances_all, vals)

    #read_counts_all = (numpy.random.multinomial(1,rel_abundances_all,size=N)==1).argmax(1)

    read_counts_all = []
    for sad_idx, sad in enumerate(rel_abundances_all.T):
        read_counts_all.append(numpy.random.multinomial(int(N[sad_idx]), sad))

    read_counts_all = numpy.asarray(read_counts_all).T

    return rel_abundances_all, read_counts_all





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

    # Normalise, to have relative abundances
    rel_abundances_all =  abundances_all/numpy.sum(abundances_all, axis=0)



    read_counts_multinomial_all = []
    #read_counts_poisson_all = []
    for sad_idx, sad in enumerate(rel_abundances_all.T):
        read_counts_multinomial_all.append(numpy.random.multinomial(int(N[sad_idx]), sad))
        #read_counts_poisson_all.append(numpy.random.poisson(lam=int(N[sad_idx])*rel_abundances_all[:,sad_idx]))

    read_counts_multinomial_all = numpy.asarray(read_counts_multinomial_all).T
    #read_counts_poisson_all = numpy.asarray(read_counts_poisson_all).T

    #, read_counts_poisson_all
    

    return rel_abundances_all, read_counts_multinomial_all



def genrate_community_from_mean_and_var_svd_multinomial_and_poisson(mean, var, N, n_sites, svd):

    Z = numpy_multivariate_normal(numpy.asarray([0]*len(mean)), svd, size=n_sites)

    U = norm.cdf(Z)

    beta = (mean**2)/var

    abundances_all = []
    for idx in range(n_sites):
        # scale = 1/rate
        G = gamma.ppf(U[idx,:], beta, scale=mean/beta)
        abundances_all.append(G)

    abundances_all = numpy.asarray(abundances_all).T

    # Normalise, to have relative abundances
    rel_abundances_all =  abundances_all/numpy.sum(abundances_all, axis=0)



    read_counts_multinomial_all = []
    read_counts_poisson_all = []
    for sad_idx, sad in enumerate(rel_abundances_all.T):
        read_counts_multinomial_all.append(numpy.random.multinomial(int(N[sad_idx]), sad))
        read_counts_poisson_all.append(numpy.random.poisson(lam=int(N[sad_idx])*rel_abundances_all[:,sad_idx]))

    read_counts_multinomial_all = numpy.asarray(read_counts_multinomial_all).T
    read_counts_poisson_all = numpy.asarray(read_counts_poisson_all).T

    #, read_counts_poisson_all
    

    return rel_abundances_all, read_counts_multinomial_all, read_counts_poisson_all









def genrate_community_from_mean_and_var_using_rvs(mean, var, N, n_sites):

    beta = (mean**2)/var

    abundances_all = []
    for idx in range(len(mean)):
        
        afd = gamma.rvs(beta[idx], scale=mean[idx]/beta[idx], size=n_sites)
        abundances_all.append(afd)

    abundances_all = numpy.asarray(abundances_all)

    # make sure relative abundances sum to 1
    rel_abundances_all =  abundances_all/numpy.sum(abundances_all, axis=0)

    read_counts_all = []
    for sad_idx, sad in enumerate(rel_abundances_all.T):
        read_counts_all.append(numpy.random.multinomial(int(N[sad_idx]), sad))

    read_counts_all = numpy.asarray(read_counts_all).T

    return rel_abundances_all, read_counts_all



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







def predict_shannon_diversity(mu, s, sigma, S, N):
    return -1 * (S/N) * (1 - (sigma/2))* ((numpy.exp(mu + (s**2)/2) * (mu + s**2 - numpy.log(N))) +  numpy.log(sigma/2) + polygamma(0, 2/sigma) + mu)


#def predict_coarse_shannon_diversity(mu, s, sigma, S, N, g):
#    return -1 * (S/N) * (1 - (sigma/2))* ((numpy.exp(mu + (s**2)/2) * (mu + s**2 - numpy.log(N) +  numpy.log(g))) +  numpy.log(sigma/2) + polygamma(0, 2/sigma) + mu)


def predict_coarse_shannon_diversity(mu, s, sigma, N, g):
    # g if vector
    return -1 * (1/N) * (1 - (sigma/2))* sum(((numpy.exp(mu + (s**2)/2) * (mu + s**2 - numpy.log(N) +  numpy.log(g))) +  numpy.log(sigma/2) + polygamma(0, 2/sigma) + mu))
