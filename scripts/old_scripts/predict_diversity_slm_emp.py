import os
#from Bio import Phylo
import random
import copy
import sys
import numpy
import random
import pickle
import scipy.stats as stats
import datetime
import signal
import threading
import time
import mle_utils


import dbd_utils

import diversity_utils
import config

from scipy.special import erf
from scipy.special import rel_entr


import simulation_utils
from sympy.core.cache import clear_cache


numpy.random.seed(123456789)


slm_emp_taxon_parameter_dict_template = config.data_directory + "slm_emp_taxon_parameter_dict/slm_emp_parameter_dict_%s%s.pickle"
slm_emp_phylo_parameter_dict_template = config.data_directory + "slm_emp_phylo_parameter_dict/slm_emp_parameter_dict_%s%s.pickle"

slm_emp_taxon_diversity_dict_template = config.data_directory + "slm_emp_taxon_diversity_dict/diversity_slm_%s.pickle"
slm_emp_phylo_diversity_dict_template = config.data_directory + "slm_emp_phylo_diversity_dict/diversity_slm_%s.pickle"


slm_gamma_emp_taxon_diversity_dict_template = config.data_directory + "slm_gamma_emp_taxon_diversity_dict/diversity_slm_%s.pickle"
slm_gamma_emp_phylo_diversity_dict_template = config.data_directory + "slm_gama_emp_phylo_diversity_dict/diversity_slm_%s.pickle"


slm_gamma_emp_diversity_dict_template = config.data_directory + "slm_gamma_emp_diversity_dict/diversity_slm_%s.pickle"


class TimeoutException(Exception):   # Custom exception class
    pass

def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException


class RunWithTimeout(object):
    def __init__(self, function, args):
        self.function = function
        self.args = args
        self.answer = None

    def worker(self):
        self.answer = self.function(*self.args)

    def run(self, timeout):
        thread = threading.Thread(target=self.worker)
        thread.start()
        thread.join(timeout)
        return self.answer




#def subsample_fixed_richness(s_by_s, N, S_obs):

#    s_by_s_sampled = numpy.zeros((S_obs, len(N)))

#    print(s_by_s.shape)
#    print(s_by_s_sampled.shape)

#    N_copy = numpy.copy(N)

#    sample_dict = {}

#    while len(sample_dict) < S_obs:

#        for N_i_idx, N_i in enumerate(N):

#            if N_i == 0:
#                continue

#            species = numpy.random.choice(s_by_s.shape[0], size=1, replace=False, p=s_by_s[:,N_i_idx])[0]
#            if species not in sample_dict:
#                sample_dict[species] = {}

#            if N_i_idx not in sample_dict[species]:
#                sample_dict[species][N_i_idx] = 0

#            sample_dict[species][N_i_idx] += 1
#            N_copy[N_i_idx] -= 1

#        print(len(sample_dict), S_obs)










def make_all_slm_parameter_dict(environment, rarefied=False, taxon=True, iter=100, fix_richness=True):

    sys.stderr.write("Subsetting samples for %s...\n" % environment)

    if fix_richness == True:
        fix_richness_label = '_fix_richness'
    else:
        fix_richness_label = ''


    if taxon == True:
        sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
        path_ = slm_emp_taxon_parameter_dict_template % (diversity_utils.get_label(environment, rarefied), fix_richness_label)

    else:
        sad_annotated_dict = dbd_utils.load_sad_annotated_phylo_dict(environment, rarefied = rarefied)
        path_ = slm_emp_phylo_parameter_dict_template % (diversity_utils.get_label(environment, rarefied), fix_richness_label)


    samples = sad_annotated_dict['samples']
    sys.stderr.write("Getting site-by-species matrix...\n")
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    N = s_by_s.sum(axis=0)
    S_obs = s_by_s.shape[0]

    diversity_species = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s)

    #mu=-16
    #s=3
    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
    mean = numpy.mean(rel_s_by_s, axis=1)

    #min_c = 1/max(s_by_s.sum(axis=0))
    min_c = 0.0000001
    # pad max_c because high c values fail to converge
    max_c = (1/min(s_by_s.sum(axis=0)))*100
    c_to_keep = []
    sigma_to_keep = []
    euc_distance_to_keep = []
    S_simulation_to_keep = []

    #count_greater = 0
    #count_less = 0
    count = 0
    #signal.signal(signal.SIGALRM, timeout_handler)
    #while (count_less < iter) or (count_greater < iter):
    while (count < iter):

        #start = time.time()
        c_i = numpy.random.uniform(numpy.log10(min_c), numpy.log10(max_c))
        c_i = 10**c_i
        sigma_i = numpy.random.uniform(0.2, 1.99)

        try:
            n = RunWithTimeout(diversity_utils.Klogn, (mean, c_i))
            mu, s = n.run(10)

        except:
            continue

        S_tot = int(2*S_obs / erf((numpy.log(c_i)-mu) / numpy.sqrt(2)*s ))

        if (S_tot < 0):
            continue

        s_by_s_null, s_by_s_null_sampled = simulation_utils.generate_community(mu, s, S_tot, N, 'constant', sigma_i, len(N))

        n_non_zero = sum(~numpy.all(s_by_s_null_sampled == 0, axis=1))

        if fix_richness == True:

            #if n_non_zero < S_obs:
            #    continue


            if n_non_zero > S_obs:
                continue

            #print(S_obs, n_non_zero)

            # ignore if n_non_zero has 5% deviation
            #if (S_obs - n_non_zero)/S_obs > 0.1:
            #    continue

            # add "missing" species back in
            #s_by_s_null_sampled_no_zeros = s_by_s_null_sampled[~(numpy.all(s_by_s_null_sampled == 0, axis=1))]
            #n_zero = S_obs - s_by_s_null_sampled_no_zeros.shape[0]
            # put zeros back in to simulated number of species is equal to observed number
            #s_by_s_zeros = numpy.zeros((n_zero, s_by_s_null_sampled_no_zeros.shape[1]))
            #s_by_s_null_sampled = numpy.vstack((s_by_s_null_sampled_no_zeros, s_by_s_zeros))


            # get random number species
            #prob_species = s_by_s_null.sum(axis=1)
            #prob_species = prob_species/sum(prob_species)
            #idx_subsample = numpy.random.choice(s_by_s_null_sampled.shape[0], size=S_obs, replace=False, p=prob_species)

            # probability of sampling species determined by their relative abundance
            #s_by_s_null_subset = s_by_s_null[idx_subsample,:]
            #rel_s_by_s_null_subset = (s_by_s_null_subset/s_by_s_null_subset.sum(axis=0))

            #read_counts_all = []
            #for sad_idx, sad in enumerate(rel_s_by_s_null_subset.T):
            #    read_counts_all.append(numpy.random.multinomial(int(N[sad_idx]), sad))
            #s_by_s_null_sampled = numpy.asarray(read_counts_all).T


        diversity_species_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_null_sampled)

        euc_distance = numpy.linalg.norm(diversity_species-diversity_species_null)
        mean_euc_distance = euc_distance/len(diversity_species)

        # clear up the cache from sympy
        clear_cache()

        c_to_keep.append(c_i)
        sigma_to_keep.append(sigma_i)
        euc_distance_to_keep.append(mean_euc_distance)
        S_simulation_to_keep.append(n_non_zero)

        count += 1



        #if n_non_zero >= S_obs:

        #    if count_greater < iter:
        #        c_to_keep.append(c_i)
        #        sigma_to_keep.append(sigma_i)
        #        euc_distance_to_keep.append(mean_euc_distance)
        #        S_simulation_to_keep.append(n_non_zero)
        #        count_greater += 1

        #    else:
        #        continue


        #if n_non_zero <= S_obs:

        #    if count_less < iter:
        #        c_to_keep.append(c_i)
        #        sigma_to_keep.append(sigma_i)
        #        euc_distance_to_keep.append(mean_euc_distance)
        #        S_simulation_to_keep.append(n_non_zero)
        #        count_less += 1

        #    else:
        #        continue


        #if count % 100 == 0:
        #    print(count)

        #print(count_less + count_greater, count_less, count_greater, c_i, sigma_i, n_non_zero, S_obs)

        print(count, c_i, sigma_i, n_non_zero, S_obs)

        #count+=1


    parameter_dict = {}
    parameter_dict['c'] = c_to_keep
    parameter_dict['sigma'] = sigma_to_keep
    parameter_dict['mean_euc_distance'] = euc_distance_to_keep
    parameter_dict['S_simulation'] = S_simulation_to_keep


    sys.stderr.write("Saving parameter dictionary...\n")
    with open(path_, 'wb') as handle:
        pickle.dump(parameter_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




def load_slm_parameter_dict(environment, taxon=True, rarefied=True, fix_richness=True):

    if fix_richness == True:
        fix_richness_label = '_fix_richness'
    else:
        fix_richness_label = ''

    if taxon == True:
        sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
        path_ = slm_emp_taxon_parameter_dict_template % (diversity_utils.get_label(environment, rarefied), fix_richness_label)

    else:
        sad_annotated_dict = dbd_utils.load_sad_annotated_phylo_dict(environment, rarefied = rarefied)
        path_ = slm_emp_phylo_parameter_dict_template % (diversity_utils.get_label(environment, rarefied), fix_richness_label)


    with open(path_, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_








def predict_diversity_taxon(environment, rarefied, fix_richness = True, iter=100):

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
    samples = sad_annotated_dict['samples']

    sys.stderr.write("Getting site-by-species matrix...\n")
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
    mean = numpy.mean(rel_s_by_s, axis=1)
    N = s_by_s.sum(axis=0)
    S_obs = s_by_s.shape[0]

    sys.stderr.write("Selecting best parameters...\n")
    parameter_dict = load_slm_parameter_dict(environment, rarefied=rarefied)

    mean_euc_distance = parameter_dict['mean_euc_distance']
    sigma = parameter_dict['sigma']
    c = parameter_dict['c']
    S_simulation = parameter_dict['S_simulation']

    mean_euc_distance = numpy.asarray(mean_euc_distance)
    sigma = numpy.asarray(sigma)
    c = numpy.asarray(c)
    S_simulation = numpy.asarray(S_simulation)

    #if fix_richness == True:

    #    mean_euc_distance = mean_euc_distance[S_simulation >= S_obs]
    #    sigma = sigma[S_simulation >= S_obs]
    #    c = c[S_simulation >= S_obs]
    #    S_simulation = S_simulation[S_simulation >= S_obs]

    #else:

    #    mean_euc_distance = mean_euc_distance[S_simulation <= S_obs]
    #    sigma = sigma[S_simulation <= S_obs]
    #    c = c[S_simulation <= S_obs]
    #    S_simulation = S_simulation[S_simulation <= S_obs]


    #idx_min = mean_euc_distance.index(min(mean_euc_distance))
    idx_min = numpy.where(mean_euc_distance == numpy.amin(mean_euc_distance))[0][0]

    #idx_min = (numpy.abs(S_simulation - S_obs)).argmin()
    #print(mean_euc_distance[idx_min])

    #print(S_obs, S_simulation[idx_min])

    sigma_best = sigma[idx_min]
    c_best = c[idx_min]

    #c = 0.00005
    mu, s = diversity_utils.Klogn(mean, c = c_best)
    S_tot = int(2*S_obs / erf((numpy.log(c_best)-mu) / numpy.sqrt(2)*s))

    diversity_species = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s)

    diversity_dict = {}
    diversity_dict['asv'] = {}
    diversity_dict['asv']['observed'] = numpy.mean(diversity_species)
    diversity_dict['asv']['null'] = []

    for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

        all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
        all_genera_idx = numpy.arange(len(all_genera))

        genus_to_taxa_dict = {}
        sad_genera_all = []
        for genus in all_genera:
            genus_to_taxa_dict[genus] = []
            for t in taxa:
                if sad_annotated_dict['taxa'][t][rank] == genus:
                    genus_to_taxa_dict[genus].append(t)

            g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
            sad_genera_all.append(s_by_s[g_taxa_idx,:].sum(axis=0))


        s_by_s_genera = numpy.stack(sad_genera_all, axis=0)

        diversity_species_rank = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_genera)

        unique_genera, counts_genera = numpy.unique([sad_annotated_dict['taxa'][t][rank] for t in taxa], return_counts=True)
        coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_genera))[:-1]

        diversity_dict[rank] = {}
        diversity_dict[rank]['observed'] = numpy.mean(diversity_species_rank)
        diversity_dict[rank]['null'] = []
        diversity_dict[rank]['idx'] = coarse_grain_by_genus_idx


    count = 0
    while count < iter:

        if S_tot < 0:
            continue

        s_by_s_null, s_by_s_null_sampled = simulation_utils.generate_community(mu, s, S_tot, N, 'constant', sigma_best, len(N))
        n_non_zero = sum(~numpy.all(s_by_s_null_sampled == 0, axis=1))

        if fix_richness == True:

            #if n_non_zero < S_obs:
            #    continue

            #print(S_obs, n_non_zero)

            if n_non_zero > S_obs:
                continue

            # ignore if n_non_zero has 5% deviation
            #print((S_obs - n_non_zero)/S_obs)
            #if (S_obs - n_non_zero)/S_obs > 0.1:
            #    continue

            # add "missing" species back in
            s_by_s_null_sampled_no_zeros = s_by_s_null_sampled[~(numpy.all(s_by_s_null_sampled == 0, axis=1))]
            n_zero = S_obs - s_by_s_null_sampled_no_zeros.shape[0]
            # put zeros back in to simulated number of species is equal to observed number
            s_by_s_zeros = numpy.zeros((n_zero, s_by_s_null_sampled_no_zeros.shape[1]))
            s_by_s_null_sampled = numpy.vstack((s_by_s_null_sampled_no_zeros, s_by_s_zeros))


            #subsample_fixed_richness(s_by_s_null, N, S_obs)

            # get random number species
            #prob_species = s_by_s_null.sum(axis=1)
            #prob_species = prob_species/sum(prob_species)
            #idx_subsample = numpy.random.choice(s_by_s_null_sampled.shape[0], size=S_obs, replace=False)#, p=prob_species)

            # probability of sampling species determined by their relative abundance

            #s_by_s_null_subset = s_by_s_null[idx_subsample,:]
            #rel_s_by_s_null_subset = (s_by_s_null_subset/s_by_s_null_subset.sum(axis=0))

            #read_counts_all = []
            #for sad_idx, sad in enumerate(rel_s_by_s_null_subset.T):
            #    read_counts_all.append(numpy.random.multinomial(int(N[sad_idx]), sad))

            #s_by_s_null_sampled = numpy.asarray(read_counts_all).T

            #n_non_zero_subsampled = sum(~numpy.all(s_by_s_null_sampled == 0, axis=1))

            # add "missing" species back in
            #if n_non_zero_subsampled < S_obs:
            #    s_by_s_null_sampled_no_zeros = s_by_s_null_sampled[~(numpy.all(s_by_s_null_sampled == 0, axis=1))]
            #    n_zero = S_obs - s_by_s_null_sampled_no_zeros.shape[0]
            #    # put zeros back in to simulated number of species is equal to observed number
            #    s_by_s_zeros = numpy.zeros((n_zero, s_by_s_null_sampled_no_zeros.shape[1]))
            #    s_by_s_null_sampled = numpy.vstack((s_by_s_null_sampled_no_zeros, s_by_s_zeros))

            # susample one at a time
            #s_by_s_
            #N_copy = numpy.copy(N)
            #S_subsampled = []
            #while len(S_subsampled) < S_obs:

            #    for N_i in range(len(N_copy)):

            #s_by_s_null_sampled_subset = s_by_s_null_sampled[idx_subsample,:]

            # get reads to resample
            #N_to_resample = N - s_by_s_null_sampled_subset.sum(axis=0)

            # normalize subsetted SADs
            #rel_s_by_s_null_subset = (s_by_s_null_subset/s_by_s_null_subset.sum(axis=0))


        else:

            # we cant coarse grain more species than we observed
            if n_non_zero > S_obs:
                continue

            # add "missing" species back in
            if n_non_zero < S_obs:

                s_by_s_null_sampled_no_zeros = s_by_s_null_sampled[~(numpy.all(s_by_s_null_sampled == 0, axis=1))]
                n_zero = S_obs - s_by_s_null_sampled_no_zeros.shape[0]
                # put zeros back in to simulated number of species is equal to observed number
                s_by_s_zeros = numpy.zeros((n_zero, s_by_s_null_sampled_no_zeros.shape[1]))
                s_by_s_null_sampled = numpy.vstack((s_by_s_null_sampled_no_zeros, s_by_s_zeros))


        # shuffle for good measure
        numpy.random.shuffle(s_by_s_null_sampled)

        diversity_species_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_null_sampled)
        diversity_dict['asv']['null'].append(numpy.mean(diversity_species_null))

        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):
            s_by_s_null_sampled_coarse = numpy.add.reduceat(s_by_s_null_sampled, diversity_dict[rank]['idx'], axis=0)
            diversity_species_coarse_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_null_sampled_coarse)
            diversity_dict[rank]['null'].append(numpy.mean(diversity_species_coarse_null))

        if count % 100 == 0:
            print(count)

        count += 1


    for rank in diversity_dict.keys():

        null = numpy.asarray(diversity_dict[rank]['null'])
        null = numpy.sort(null)

        lower_ci = null[int(iter*0.025)]
        upper_ci = null[int(iter*0.975)]

        diversity_dict[rank]['lower_ci'] = lower_ci
        diversity_dict[rank]['upper_ci'] = upper_ci
        diversity_dict[rank]['mean_null'] = numpy.mean(null)

        del diversity_dict[rank]['null']

        if rank != 'asv':

            del diversity_dict[rank]['idx']

        print(rank, diversity_dict[rank]['observed'], lower_ci, upper_ci)


    diversity_slm_taxon_dict_path = slm_emp_taxon_diversity_dict_template % (diversity_utils.get_label(environment, rarefied))
    sys.stderr.write("Saving diversity dictionary...\n")
    with open(diversity_slm_taxon_dict_path, 'wb') as handle:
        pickle.dump(diversity_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_predict_diversity_taxon_dict(environment, rarefied=False):

    diversity_slm_taxon_dict_path = slm_emp_taxon_diversity_dict_template % (diversity_utils.get_label(environment, rarefied))

    with open(diversity_slm_taxon_dict_path, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_



def predict_diversity_phylo(environment, rarefied, fix_richness = True, iter=1000):

    parameter_dict = load_slm_parameter_dict(environment, rarefied=rarefied, taxon=False)

    mean_euc_distance = parameter_dict['mean_euc_distance']
    sigma = parameter_dict['sigma']
    c = parameter_dict['c']

    idx_min = mean_euc_distance.index(min(mean_euc_distance))
    sigma_best = sigma[idx_min]
    c_best = c[idx_min]

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_phylo_dict(environment, rarefied = rarefied)
    samples = sad_annotated_dict['samples']

    sys.stderr.write("Getting site-by-species matrix...\n")
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    #diversity_species = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s)

    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
    mean = numpy.mean(rel_s_by_s, axis=1)
    #c = 0.00005
    mu, s = diversity_utils.Klogn(mean, c = c_best)

    N = s_by_s.sum(axis=0)
    S_obs = s_by_s.shape[0]
    S_tot = int(2*S_obs / erf((numpy.log(c_best)-mu) / numpy.sqrt(2)*s))

    sad_annotated_phylo_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)
    distances = list(sad_annotated_phylo_dict.keys())
    distances.sort()

    # get observed mean diversity
    diversity_dict = {}
    # save coarse grained indices for null simulation
    coarse_grained_dict = {}
    for distance in distances:
        sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))

        coarse_grained_list = sad_annotated_phylo_dict[distance]
        coarse_grained_n = [len(i) for i in coarse_grained_list]

        # get indexes for each clade
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
        # coarse grain s-by-s for all clades
        s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
        diversity_clades = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_all_clades)

        diversity_dict[distance] = {}
        diversity_dict[distance]['observed'] = numpy.mean(diversity_clades)
        diversity_dict[distance]['null'] = []

        coarse_grain_list = sad_annotated_phylo_dict[distance]
        coarse_grain_n = [len(i) for i in coarse_grain_list]
        coarse_grain_idx = numpy.append([0], numpy.cumsum(coarse_grain_n))[:-1]
        coarse_grained_dict[distance] = coarse_grain_idx


    count = 0
    while count < iter:

        s_by_s_null, s_by_s_null_sampled = simulation_utils.generate_community(mu, s, S_tot, N, 'constant', sigma_best, len(N))

        n_non_zero = sum(~numpy.all(s_by_s_null_sampled == 0, axis=1))

        if fix_richness == True:

            #if n_non_zero < S_obs:
            #    continue

            if n_non_zero > S_obs:
                continue

            # add "missing" species back in
            if n_non_zero < S_obs:
                s_by_s_null_sampled_no_zeros = s_by_s_null_sampled[~(numpy.all(s_by_s_null_sampled == 0, axis=1))]
                n_zero = S_obs - s_by_s_null_sampled_no_zeros.shape[0]
                # put zeros back in to simulated number of species is equal to observed number
                s_by_s_zeros = numpy.zeros((n_zero, s_by_s_null_sampled_no_zeros.shape[1]))
                s_by_s_null_sampled = numpy.vstack((s_by_s_null_sampled_no_zeros, s_by_s_zeros))


        # randomize rows
        numpy.random.shuffle(s_by_s_null_sampled)

        for distance in distances:
            s_by_s_null_sampled_coarse = numpy.add.reduceat(s_by_s_null_sampled, coarse_grained_dict[distance], axis=0)
            diversity_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_null_sampled_coarse)
            diversity_dict[distance]['null'].append(numpy.mean(diversity_null))


        if count % 100 == 0:
            print(count)

        count += 1


    #for distance in diversity_dict.keys():
    for distance in distances:

        null = numpy.asarray(diversity_dict[distance]['null'])
        null = numpy.sort(null)

        lower_ci = null[int(iter*0.025)]
        upper_ci = null[int(iter*0.975)]

        diversity_dict[distance]['lower_ci'] = lower_ci
        diversity_dict[distance]['upper_ci'] = upper_ci
        diversity_dict[distance]['mean_null'] = numpy.mean(null)

        del diversity_dict[distance]['null']

        print(diversity_dict[distance]['observed'], lower_ci, upper_ci)


    diversity_slm_phylo_dict_path = slm_emp_phylo_diversity_dict_template % (diversity_utils.get_label(environment, rarefied))
    sys.stderr.write("Saving diversity dictionary...\n")
    with open(diversity_slm_phylo_dict_path, 'wb') as handle:
        pickle.dump(diversity_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_predict_diversity_phylo_dict(environment, rarefied=False):

    diversity_slm_taxon_dict_path = slm_emp_phylo_diversity_dict_template % (diversity_utils.get_label(environment, rarefied))

    with open(diversity_slm_taxon_dict_path, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_







def predict_diversity_gamma(environment, rarefied=False, iter=100):

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    #sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)
    samples = sad_annotated_dict['samples']

    sys.stderr.write("Getting site-by-species matrix...\n")
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
    mean = numpy.mean(rel_s_by_s, axis=1)
    var = numpy.var(rel_s_by_s, axis=1)
    N_reads = s_by_s.sum(axis=0)

    # get means and variances from MLE estimator.

    mle_dict = mle_utils.load_run_mle_gamma_taxon_dict()

    mean_mle = [mle_dict[environment][t]['mu_mle'] for t in taxa]
    var_mle = [mle_dict[environment][t]['var_mle'] for t in taxa]
    mean_mle = numpy.asarray(mean_mle)
    var_mle = numpy.asarray(var_mle)


    corr_bw_asvs = numpy.corrcoef(rel_s_by_s)
    corr_bw_sites = numpy.corrcoef(rel_s_by_s.T)
    I_bw_asvs = numpy.identity(len(mean))
    I_bw_sites = numpy.identity(len(N_reads))

    diversity = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s)
    richness = numpy.apply_along_axis(diversity_utils.calculate_richness, 0, s_by_s)
    diversity_sim_all = []
    #diversity_sim_corr_all = []
    for i in range(iter):

        #s_by_s_sim, s_by_s_reads_sim = simulation_utils.genrate_community_from_mean_and_var(mean, var, N_reads, len(N_reads), rhogamma=0, corr_matrix=I)

        #s_by_s_sim, s_by_s_reads_sim = simulation_utils.genrate_community_from_mean_and_var_using_rvs(mean, var*1.3, N_reads, len(N_reads))
        s_by_s_sim, s_by_s_reads_sim = simulation_utils.genrate_community_from_mean_and_var_using_rvs(mean_mle, var_mle, N_reads, len(N_reads))

        #s_by_s_sim_corr_bw_sites, s_by_s_reads_sim_corr_bw_sites = simulation_utils.genrate_community_from_mean_and_var_random_sites(mean, var, N_reads, len(N_reads), corr_matrix=corr_bw_sites)
        #s_by_s_sim, s_by_s_reads_sim = simulation_utils.genrate_community_from_mean_and_var_random_sites(mean, var, N_reads, len(N_reads), corr_matrix=I_bw_sites)

        diversity_sim = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_reads_sim)
        richness_sim = numpy.apply_along_axis(diversity_utils.calculate_richness, 0, s_by_s_reads_sim)
        #diversity_sim_corr = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_reads_sim_corr_bw_sites)
        #diversity_sim_corr = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_reads_sim_corr)
        #diversity_dict[rank]['null'].append(numpy.mean(diversity_coarse))

        diversity_sim_all.append(numpy.mean(diversity_sim))
        #diversity_sim_corr_all.append(diversity_sim_corr)

        print(numpy.mean(richness), numpy.mean(richness_sim))


    diversity_sim_all = numpy.asarray(diversity_sim_all)
    mean_diverity_slm = numpy.mean(diversity_sim_all, axis=0)

    print(numpy.mean(diversity), numpy.mean(diversity_sim_all))

    #diversity_sim_corr_all = numpy.asarray(diversity_sim_corr_all)
    #mean_diverity_corr_slm = numpy.mean(diversity_sim_corr_all, axis=0)


    #print(numpy.mean(diversity), mean_diverity_slm)



    #print(numpy.absolute(mean_diverity_corr_slm - diversity)/diversity)


    diversity_dict = {}
    diversity_dict['observed'] = diversity.tolist()
    diversity_dict['mean_null'] = mean_diverity_slm.tolist()
    #diversity_dict['mean_corr_null'] = mean_diverity_corr_slm.tolist()


    diversity_slm_dict_path = slm_gamma_emp_diversity_dict_template % (diversity_utils.get_label(environment, rarefied))
    sys.stderr.write("Saving diversity dictionary...\n")
    with open(diversity_slm_dict_path, 'wb') as handle:
        pickle.dump(diversity_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_diversity_gamma(environment, rarefied=False):

    diversity_slm_dict_path = slm_gamma_emp_diversity_dict_template % (diversity_utils.get_label(environment, rarefied))

    with open(diversity_slm_dict_path, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_




def predict_diversity_gamma_taxon(environment, rarefied, fix_richness = True, iter=100):

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    #sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = rarefied)
    samples = sad_annotated_dict['samples']

    sys.stderr.write("Getting site-by-species matrix...\n")
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
    mean = numpy.mean(rel_s_by_s, axis=1)
    var = numpy.var(rel_s_by_s, axis=1)
    N_reads = s_by_s.sum(axis=0)

    corr_bw_asvs = numpy.corrcoef(rel_s_by_s)
    corr_bw_sites = numpy.corrcoef(rel_s_by_s.T)
    I_bw_asvs = numpy.identity(len(mean))
    I_bw_sites = numpy.identity(len(N_reads))

    diversity = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s)


    #cov = cov[numpy.triu_indices(cov.shape[0], k = 1)]

    #print(max(cov))

    #mle_dict = mle_utils.load_run_mle_gamma_taxon_dict()

    #mean_mle = [mle_dict[environment][t]['mu_mle'] for t in taxa]
    #var_mle = [mle_dict[environment][t]['var_mle'] for t in taxa]
    #mean_mle = numpy.asarray(mean_mle)
    #var_mle = numpy.asarray(var_mle)

    diversity_dict = {}
    for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

        all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
        all_genera_idx = numpy.arange(len(all_genera))

        genus_to_taxa_dict = {}
        sad_genera_all = []
        g_taxa_idx_all = []
        g_taxa_idx_flat = []
        for genus in all_genera:
            genus_to_taxa_dict[genus] = []
            for t in taxa:
                if sad_annotated_dict['taxa'][t][rank] == genus:
                    genus_to_taxa_dict[genus].append(t)

            g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
            g_taxa_idx_all.append(g_taxa_idx)
            g_taxa_idx_flat.extend(g_taxa_idx.tolist())
            sad_genera_all.append(s_by_s[g_taxa_idx,:].sum(axis=0))



        s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
        diversity_species_rank = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_genera)

        #g_taxa_idx_flat = numpy.asarray(g_taxa_idx_flat)
        #n_coarse_grained_all = numpy.asarray([len(c) for c in g_taxa_idx_all])
        #coarse_grain_idx = numpy.append([0], numpy.cumsum(n_coarse_grained_all))[:-1]

        # sort mean
        #mean_sort = mean[g_taxa_idx_flat]
        #var_sort = var[g_taxa_idx_flat]
        #N_reads_sort = N_reads[g_taxa_idx_flat]


        diversity_dict[rank] = {}
        diversity_dict[rank]['observed'] = numpy.mean(diversity_species_rank)
        diversity_dict[rank]['null'] = []
        diversity_dict[rank]['null_corr'] = []
        #diversity_dict[rank]['null_mle'] = []

        #diversity_dict[rank]['coarse_grain_cumulative_n'] = coarse_grain_idx
        diversity_dict[rank]['idx'] = g_taxa_idx_all
        #diversity_dict[rank]['mean'] = mean_sort
        #diversity_dict[rank]['var'] = var_sort


    for i in range(iter):

        #s_by_s_sim, s_by_s_reads_sim = simulation_utils.genrate_community_from_mean_and_var(mean, var, N_reads, len(N_reads), rhogamma=0, corr_matrix=I)

        s_by_s_sim, s_by_s_reads_sim = simulation_utils.genrate_community_from_mean_and_var(mean, var, N_reads, len(N_reads), rhogamma=0, corr_matrix=corr_bw_asvs)

        #s_by_s_sim, s_by_s_reads_sim = simulation_utils.genrate_community_from_mean_and_var_using_rvs(mean, var, N_reads, len(N_reads))
        #s_by_s_sim_corr_bw_sites, s_by_s_reads_sim_corr_bw_sites = simulation_utils.genrate_community_from_mean_and_var_random_sites(mean, var, N_reads, len(N_reads), corr_matrix=corr_bw_sites)
        #s_by_s_sim, s_by_s_reads_sim = simulation_utils.genrate_community_from_mean_and_var_random_sites(mean, var, N_reads, len(N_reads), corr_matrix=I_bw_sites)

        diversity_sim = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_reads_sim)
        #diversity_sim_corr = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_reads_sim_corr)
        #diversity_dict[rank]['null'].append(numpy.mean(diversity_coarse))

        print(numpy.mean(diversity), numpy.mean(diversity_sim))

        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            s_by_s_reads_sim_coarse = numpy.stack([numpy.sum(s_by_s_reads_sim[coarse_grained_idx,:], axis=0) for coarse_grained_idx in diversity_dict[rank]['idx']], axis=0)
            diversity_sim_coarse = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_reads_sim_coarse)
            diversity_dict[rank]['null'].append(numpy.mean(diversity_sim_coarse))


            #s_by_s_reads_sim_corr_coarse = numpy.stack([numpy.sum(s_by_s_reads_sim_corr[coarse_grained_idx,:], axis=0) for coarse_grained_idx in diversity_dict[rank]['idx']], axis=0)
            #diversity_corr_coarse = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_reads_sim_corr_coarse)
            #diversity_dict[rank]['null_corr'].append(numpy.mean(diversity_corr_coarse))

            #print(diversity_dict[rank]['observed'], numpy.mean(diversity_sim_coarse), numpy.mean(diversity_corr_coarse))

            error = numpy.absolute(diversity_dict[rank]['observed'] - numpy.mean(diversity_sim_coarse))/diversity_dict[rank]['observed']
            #error_corr = numpy.absolute(diversity_dict[rank]['observed'] - numpy.mean(diversity_corr_coarse))/diversity_dict[rank]['observed']

            print(rank, diversity_dict[rank]['observed'], numpy.mean(diversity_sim_coarse), error)

            #print(rank, error, error_corr)

            #print(diversity_dict[rank]['observed'], numpy.mean(diversity_sim_coarse))

            #s_by_s_reads_sim_mle_coarse = numpy.stack([numpy.sum(s_by_s_reads_mle_sim[coarse_grained_idx,:], axis=0) for coarse_grained_idx in diversity_dict[rank]['idx']], axis=0)
            #diversity_coarse_mle = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_reads_sim_mle_coarse)
            #diversity_dict[rank]['null_mle'].append(numpy.mean(diversity_coarse_mle))






        #for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

        #mean_rank = diversity_dict[rank]['mean']
        #var_rank = diversity_dict[rank]['var']
        #coarse_grain_rank_idx = diversity_dict[rank]['idx']

        #for i in range(iter):
        #    s_by_s_sim, s_by_s_reads_sim = simulation_utils.genrate_community_from_mean_and_var(mean, var, N_reads, len(N_reads))

        #    rel_s_by_s_constrain_mad_coarse = numpy.add.reduceat(s_by_s_reads_sim, coarse_grain_rank_idx, axis=0)

        #    diversity_coarse = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_s_by_s_constrain_mad_coarse)
        #    diversity_dict[rank]['null'].append(numpy.mean(diversity_coarse))


    for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

        null = numpy.asarray(diversity_dict[rank]['null'])
        null = numpy.sort(null)

        #null_mle = numpy.asarray(diversity_dict[rank]['null_mle'])
        #null_mle = numpy.sort(null_mle)

        lower_ci = null[int(iter*0.025)]
        upper_ci = null[int(iter*0.975)]

        #lower_ci_mle = null_mle[int(iter*0.025)]
        #upper_ci_mle = null_mle[int(iter*0.975)]

        diversity_dict[rank]['lower_ci'] = lower_ci
        diversity_dict[rank]['upper_ci'] = upper_ci
        diversity_dict[rank]['mean_null'] = numpy.mean(null)



        null_corr = numpy.asarray(diversity_dict[rank]['null_corr'])
        null_corr = numpy.sort(null_corr)

        corr_lower_ci = null_corr[int(iter*0.025)]
        corr_upper_ci = null_corr[int(iter*0.975)]

        diversity_dict[rank]['lower_ci_corr'] = corr_lower_ci
        diversity_dict[rank]['upper_ci_corr'] = corr_upper_ci
        diversity_dict[rank]['mean_null_corr'] = numpy.mean(null_corr)



        #diversity_dict[rank]['lower_ci_mle'] = lower_ci_mle
        #diversity_dict[rank]['upper_ci_mle'] = upper_ci_mle
        #diversity_dict[rank]['mean_null_mle'] = numpy.mean(null_mle)

        #del diversity_dict[rank]['null']
        #del diversity_dict[rank]['mean_null_mle']

        if rank != 'asv':

            del diversity_dict[rank]['idx']
            #del diversity_dict[rank]['mean']
            #del diversity_dict[rank]['var']
            #del diversity_dict[rank]['coarse_grain_cumulative_n']


        print(rank, diversity_dict[rank]['observed'], lower_ci_mle, upper_ci_mle)


    diversity_slm_taxon_dict_path = slm_gamma_emp_taxon_diversity_dict_template % (diversity_utils.get_label(environment, rarefied))
    sys.stderr.write("Saving diversity dictionary...\n")
    with open(diversity_slm_taxon_dict_path, 'wb') as handle:
        pickle.dump(diversity_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_predict_diversity_gamma_taxon(environment, rarefied=False):

    diversity_slm_taxon_dict_path = slm_gamma_emp_taxon_diversity_dict_template % (diversity_utils.get_label(environment, rarefied))

    with open(diversity_slm_taxon_dict_path, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_
