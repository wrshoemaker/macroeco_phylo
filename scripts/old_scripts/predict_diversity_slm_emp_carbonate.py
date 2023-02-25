import os
#from Bio import Phylo
import random
import copy
import sys
import numpy
import random
import pickle
import scipy.stats as stats

import dbd_utils

import diversity_utils
import config

from scipy.special import erf
from scipy.special import rel_entr

import simulation_utils

import datetime
import signal
import threading
import time





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



numpy.random.seed(123456789)



slm_emp_parameter_dict_template = config.data_directory + "slm_emp_parameter_dict/slm_emp_parameter_dict_%s.pickle"


environment = str(sys.argv[1])
environment = environment.replace('_', ' ')

rarefied = True
iter=1000
# max value of 2


sys.stderr.write("Subsetting samples for %s...\n" % environment)
sad_annotated_dict = dbd_utils.load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
samples = sad_annotated_dict['samples']

sys.stderr.write("Getting site-by-species matrix...\n")

taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

diversity_species = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s)

#mu=-16
#s=3
rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
mean = numpy.mean(rel_s_by_s, axis=1)

#min_c = 1/max(s_by_s.sum(axis=0))
min_c = 0.000001
max_c = (1/min(s_by_s.sum(axis=0)))*0.1

#c_range = numpy.logspace(numpy.log10(min_c), numpy.log10(max_c), num=50, endpoint=True, base=10)
#sigma_range = numpy.linspace(0.1, 2, num=50, endpoint=False)

#c_all = numpy.random.uniform(numpy.log10(min_c), numpy.log10(max_c), size=iter)
# back transform
#c_all = 10**c_all

#sigma_all = numpy.random.uniform(0.1, 1.99, size=iter)

c_to_keep = []
sigma_to_keep = []
euc_distance_to_keep = []

count = 0
#for i in range(iter):
while count < iter:

    #start = time.time()
    c_i = numpy.random.uniform(numpy.log10(min_c), numpy.log10(max_c))
    c_i = 10**c_i
    sigma_i = numpy.random.uniform(0.1, 1.99)

    #c_i = c_all[i]
    #sigma_i = sigma_all[i]

    #sys.stdout.write("%s\t%s\n" % (c_i, sigma_i))
    #sys.stdout.flush()

    #signal.alarm(10)
    #start = datetime.datetime.now()

    #n = RunWithTimeout(foo, (5,3))

    #mu, s = diversity_utils.Klogn(mean, c = c_i)

    try:
        n = RunWithTimeout(diversity_utils.Klogn, (mean, c_i))
        mu, s = n.run(10)

    except:
        continue


    #try:
    #    # Whatever your function that might hang
    #    # use S
    #    mu, s = diversity_utils.Klogn(mean, c = c_i)

    #    end = datetime.datetime.now()
    #    seconds = end - start
    #    #print str(c.seconds) + " seconds"
    #except:
    #    continue # continue the for loop if function takes more than x seconds
    #else:
    #    # Reset the alarm
    #    signal.alarm(0)


    #try:
    #    mu, s = diversity_utils.Klogn(mean, c = c_i)
    #except:
    #    continue

    N = s_by_s.sum(axis=0)
    #S = s_by_s.shape[0]
    S_obs = s_by_s.shape[0]
    S_tot = int(2*S_obs / erf((numpy.log(c_i)-mu) / numpy.sqrt(2)*s ))

    if S_tot < 0:
        continue

    s_by_s_null = simulation_utils.generate_community(mu, s, S_tot, N, 'constant', sigma_i, len(N))
    diversity_species_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_null)
    #diversity_species_null_all.append(diversity_species_null)

    #idx_to_keep = (diversity_species>0) & (diversity_species_null>0)
    #rel_entr_sum = sum(numpy.absolute(rel_entr(diversity_species[idx_to_keep], diversity_species_null[idx_to_keep])))
    euc_distance = numpy.linalg.norm(diversity_species-diversity_species_null)
    mean_euc_distance = euc_distance/len(diversity_species)

    c_to_keep.append(c_i)
    sigma_to_keep.append(sigma_i)
    euc_distance_to_keep.append(mean_euc_distance)

    #print(seconds)

    #end = time.time()
    #sys.stderr.write("%s\t%s\n" % (c_i, sigma_i))
    #sys.stderr.flush()

    if count % 100 == 0:
        print(count)

    count+=1



parameter_dict = {}
parameter_dict['c'] = c_to_keep
parameter_dict['sigma'] = sigma_to_keep
parameter_dict['mean_euc_distance'] = euc_distance_to_keep


path_ = slm_emp_parameter_dict_template % (diversity_utils.get_label(environment, rarefied))
sys.stderr.write("Saving slope dictionary...\n")
with open(path_, 'wb') as handle:
    pickle.dump(parameter_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
