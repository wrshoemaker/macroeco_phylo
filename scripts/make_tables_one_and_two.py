from __future__ import division
import config
import os
import sys
import subprocess
import random
import re
import itertools
import pickle
from collections import Counter

import numpy
import diversity_utils
import simulation_utils

import scipy.stats as stats
import scipy.integrate as integrate

#import ete3
import tree_utils
import random
import dbd_utils



# table for taxonomic coarse-graining

output_file_taxon = open(config.data_directory + "metadata_taxon.csv", "w")
output_file_taxon.write(", ".join(["Environment", "Total #OTUs", "Mean #OTUs", "Mean #reads"]))


for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    sys.stderr.write("Loading taxon dict for %s...\n" % environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied=False)
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    n_otus_taxon = s_by_s.shape[0]
    mean_n_otus = numpy.mean(numpy.sum(s_by_s>0, axis=0))
    mean_n_reads = numpy.mean(numpy.sum(s_by_s, axis=0))

    print(numpy.sum(s_by_s, axis=0))

    output_file_taxon.write("\n")
    output_file_taxon.write("%s, %d, %0.2f, %0.2f" % (environment, n_otus_taxon, mean_n_otus, mean_n_reads))




output_file_taxon.close()


output_file_phylo = open(config.data_directory + "metadata_phylo.csv", "w")
output_file_phylo.write(", ".join(["Environment", "Total #OTUs", "Mean #OTUs", "Mean #reads"]))

# table for phylogenetic coarse-graining
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    sys.stderr.write("Loading taxon dict for %s...\n" % environment)
    sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied=False)
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    n_otus_taxon = s_by_s.shape[0]
    mean_n_otus = numpy.mean(numpy.sum(s_by_s>0, axis=0))
    mean_n_reads = numpy.mean(numpy.sum(s_by_s, axis=0))

    output_file_phylo.write("\n")
    output_file_phylo.write("%s, %d, %0.2f, %0.2f" % (environment, n_otus_taxon, mean_n_otus, mean_n_reads))



output_file_phylo.close()