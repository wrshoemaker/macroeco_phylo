from __future__ import division
import config
import os
import sys
import subprocess
import re
import itertools
import pickle
import numpy
import diversity_utils
import scipy.stats as stats
import dbd_utils



rarefied = False




for environment in diversity_utils.environments_to_keep:

    sys.stderr.write("Running analyses for %s...\n" % environment)

    # subsample EMP data
    sys.stderr.write("Subsampling EMP data...\n")
    dbd_utils.make_sad_annotated_taxon_dict(environment=environment, rarefied=False, n_subsample=100, n_reads=5000)
    

    # get full dataset using set of samples chosen in make_sad_annotated_taxon_dict
    sys.stderr.write("Genrating species abundance distributions...\n")
    dbd_utils.make_sad_annotated_no_subsampling_taxon_dict(environment)
    dbd_utils.make_sad_annotated_no_subsampling_phylo_dict(environment)


    # generate subsampled dict using all ASVs using same sites used in taxon SADs
    sys.stderr.write("Genrating species abundance distributions...\n")
    dbd_utils.make_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=False)


    # estimate the expected contribution towards diversity of each taxon in each site via numerical integration
    sys.stderr.write("Estimating integrals for Shannon's diversity predictions...\n")
    dbd_utils.make_diversity_slm_dbd_integral_dict(environment)
    

    # predict statistical moments of richness and diversity for phylogenetic coarse-graining
    sys.stderr.write("Predicting within-phylogenetic scale richness and diversity...\n")
    dbd_utils.make_richness_diversity_prediction_phylo_dict(environment, iter_=1000)


    # predict DBD slope of diversity estimates using taxonomic and phylogenetic coarse-graining
    sys.stderr.write("Predicting DBD slopes for diversity estimates...\n")
    dbd_utils.make_diversity_slm_dbd_dict(environment)


    # simulate null site-by-OTU read count matrices using correlatd gamma distributions
    # for phylogenetic coarse-graining AND generate DBD prediction
    sys.stderr.write("Predicting phylogenetic DBD slope for diversity estimates via simulations...\n")
    dbd_utils.run_svd_and_make_diversity_slm_dbd_simulation_phylo_dict(environment, iter_=1000)


    


####################################
# functions that run all environments
#####################################


# generate singular value decomposition dictionary for taxonomic coarse-graining
sys.stderr.write("Building SVD dictionaries for taxonomic data...\n")
dbd_utils.make_svd_taxon_dict()


# predict statistical moments of richness and diversity for taxonomic coarse-graining
sys.stderr.write("Predicting within-taxonomic scale richness and diversity...\n")
dbd_utils.make_richness_diversity_prediction_taxon_dict()


# predict DBD slope of richness estimates using taxonomic and phylogenetic coarse-graining
sys.stderr.write("Predicting DBD slope for richness estimates...\n")
dbd_utils.make_richness_slm_dbd_dict()


# simulate null site-by-OTU read count matrices using correlatd gamma distributions
# for taxonomic coarse-graining AND generate DBD prediction 
sys.stderr.write("Predicting taxonomic DBD slope for diversity estimates via simulations...\n")
dbd_utils.make_diversity_slm_dbd_simulation_taxon_dict()


