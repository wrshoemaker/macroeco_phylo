from __future__ import division
import config
import os
import sys
import re
import itertools
import pickle
import numpy
import diversity_utils
import scipy.stats as stats
#import  matplotlib.pyplot as plt
import ete3
import tree_utils


#sad_path_template = config.data_directory + "sad_annotated_%s.pickle"
slope_path_template = config.data_directory + "dbd_slope_emp_phylo/dbd_slope%s.pickle"
coarse_grained_template = config.data_directory + "coarse_grained_%s.pickle"


environment='human gut metagenome'
iter = 1000
distance = float(sys.argv[1])


def load_sad_annotated_phylo_dict():

    coarse_grained_path = coarse_grained_template % environment.replace(' ', '_')

    with open(coarse_grained_path, 'rb') as handle:
        distance_collapsed_dict = pickle.load(handle)
    return distance_collapsed_dict




sad_annotated_phylo_dict = load_sad_annotated_phylo_dict()

distances = list(sad_annotated_phylo_dict.keys())
distances.sort()

distances = [distance]
#dbd_slope_dict = {}
sys.stderr.write("Subsetting samples...\n")
samples = diversity_utils.subset_observations(environment=environment)

sys.stderr.write("Getting site-by-species matrix...\n")
s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)
for distance in distances:

    dbd_slope_dict = {}

    dbd_slope_dict = {}
    dbd_slope_dict['slope'] = []
    dbd_slope_dict['collapsed_nodes'] = []
    dbd_slope_dict['null_dist'] = []
    dbd_slope_dict['p_value'] = []
    dbd_slope_dict['percentile'] = []
    dbd_slope_dict['lower_ci'] = []
    dbd_slope_dict['upper_ci'] = []
    dbd_slope_dict['n_asvs_focal_genus'] = []
    dbd_slope_dict['standard_score'] = []
    dbd_slope_dict['asv'] = []

    sys.stderr.write("Phylo distance = %s \n" % str(round(distance, 7)))

    coarse_grained_list = sad_annotated_phylo_dict[distance]
    coarse_grained_n = [len(i) for i in coarse_grained_list]

    # get indexes for each clade
    coarse_grained_idx_all = [numpy.asarray([numpy.where(taxonomy_names==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
    # coarse grain s-by-s for all clades
    s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)

    idx_all = numpy.arange(len(coarse_grained_idx_all))
    for i in idx_all:

        if len(coarse_grained_idx_all[i]) < 10:
            continue

        s_by_s_focal = s_by_s[coarse_grained_idx_all[i],:]
        s_by_s_non_focal = s_by_s_all_clades[idx_all!=i,:]

        diversity_focal = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_focal)
        diversity_non_focal = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_non_focal)

        idx_to_keep = (diversity_focal>0) & (diversity_non_focal>0) & (numpy.isfinite(diversity_focal)) & (numpy.isfinite(diversity_non_focal))
        diversity_focal = diversity_focal[idx_to_keep]
        diversity_non_focal = diversity_non_focal[idx_to_keep]

        # ignore cut-offs with low richness
        if len(diversity_focal) < 10:
            continue

        # ignore cut-offs where all diversity estimates are the same
        if (sum(diversity_non_focal==1) == len(diversity_non_focal)) or (sum(diversity_focal==1) == len(diversity_focal)):
            continue

        slope, intercept, r_value, p_value, std_err = stats.linregress(diversity_non_focal, diversity_focal)
        s_by_s_copy = numpy.copy(s_by_s)

        #print(i, len(idx_all))

        slope_null_all = []
        while len(slope_null_all) < iter:
            numpy.random.shuffle(s_by_s_copy)

            s_by_s_focal_null = s_by_s_copy[:coarse_grained_n[i],:]
            s_by_s_non_focal_null = s_by_s_copy[coarse_grained_n[i]:,:]

            diversity_focal_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_focal_null)
            diversity_non_focal_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_non_focal_null)

            idx_to_keep_null = (diversity_focal_null>0) & (diversity_non_focal_null>0) & (numpy.isfinite(diversity_focal_null)) & (numpy.isfinite(diversity_non_focal_null))
            diversity_focal_null = diversity_focal_null[idx_to_keep_null]
            diversity_non_focal_null = diversity_non_focal_null[idx_to_keep_null]

            # not enough
            if len(diversity_non_focal_null) < len(diversity_non_focal):
                continue

            # randomly sample indexes
            idx_subsample = numpy.random.choice(numpy.arange(len(diversity_non_focal_null)), size=len(diversity_non_focal_null), replace=True)

            diversity_focal_null_subsample = diversity_focal_null[idx_subsample]
            diversity_non_focal_null_subsample = diversity_non_focal_null[idx_subsample]

            slope_null, intercept_null, r_value_null, p_value_null, std_err_null = stats.linregress(diversity_non_focal_null_subsample, diversity_focal_null_subsample)
            slope_null_all.append(slope_null)

            if (len(slope_null_all) % 100 == 0):
                sys.stderr.write("Iteration %d...\n" % len(slope_null_all))


        slope_null_all = numpy.asarray(slope_null_all)
        slope_null_all = numpy.sort(slope_null_all)

        slope_null_median = numpy.median(slope_null_all)
        if slope >= slope_null_median:
            p_value = sum(slope_null_all > slope)/iter
        else:
            p_value = sum(slope_null_all < slope)/iter

        percentile = sum(slope_null_all > slope)/iter

        lower_ci = slope_null_all[int(iter*0.025)]
        upper_ci = slope_null_all[int(iter*0.975)]

        n_asvs_focal_genus = s_by_s_focal.shape[0]

        dbd_slope_dict['slope'].append(slope)
        dbd_slope_dict['collapsed_nodes'].append(coarse_grained_list[i])
        dbd_slope_dict['null_dist'].append(slope_null_all.tolist())
        dbd_slope_dict['p_value'].append(p_value)
        dbd_slope_dict['percentile'].append(percentile)
        dbd_slope_dict['lower_ci'].append(lower_ci)
        dbd_slope_dict['upper_ci'].append(upper_ci)
        dbd_slope_dict['n_asvs_focal_genus'].append(n_asvs_focal_genus)
        dbd_slope_dict['standard_score'].append((slope - numpy.mean(slope_null_all)) / numpy.std(slope_null_all))
        dbd_slope_dict['asv'].append(coarse_grained_list[i])

    #if distance in dbd_slope_dict:
    #    print('slope = ', numpy.mean(dbd_slope_dict[distance]['slope']))

slope_path = slope_path_template % ('_' + environment.replace(' ', '_') + '_phylo_' + str(distance))
sys.stderr.write("Saving slope dictionary...\n")
with open(slope_path, 'wb') as handle:
    pickle.dump(dbd_slope_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
