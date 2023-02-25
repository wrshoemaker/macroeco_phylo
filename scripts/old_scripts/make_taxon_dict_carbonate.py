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


sad_path_template = config.data_directory + "sad_annotated/sad_annotated_%s.pickle"
slope_emp_taxa_template = config.data_directory + "dbd_slope_emp_taxa/dbd_slope_%s.pickle"
#coarse_grained_template = config.data_directory + "coarse_grained_%s.pickle"


iter = 1000
environment = str(sys.argv[1])
environment = environment.replace('_', ' ')
rarefied = True


def load_sad_annotated_taxon_dict(environment, rarefied=False):

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    sad_path = sad_path_template % (environment.replace(' ', '_') + rarefied_label)

    with open(sad_path, 'rb') as handle:
        presence_absence_dict = pickle.load(handle)
    return presence_absence_dict




sys.stderr.write("Loading presence-absence dict...\n")
pres_abs_dict = load_sad_annotated_taxon_dict(environment, rarefied=rarefied)

dbd_slope_dict = {}
for rank in diversity_utils.taxa_ranks:
    dbd_slope_dict[rank] = {}

    sys.stderr.write("Starting %s level analysis...\n" % rank)

    sys.stderr.write("Estimating dbd slope...\n")
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
    all_genera = numpy.asarray(list(set([pres_abs_dict['taxa'][t][rank] for t in taxa])))
    all_genera_idx = numpy.arange(len(all_genera))

    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

    genus_to_taxa_dict = {}
    sad_genera_all = []

    for genus in all_genera:
        genus_to_taxa_dict[genus] = []
        for t in taxa:
            if pres_abs_dict['taxa'][t][rank] == genus:
                genus_to_taxa_dict[genus].append(t)

        g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
        sad_genera_all.append(s_by_s[g_taxa_idx,:].sum(axis=0))

    s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
    # remove sites where there are no observations
    s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]

    for focal_genus in all_genera:

        #presence_absence_all = []
        sad_all = []
        s_by_s_focal_genus = []
        sad_non_focal_genus = []
        non_focal_genera_all = []

        for t in taxa:
            genus_t = pres_abs_dict['taxa'][t][rank]

            if genus_t == focal_genus:
                s_by_s_focal_genus.append(pres_abs_dict['taxa'][t]['abundance'])

            else:
                non_focal_genera_all.append(genus_t)

        # dont look at genera with less than 10 ASVs
        n_asvs_focal_genus = len(s_by_s_focal_genus)
        if n_asvs_focal_genus < 10:
            continue

        s_by_s_focal_genus = numpy.asarray(s_by_s_focal_genus)
        # get s-by-s for non-focal
        focal_genus_idx = numpy.where(all_genera==focal_genus)[0][0]
        s_by_s_non_focal_genera = s_by_s_genera[all_genera_idx != focal_genus_idx,: ]

        # remove sites where there are no observations in either focal or non-focal
        idx_to_remove = numpy.all(s_by_s_focal_genus == 0, axis=0) | numpy.all(s_by_s_non_focal_genera == 0, axis=0)
        s_by_s_focal_genus = s_by_s_focal_genus[:,~idx_to_remove]
        s_by_s_non_focal_genera = s_by_s_non_focal_genera[:,~idx_to_remove]

        diversity_focal_genus = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_focal_genus)
        diversity_non_focal_genera = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_non_focal_genera)

        idx_diversity_to_keep = (diversity_focal_genus>0) & (diversity_non_focal_genera>0)
        diversity_focal_genus = diversity_focal_genus[idx_diversity_to_keep]
        diversity_non_focal_genera = diversity_non_focal_genera[idx_diversity_to_keep]

        if len(diversity_focal_genus) < 20:
            continue

        slope, intercept, r_value, p_value, std_err = stats.linregress(diversity_non_focal_genera, diversity_focal_genus)
        # get indices for null genera coarse graining
        unique_non_focal_genera, counts_non_focal_genera = numpy.unique(non_focal_genera_all, return_counts=True)
        coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_non_focal_genera))[:-1]

        s_by_s_copy = numpy.copy(s_by_s)
        slope_null_all = []
        while len(slope_null_all) < iter:

            numpy.random.shuffle(s_by_s_copy)
            # each element is the number of ASVS belonging to the non-focal genus in a given site
            s_by_s_focal_null = s_by_s_copy[:n_asvs_focal_genus,:]
            s_by_s_non_focal_null = s_by_s_copy[n_asvs_focal_genus:,:]

            s_by_s_non_focal_coarse_null = numpy.add.reduceat(s_by_s_non_focal_null, coarse_grain_by_genus_idx, axis=0)

            # remove sites where there are no observations in either focal or non-focal
            idx_to_remove_null = numpy.all(s_by_s_focal_null == 0, axis=0) | numpy.all(s_by_s_non_focal_coarse_null == 0, axis=0)
            s_by_s_focal_null = s_by_s_focal_null[:,~idx_to_remove_null]
            s_by_s_non_focal_coarse_null = s_by_s_non_focal_coarse_null[:,~idx_to_remove_null]

            diversity_focal_genus_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_focal_null)
            diversity_non_focal_genera_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_non_focal_coarse_null)

            idx_diversity_to_keep_null = (diversity_focal_genus_null>0) & (diversity_non_focal_genera_null>0)
            diversity_focal_genus_null = diversity_focal_genus_null[idx_diversity_to_keep_null]
            diversity_non_focal_genera_null = diversity_non_focal_genera_null[idx_diversity_to_keep_null]

            # not enough
            if len(diversity_non_focal_genera_null) < len(diversity_non_focal_genera):
                continue

            # randomly sample indexes
            idx_subsample = numpy.random.choice(numpy.arange(len(diversity_focal_genus_null)), size=len(diversity_non_focal_genera), replace=False)

            diversity_focal_genus_null_subsample = diversity_focal_genus_null[idx_subsample]
            diversity_non_focal_genera_null_subsample = diversity_non_focal_genera_null[idx_subsample]

            slope_null, intercept_null, r_value_null, p_value_null, std_err_null = stats.linregress(diversity_non_focal_genera_null_subsample, diversity_focal_genus_null_subsample)
            slope_null_all.append(slope_null)


        slope_null_all = numpy.asarray(slope_null_all)
        slope_null_all = numpy.sort(slope_null_all)

        slope_null_median = numpy.median(slope_null_all)
        if slope >= slope_null_median:
            p_value = sum(slope_null_all > slope)/iter
        else:
            p_value = sum(slope_null_all < slope)/iter

        percentile = sum(slope_null_all < slope)/iter

        #p_value = sum(slope_null_all > slope)/iter
        lower_ci = slope_null_all[int(iter*0.025)]
        upper_ci = slope_null_all[int(iter*0.975)]

        dbd_slope_dict[rank][focal_genus] = {}
        dbd_slope_dict[rank][focal_genus]['slope'] = slope
        dbd_slope_dict[rank][focal_genus]['null_dist'] = slope_null_all.tolist()
        dbd_slope_dict[rank][focal_genus]['p_value'] = p_value
        dbd_slope_dict[rank][focal_genus]['percentile'] = percentile
        dbd_slope_dict[rank][focal_genus]['lower_ci'] = lower_ci
        dbd_slope_dict[rank][focal_genus]['upper_ci'] = upper_ci
        dbd_slope_dict[rank][focal_genus]['n_asvs_focal_genus'] = n_asvs_focal_genus
        dbd_slope_dict[rank][focal_genus]['standard_score'] =  (slope - numpy.mean(slope_null_all)) / numpy.std(slope_null_all)

        print(focal_genus, slope, lower_ci, upper_ci, p_value)



if rarefied == True:
    rarefied_label = '_rarefied'
else:
    rarefied_label = ''

slope_path = slope_emp_taxa_template % (environment.replace(' ', '_') + rarefied_label)
sys.stderr.write("Saving slope dictionary...\n")
with open(slope_path, 'wb') as handle:
    pickle.dump(dbd_slope_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
