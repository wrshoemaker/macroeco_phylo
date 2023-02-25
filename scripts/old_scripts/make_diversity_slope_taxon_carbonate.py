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



environment = str(sys.argv[1])
environment = environment.replace('_', ' ')

measure = str(sys.argv[2])

rarefied = True
iter = 1000


sad_path_template = config.data_directory + "sad_annotated/sad_annotated_%s.pickle"
slope_emp_taxa_template = config.data_directory + "dbd_slope_emp_taxa/dbd_slope_%s.pickle"




def load_sad_annotated_taxon_dict(environment, rarefied=False):

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    sad_path = sad_path_template % (environment.replace(' ', '_') + rarefied_label)

    with open(sad_path, 'rb') as handle:
        presence_absence_dict = pickle.load(handle)
    return presence_absence_dict




#sys.stderr.write("Getting site-by-species matrix...\n")


sys.stderr.write("Loading taxon dict for %s...\n" % environment)
sad_annotated_dict = load_sad_annotated_taxon_dict(environment, rarefied=rarefied)
taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

dbd_slope_dict = {}
for coarse_rank_idx, coarse_rank in enumerate(diversity_utils.taxa_ranks):

    sys.stderr.write("Starting %s level analysis...\n" % coarse_rank)
    dbd_slope_dict[coarse_rank] = {}

    if coarse_rank == 'genus':
        fine_rank = 'asv'
        all_fine = numpy.copy(taxa).tolist()
    else:
        fine_rank = diversity_utils.taxa_ranks[coarse_rank_idx-1]
        all_fine = list([sad_annotated_dict['taxa'][t][fine_rank] for t in taxa])

    sys.stderr.write("Estimating dbd slope...\n")

    set_fine = list(set(all_fine))
    all_fine = numpy.asarray(all_fine)
    set_fine =  numpy.asarray(set_fine)

    fine_coarse_dict = {}
    for t in taxa:

        if fine_rank == 'asv':
            fine_coarse_dict[t] = sad_annotated_dict['taxa'][t][coarse_rank]

        else:
            fine_rank_t = sad_annotated_dict['taxa'][t][fine_rank]
            if fine_rank_t not in fine_coarse_dict:
                fine_coarse_dict[fine_rank_t] = sad_annotated_dict['taxa'][t][coarse_rank]


    if fine_rank == 'asv':
        s_by_s_fine = numpy.copy(s_by_s)
        all_coarse = []
        for fine in set_fine:
            all_coarse.append(fine_coarse_dict[fine])
        all_coarse = numpy.asarray(all_coarse)

    else:
        sad_fine_all = []
        all_coarse = []
        for fine in set_fine:
            taxa_in_fine = []
            for t in taxa:

                if sad_annotated_dict['taxa'][t][fine_rank] == fine:
                    taxa_in_fine.append(t)

            # numpy.where() is a major bottleneck
            g_taxa_idx = numpy.asarray([numpy.where(taxa == taxa_in_fine_i)[0][0] for taxa_in_fine_i in taxa_in_fine])
            sad_fine_all.append(s_by_s[g_taxa_idx,:].sum(axis=0))
            all_coarse.append(fine_coarse_dict[fine])

        all_coarse = numpy.asarray(all_coarse)
        s_by_s_fine = numpy.stack(sad_fine_all, axis=0)


    # remove sites where there are no observations
    s_by_s_fine = s_by_s_fine[:,~(numpy.all(s_by_s_fine == 0, axis=0))]

    # sort s_by_s_fine
    # so we can use cumsum
    set_coarse = list(set(all_coarse))
    idx_to_sort = []
    counts_coarse = []
    for coarse in set_coarse:
        idx_coarse_i = numpy.where(all_coarse==coarse)[0].tolist()
        counts_coarse.append(len(idx_coarse_i))
        idx_to_sort.append(idx_coarse_i)

    idx_to_sort_flat =  list(itertools.chain(*idx_to_sort))
    idx_to_sort_flat = numpy.asarray(idx_to_sort_flat)
    s_by_s_fine = s_by_s_fine[idx_to_sort_flat,:]

    fine_rel_s_by_s =  s_by_s_fine/numpy.sum(s_by_s_fine, axis=0)
    coarse_idx = numpy.append([0], numpy.cumsum(counts_coarse))[:-1]
    coarse_rel_s_by_s = numpy.add.reduceat(fine_rel_s_by_s, coarse_idx, axis=0)

    # get index for fine for null
    unique_fine, counts_fine = numpy.unique(all_fine, return_counts=True)
    fine_idx = numpy.append([0], numpy.cumsum(counts_fine))[:-1]
    for focal_coarse_idx, focal_coarse in enumerate(set_coarse):

        # ignore coarse-grained taxa with less than five fine-grained taxa

        if counts_coarse[focal_coarse_idx] < 5:
            continue

        # all the fine-scale indices for the focal
        focal_coarse_s_by_s_idx = numpy.asarray(idx_to_sort[focal_coarse_idx])
        focal_fine_rel_s_by_s = fine_rel_s_by_s[focal_coarse_s_by_s_idx,:]
        nonfocal_coarse_rel_s_by_s = numpy.delete(coarse_rel_s_by_s, focal_coarse_idx, axis=0)

        n_fine_focal_coarse = len(focal_coarse_s_by_s_idx)

        if measure == 'richness':
            measure_coarse = numpy.sum(nonfocal_coarse_rel_s_by_s > 0, axis=0)
            measure_fine = numpy.sum(focal_fine_rel_s_by_s > 0, axis=0)

        else:
            measure_coarse = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, nonfocal_coarse_rel_s_by_s)
            measure_fine = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, focal_fine_rel_s_by_s)

        # remove sites where there are no observations in either focal or non-focal
        idx_to_remove = (measure_fine>0) | (measure_coarse>0)

        measure_coarse = measure_coarse[idx_to_remove]
        measure_fine = measure_fine[idx_to_remove]

        if len(measure_fine) < 10:
            continue

        slope, intercept, r_value, p_value, std_err = stats.linregress(measure_coarse, measure_fine)

        s_by_s_copy = numpy.copy(s_by_s)
        slope_null_all = []
        while len(slope_null_all) < iter:

            numpy.random.shuffle(s_by_s_copy)

            rel_s_by_s_null =  s_by_s_copy/numpy.sum(s_by_s_copy, axis=0)
            # each element is the number of ASVS belonging to the non-focal genus in a given site
            fine_rel_s_by_s_null = numpy.add.reduceat(rel_s_by_s_null, fine_idx, axis=0)
            coarse_rel_s_by_s_null = numpy.add.reduceat(fine_rel_s_by_s_null, coarse_idx, axis=0)

            focal_fine_rel_s_by_s_null = fine_rel_s_by_s_null[focal_coarse_s_by_s_idx,:]
            nonfocal_coarse_rel_s_by_s_null = numpy.delete(coarse_rel_s_by_s_null, focal_coarse_idx, axis=0)
            focal_fine_rel_s_by_s_null = fine_rel_s_by_s_null[focal_coarse_s_by_s_idx,:]

            if measure == 'richness':
                measure_coarse_null = numpy.sum(nonfocal_coarse_rel_s_by_s_null > 0, axis=0)
                measure_fine_null = numpy.sum(focal_fine_rel_s_by_s_null > 0, axis=0)

            else:
                measure_coarse_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, nonfocal_coarse_rel_s_by_s_null)
                measure_fine_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, focal_fine_rel_s_by_s_null)


            idx_to_remove_null = (measure_fine_null>0) | (measure_coarse_null>0)
            measure_coarse_null = measure_coarse_null[idx_to_remove_null]
            measure_fine_null = measure_fine_null[idx_to_remove_null]

            if len(measure_fine_null) < len(measure_fine):
                continue

            #focal_coarse_idx
            # randomly sample indexes
            idx_subsample = numpy.random.choice(numpy.arange(len(measure_coarse_null)), size=len(measure_coarse), replace=False)

            measure_coarse_null = measure_coarse_null[idx_subsample]
            measure_fine_null = measure_fine_null[idx_subsample]

            slope_null, intercept_null, r_value_null, p_value_null, std_err_null = stats.linregress(measure_coarse_null, measure_fine_null)
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

        print(coarse_rank, focal_coarse,  n_fine_focal_coarse, slope, lower_ci, upper_ci)

        dbd_slope_dict[coarse_rank][focal_coarse] = {}
        dbd_slope_dict[coarse_rank][focal_coarse]['slope'] = slope
        dbd_slope_dict[coarse_rank][focal_coarse]['null_dist'] = slope_null_all.tolist()
        dbd_slope_dict[coarse_rank][focal_coarse]['p_value'] = p_value
        dbd_slope_dict[coarse_rank][focal_coarse]['percentile'] = percentile
        dbd_slope_dict[coarse_rank][focal_coarse]['lower_ci'] = lower_ci
        dbd_slope_dict[coarse_rank][focal_coarse]['upper_ci'] = upper_ci
        dbd_slope_dict[coarse_rank][focal_coarse]['n_fine_focal_coarse'] = n_fine_focal_coarse
        dbd_slope_dict[coarse_rank][focal_coarse]['standard_score'] =  (slope - numpy.mean(slope_null_all)) / numpy.std(slope_null_all)



slope_path = slope_emp_taxa_template % (diversity_utils.get_label(environment, rarefied) + '_' + measure)
sys.stderr.write("Saving slope dictionary...\n")
with open(slope_path, 'wb') as handle:
    pickle.dump(dbd_slope_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
