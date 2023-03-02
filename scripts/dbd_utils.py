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




n_fna_characters = 80

numpy.random.seed(123456789)
random.seed(123456789)


#coarse_grained_phylo_template = config.data_directory + "coarse_grained_%s.pickle"


permutation_rescaled_abundance_taxon_emp_template = config.data_directory + "test_constrain_occupancy_taxon_%d.pickle"
permutation_rescaled_abundance_phylo_emp_template = config.data_directory + "test_constrain_occupancy_phylo_%d.pickle"



sad_taxon_path_template = config.data_directory + "sad_annotated_taxon/sad_annotated_%s.pickle"
sad_phylo_path_template = config.data_directory + "sad_annotated_phylo/sad_annotated_%s.pickle"

sad_taxon_no_subsampling_path_template = config.data_directory + "sad_annotated_taxon_no_subsampling/sad_annotated_%s.pickle"
sad_phylo_no_subsampling_path_template = config.data_directory + "sad_annotated_phylo_no_subsampling/sad_annotated_%s.pickle"



coarse_grained_tree_path_template = config.data_directory + "coarse_grained_tree/%s.pickle"
coarse_grained_tree_no_subsampling_path_template = config.data_directory + "coarse_grained_tree_no_subsampling/%s.pickle"

#slope_emp_taxa_template = config.data_directory + "dbd_slope_emp_taxa/dbd_slope_%s.pickle"
#slope_emp_taxa_template = config.data_directory + "dbd_slope_emp_taxa/dbd_slope_%s.pickle"

slope_emp_phylo_template = config.data_directory + "dbd_slope_emp_phylo/coarse_grained_%s.pickle"
slope_emp_taxa_template = config.data_directory + "dbd_slope_emp_taxa/dbd_slope_%s.pickle"

slope_resource_taxa_path = "%sdbd_slope_resource_taxa.pickle" % config.data_directory
slope_resource_phylo_path = "%sdbd_slope_resource_phylo.pickle" % config.data_directory

diversity_vs_diversity_taxon_template = config.data_directory + "diversity_vs_diversity_taxon%s.pickle"
diversity_vs_diversity_phylo_resource_template = config.data_directory + "diversity_vs_diversity_phylo_resource%s.pickle"


diversity_vs_diversity_phylo_emp_template = config.data_directory + "diversity_vs_diversity_phylo%s.pickle"


richness_dbd_dict_path = config.data_directory + "richness_dbd_dict.pickle"
diversity_dbd_dict_path = config.data_directory + "diversity_dbd_dict/diversity_dbd_%s_dict.pickle"

diversity_dbd_simulation_taxon_dict_path = config.data_directory + "diversity_dbd_simulation_taxon_dict.pickle"
diversity_dbd_simulation_phylo_dict_path = config.data_directory + "diversity_dbd_simulation_phylo_dict/diversity_dbd_simulation_phylo_%s_dict.pickle"


richness_diversity_prediction_taxon_dict_path = config.data_directory + "richness_diversity_prediction_taxon_dict.pickle"
richness_diversity_prediction_phylo_dict_path = config.data_directory + "richness_diversity_prediction_phylo_dict/richness_diversity_prediction_phylo_%s_dict.pickle"


cov_phylo_dict_path = config.data_directory + "cov_phylo_dict/cov_phylo_%s_dict.pickle"



diversity_slm_phylo_integral_dict_path = config.data_directory + "diversity_slm_phylo_integral_dict/diversity_slm_phylo_integral_dict_%s.pickle"
diversity_slm_dbd_integral_dict_path = config.data_directory + "diversity_slm_dbd_integral_dict/diversity_slm_dbd_integral_dict_%s.pickle"


diversity_slm_phylo_integral_otu_level_dict_path = config.data_directory + "diversity_slm_phylo_integral_otu_level_dict/diversity_slm_phylo_integral_otu_level_dict_%s.pickle"


svd_taxon_dict_path = config.data_directory + "svd_taxon/svd_dict_%s.pickle"
svd_phylo_dict_path = config.data_directory + "svd_phylo/svd_dict_%s.pickle"


richness_diversity_phylo_mean_simulation_path = config.data_directory + "richness_diversity_phylo_mean_simulation_dict/richness_diversity_phylo_mean_simulation_dict_%s.pickle"



def make_coarse_grained_tree_dict(environment='soil metagenome', rarefied=False):

    # makes a dictionary containing a list of lists of taxa in each collapsed node
    tree = tree_utils.get_emp_tree()
    tree_copy = tree.copy()

    #sys.stderr.write("Subsetting samples...\n")
    #samples = diversity_utils.subset_observations(environment=environment)

    # get samples that were already subsetted
    #sys.stderr.write("Loading taxon dict...\n")
    #sad_annotated_dict = load_sad_annotated_taxon_dict(environment, rarefied=rarefied)
    #samples = sad_annotated_dict['samples']

    #sys.stderr.write("Getting site-by-species matrix...\n")
    #SADs, taxonomy_names, samples_keep = diversity_utils.get_SADs(samples)
    #biom_table_subset_values, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples, rarefied)
    #taxonomy_names = numpy.asarray(taxonomy_names)

    sys.stderr.write("Loading SAD phylo dictionary...\n")
    sad_annotated_phylo_dict = load_sad_annotated_phylo_dict(environment, rarefied=False)
    taxonomy_names = list(sad_annotated_phylo_dict['taxa'].keys())
    # need to convert from numpy str to Python str to work with ete3
    taxonomy_names = [str(s) for s in taxonomy_names]

    sys.stderr.write("Collapsing tree by distance...\n")
    #taxonomy_names_union = list(set(itertools.chain(*taxonomy_names)))
    tree_subset = tree_utils.subset_tree(taxonomy_names, tree_copy)


    distances = numpy.logspace(-2, 0.1, num=20)
    distance_collapsed_dict = tree_utils.make_distance_collapsed_dict(tree_subset, distances)

    sys.stderr.write("Saving coarse-grained dictionary...\n")

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    coarse_grained_tree_path = coarse_grained_tree_path_template % (environment.replace(' ', '_') + rarefied_label)
    with open(coarse_grained_tree_path, 'wb') as handle:
        pickle.dump(distance_collapsed_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)





def make_coarse_grained_tree_no_subsampling_dict(environment='soil metagenome', rarefied=False):

    # makes a dictionary containing a list of lists of taxa in each collapsed node
    tree = tree_utils.get_emp_tree()
    tree_copy = tree.copy()

    sys.stderr.write("Loading SAD phylo dictionary...\n")
    sad_annotated_phylo_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied=False)
    taxonomy_names = list(sad_annotated_phylo_dict['taxa'].keys())
    # need to convert from numpy str to Python str to work with ete3
    taxonomy_names = [str(s) for s in taxonomy_names]

    sys.stderr.write("Collapsing tree by distance...\n")
    #taxonomy_names_union = list(set(itertools.chain(*taxonomy_names)))
    tree_subset = tree_utils.subset_tree(taxonomy_names, tree_copy)


    distances = numpy.logspace(-2, 0.1, num=20)
    distance_collapsed_dict = tree_utils.make_distance_collapsed_dict(tree_subset, distances)

    sys.stderr.write("Saving coarse-grained dictionary...\n")

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    coarse_grained_tree_path = coarse_grained_tree_no_subsampling_path_template % (environment.replace(' ', '_') + rarefied_label)
    with open(coarse_grained_tree_path, 'wb') as handle:
        pickle.dump(distance_collapsed_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)











def load_coarse_grained_tree_dict(environment='human gut metagenome', rarefied=False):

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    coarse_grained_tree_path = coarse_grained_tree_path_template % (environment.replace(' ', '_') + rarefied_label)

    with open(coarse_grained_tree_path, 'rb') as handle:
        distance_collapsed_dict = pickle.load(handle)
    return distance_collapsed_dict




def load_coarse_grained_tree_no_subsampling_dict(environment='human gut metagenome', rarefied=False):

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    coarse_grained_tree_path = coarse_grained_tree_no_subsampling_path_template % (environment.replace(' ', '_') + rarefied_label)

    with open(coarse_grained_tree_path, 'rb') as handle:
        distance_collapsed_dict = pickle.load(handle)
    return distance_collapsed_dict



def make_diversity_slope_phylo_dict(environment='human gut metagenome', iter = 1000):


    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]
    pres_abs_dict = load_sad_annotated_phylo_dict(environment, rarefied = rarefied)

    samples = pres_abs_dict['samples']
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))

    sys.stderr.write("Getting site-by-species matrix...\n")
    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])


    sad_annotated_phylo_dict = load_coarse_grained_tree_dict(environment=environment)

    distances = list(coarse_grained_tree_dict.keys())
    distances.sort()

    sys.stderr.write("Subsetting samples...\n")
    #samples = diversity_utils.subset_observations(environment=environment)

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

        sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))

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


            #if distance not in dbd_slope_dict:
            #    dbd_slope_dict[distance] = {}
            #    dbd_slope_dict[distance]['slope'] = []
            #    dbd_slope_dict[distance]['collapsed_nodes'] = []
            #    dbd_slope_dict[distance]['null_dist'] = []
            #    dbd_slope_dict[distance]['p_value'] = []
            #    dbd_slope_dict[distance]['percentile'] = []
            #    dbd_slope_dict[distance]['lower_ci'] = []
            #    dbd_slope_dict[distance]['upper_ci'] = []
            #    dbd_slope_dict[distance]['n_asvs_focal_genus'] = []
            #    dbd_slope_dict[distance]['standard_score'] = []

            slope_null_all = numpy.asarray(slope_null_all)
            slope_null_all = numpy.sort(slope_null_all)

            slope_null_median = numpy.median(slope_null_all)
            if slope >= slope_null_median:
                p_value = sum(slope_null_all > slope)/iter
            else:
                p_value = sum(slope_null_all < slope)/iter

            percentile = sum(slope_null_all < slope)/iter

            lower_ci = slope_null_all[int(iter*0.025)]
            upper_ci = slope_null_all[int(iter*0.975)]

            #dbd_slope_dict[distance]['slope'].append(slope)
            #dbd_slope_dict[distance]['collapsed_nodes'].append(coarse_grained_list[i])
            #dbd_slope_dict[distance]['null_dist'].append(slope_null_all.tolist())
            #dbd_slope_dict[distance]['p_value'].append(p_value)
            #dbd_slope_dict[distance]['percentile'].append(percentile)
            #dbd_slope_dict[distance]['lower_ci'].append(lower_ci)
            #dbd_slope_dict[distance]['upper_ci'].append(upper_ci)
            #dbd_slope_dict[distance]['n_asvs_focal_genus'].append(n_asvs_focal_genus)
            #dbd_slope_dict[distance]['standard_score'].append((slope - numpy.mean(slope_null_all)) / numpy.std(slope_null_all))

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


        slope_path = slope_emp_phylo_template % (environment.replace(' ', '_') + '_phylo_' + str(distance))
        sys.stderr.write("Saving slope dictionary...\n")
        with open(slope_path, 'wb') as handle:
            pickle.dump(dbd_slope_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




def load_diversity_slope_phylo_dict(environment):

    slope_path = slope_emp_phylo_template % (environment.replace(' ', '_') + '_phylo')

    with open(slope_path, 'rb') as handle:
        dbd_slope_dict = pickle.load(handle)
    return dbd_slope_dict





def get_subsampled_sad(sad, replace=True, n_subsample = 10):

    s = len(sad)

    relative_sad = sad / sum(sad)

    if replace == True:

        # sampling WITH replacement
        # write function to sample SAD without replacement
        individuals_subsample = numpy.random.choice(range(s), size=n_subsample, replace=True, p=relative_sad)

        species_subsample, sad_subsample = numpy.unique(individuals_subsample, return_counts=True)


    else:

        all_individuals = []
        for idx, i in enumerate(range(len(sad))):
            count_i = int(sad[idx])
            if count_i == 0:
                continue
            all_individuals.extend( [i]* count_i)

        subsample_individuals = random.sample(all_individuals, n_subsample)
        species_counts = Counter(subsample_individuals)

        sad_subsample = []
        # put zeros back in so we can map to ASV labels
        for i in range(len(sad)):

            if i in species_counts:
                sad_subsample.append(species_counts[i])

            else:
                sad_subsample.append(0)

        #sad_subsample = list(species_counts.values())
        #sad_subsample.sort(reverse=True)
        #sad_subsample = numpy.asarray(sad_subsample)

        sad_subsample = numpy.asarray(sad_subsample)
    return sad_subsample







def make_sad_annotated_taxon_dict(environment='human gut metagenome', n_subsample=100, rarefied=False, n_reads=10**3):

    # subsamples all the data

    sys.stderr.write("Subsetting samples...\n")
    samples = diversity_utils.subset_observations(environment=environment)

    sys.stderr.write("Getting SADs...\n")
    SADs, taxonomy_names, samples_to_keep = diversity_utils.get_SADs(samples, rarefied=rarefied)

    sys.stderr.write("Making SAD dictionary...\n")
    # D_5
    taxonomy_names_flat = numpy.concatenate(taxonomy_names).ravel()
    taxonomy_names_unique = list(set(taxonomy_names_flat.tolist()))

    afd_dict = {}
    #afd_dict['samples'] = samples
    afd_dict['taxa'] = {}

    for line in open("%semp/otu_info/silva_123/taxonomy/consensus_taxonomy_7_levels.txt" % config.data_directory, 'r'):
        line_split = line.strip().split('\t')
        taxon = line_split[0]
        phylum_ = line_split[1].split(';')[1].split('D_1__')[-1]
        class_ = line_split[1].split(';')[2].split('D_2__')[-1]
        order_ = line_split[1].split(';')[3].split('D_3__')[-1]
        family_ = line_split[1].split(';')[4].split('D_4__')[-1]
        genus_ = line_split[1].split(';')[5].split('D_5__')[-1]

        if genus_ in diversity_utils.genera_to_ignore:
            continue

        if '[' in genus_:
            genus_ = re.findall('\[(.*?)\]', genus_)[0]

        # lower case to make finding bugs easier
        phylum_ = phylum_.lower()
        class_ = class_.lower()
        order_ = order_.lower()
        family_ = family_.lower()
        genus_ = genus_.lower()

        if any(s in phylum_ for s in diversity_utils.phylum_partial_match_to_ignore):
            continue

        if any(s in class_ for s in diversity_utils.class_partial_match_to_ignore):
            continue

        if any(s in order_ for s in diversity_utils.order_partial_match_to_ignore):
            continue

        if any(s in family_ for s in diversity_utils.family_partial_match_to_ignore):
            continue

        if any(s in genus_ for s in diversity_utils.genera_partial_match_to_ignore):
            continue

        genus_split = genus_.split(' ')
        family_split = family_.split(' ')
        if len(genus_split) == 2:
            if genus_split[1].isdigit() == True:
                genus_ = genus_split[0]

        if 'sp.' in genus_:
            genus_ = genus_split[genus_split.index('sp.') -1]

        if (len(genus_split) > 1) and ('candidatus' not in genus_):
            genus_ = genus_split[0]

        if (len(family_split) > 1) and ('candidatus' not in genus_):
            family_ = family_split[0]

        if family_ == 'clostridiaceae 1':
            family_ = 'clostridiaceae'

        if family_ == 'acidobacteriaceae (subgroup 1)':
            family_ = 'acidobacteriaceae'

        if taxon in taxonomy_names_unique:
            afd_dict['taxa'][taxon] = {}
            afd_dict['taxa'][taxon]['phylum'] = phylum_
            afd_dict['taxa'][taxon]['class'] = class_
            afd_dict['taxa'][taxon]['order'] = order_
            afd_dict['taxa'][taxon]['family'] = family_
            afd_dict['taxa'][taxon]['genus'] = genus_
            afd_dict['taxa'][taxon]['abundance'] = []


    taxa_to_keep = list(afd_dict['taxa'].keys())
    for i in range(len(SADs)):

        SAD_i = SADs[i]
        taxonomy_names_i = taxonomy_names[i]

        for t in taxa_to_keep:
            if t in taxonomy_names_i:
                afd_dict['taxa'][t]['abundance'].append(SAD_i[numpy.where(taxonomy_names_i==t)[0][0]])
            else:
                afd_dict['taxa'][t]['abundance'].append(float(0))

    # re-rarify
    #turn dictionary into matrix
    all_taxa = numpy.asarray(list(afd_dict['taxa'].keys()))
    s_by_s = numpy.asarray([afd_dict['taxa'][t]['abundance'] for t in all_taxa])
    # remove low coverage sites
    idx_samples_to_keep = (s_by_s.sum(axis=0) > n_reads)
    samples_to_keep = numpy.asarray(samples_to_keep)
    samples_to_keep = samples_to_keep[idx_samples_to_keep]
    s_by_s = s_by_s[:,idx_samples_to_keep]

    # subsample without replacement
    subsample_sad_all = []
    for sad in s_by_s.T:
        subsample_sad = get_subsampled_sad(sad, replace=False, n_subsample = n_reads)
        subsample_sad_all.append(subsample_sad)

    subsample_s_by_s = numpy.array(subsample_sad_all)

    # sample sites
    numpy.random.shuffle(subsample_s_by_s)
    idx_samples_to_keep_subset = numpy.random.choice(len(samples_to_keep), size=n_subsample, replace=False)
    samples_to_keep_subset = samples_to_keep[idx_samples_to_keep_subset]
    subsample_s_by_s_subset = subsample_s_by_s[idx_samples_to_keep_subset,:]

    # remove absent species, species are columns
    idx_no_zeros = ~(numpy.all(subsample_s_by_s_subset == 0, axis=0))
    taxa_to_keep = all_taxa[idx_no_zeros]
    subsample_s_by_s_subset = subsample_s_by_s_subset[:,idx_no_zeros]

    subsample_afd_dict = {}
    subsample_afd_dict['samples'] = samples_to_keep_subset
    subsample_afd_dict['taxa'] = {}

    for t_idx, t in enumerate(taxa_to_keep):

        subsample_afd_dict['taxa'][t] = {}
        subsample_afd_dict['taxa'][t]['phylum'] = afd_dict['taxa'][t]['phylum']
        subsample_afd_dict['taxa'][t]['class'] = afd_dict['taxa'][t]['class']
        subsample_afd_dict['taxa'][t]['order'] = afd_dict['taxa'][t]['order']
        subsample_afd_dict['taxa'][t]['family'] = afd_dict['taxa'][t]['family']
        subsample_afd_dict['taxa'][t]['genus'] = afd_dict['taxa'][t]['genus']
        subsample_afd_dict['taxa'][t]['abundance'] = subsample_s_by_s_subset[:,t_idx].tolist()

    # get number of taxa without annotation
    n_excluded_taxa = len(taxonomy_names_unique) - len(afd_dict['taxa'].keys())
    subsample_afd_dict['n_taxa_without_annotation'] = n_excluded_taxa

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    sad_path = sad_taxon_path_template % (environment.replace(' ', '_') + rarefied_label)

    sys.stderr.write("Saving SAD dictionary...\n")
    with open(sad_path, 'wb') as handle:
        pickle.dump(subsample_afd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)









def load_sad_annotated_taxon_dict(environment, rarefied=False):

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    sad_path = sad_taxon_path_template % (environment.replace(' ', '_') + rarefied_label)

    with open(sad_path, 'rb') as handle:
        presence_absence_dict = pickle.load(handle)
    return presence_absence_dict





def make_sad_annotated_no_subsampling_taxon_dict(environment):

    sys.stderr.write("Getting samples...\n")

    sad_annotated_dict = load_sad_annotated_taxon_dict(environment)
    samples = sad_annotated_dict['samples']

    sys.stderr.write("Getting SADs...\n")
    SADs, taxonomy_names, samples_to_keep = diversity_utils.get_SADs(samples, rarefied=False)

    sys.stderr.write("Making SAD dictionary...\n")
    # D_5
    taxonomy_names_flat = numpy.concatenate(taxonomy_names).ravel()
    taxonomy_names_unique = list(set(taxonomy_names_flat.tolist()))


    afd_dict = {}
    afd_dict['samples'] = samples_to_keep
    afd_dict['taxa'] = {}

    for line in open("%semp/otu_info/silva_123/taxonomy/consensus_taxonomy_7_levels.txt" % config.data_directory, 'r'):
        line_split = line.strip().split('\t')
        taxon = line_split[0]
        phylum_ = line_split[1].split(';')[1].split('D_1__')[-1]
        class_ = line_split[1].split(';')[2].split('D_2__')[-1]
        order_ = line_split[1].split(';')[3].split('D_3__')[-1]
        family_ = line_split[1].split(';')[4].split('D_4__')[-1]
        genus_ = line_split[1].split(';')[5].split('D_5__')[-1]

        if genus_ in diversity_utils.genera_to_ignore:
            continue

        if '[' in genus_:
            genus_ = re.findall('\[(.*?)\]', genus_)[0]

        # lower case to make finding bugs easier
        phylum_ = phylum_.lower()
        class_ = class_.lower()
        order_ = order_.lower()
        family_ = family_.lower()
        genus_ = genus_.lower()

        if any(s in phylum_ for s in diversity_utils.phylum_partial_match_to_ignore):
            continue

        if any(s in class_ for s in diversity_utils.class_partial_match_to_ignore):
            continue

        if any(s in order_ for s in diversity_utils.order_partial_match_to_ignore):
            continue

        if any(s in family_ for s in diversity_utils.family_partial_match_to_ignore):
            continue

        if any(s in genus_ for s in diversity_utils.genera_partial_match_to_ignore):
            continue

        genus_split = genus_.split(' ')
        family_split = family_.split(' ')
        if len(genus_split) == 2:
            if genus_split[1].isdigit() == True:
                genus_ = genus_split[0]

        if 'sp.' in genus_:
            genus_ = genus_split[genus_split.index('sp.') -1]

        if (len(genus_split) > 1) and ('candidatus' not in genus_):
            genus_ = genus_split[0]

        if (len(family_split) > 1) and ('candidatus' not in genus_):
            family_ = family_split[0]

        if family_ == 'clostridiaceae 1':
            family_ = 'clostridiaceae'

        if family_ == 'acidobacteriaceae (subgroup 1)':
            family_ = 'acidobacteriaceae'

        if taxon in taxonomy_names_unique:
            afd_dict['taxa'][taxon] = {}
            afd_dict['taxa'][taxon]['phylum'] = phylum_
            afd_dict['taxa'][taxon]['class'] = class_
            afd_dict['taxa'][taxon]['order'] = order_
            afd_dict['taxa'][taxon]['family'] = family_
            afd_dict['taxa'][taxon]['genus'] = genus_
            afd_dict['taxa'][taxon]['abundance'] = []


    taxa_to_keep = list(afd_dict['taxa'].keys())
    for i in range(len(SADs)):

        SAD_i = SADs[i]
        taxonomy_names_i = taxonomy_names[i]

        for t in taxa_to_keep:
            if t in taxonomy_names_i:
                afd_dict['taxa'][t]['abundance'].append(SAD_i[numpy.where(taxonomy_names_i==t)[0][0]])
            else:
                afd_dict['taxa'][t]['abundance'].append(float(0))


    sad_path = sad_taxon_no_subsampling_path_template % (environment.replace(' ', '_'))

    sys.stderr.write("Saving SAD dictionary...\n")
    with open(sad_path, 'wb') as handle:
        pickle.dump(afd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)






def load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied=False):

    sad_path = sad_taxon_no_subsampling_path_template % (environment.replace(' ', '_'))

    with open(sad_path, 'rb') as handle:
        presence_absence_dict = pickle.load(handle)
    return presence_absence_dict









def make_sad_annotated_no_subsampling_phylo_dict(environment):

    sys.stderr.write("Getting samples...\n")

    sad_annotated_dict = load_sad_annotated_phylo_dict(environment)
    samples = sad_annotated_dict['samples']

    sys.stderr.write("Getting SADs...\n")
    SADs, taxonomy_names, samples_to_keep = diversity_utils.get_SADs(samples, rarefied=False)

    sys.stderr.write("Making SAD dictionary...\n")

    afd_dict = {}
    afd_dict['samples'] = samples_to_keep
    afd_dict['taxa'] = {}

    # make sure OTU is also in phylogeny
    tree = tree_utils.get_emp_tree()
    tree_copy = tree.copy()
    leaf_names = tree_copy.get_leaf_names()



    taxonomy_names_all = numpy.concatenate(taxonomy_names)
    taxonomy_names_set_old = list(set(taxonomy_names_all.tolist()))
    taxonomy_names_set = list(set(taxonomy_names_all.tolist()) & set(leaf_names))
    print(len(taxonomy_names_set_old), len(taxonomy_names_set))

    # initialize abundance dictionary
    for t in taxonomy_names_set:
        afd_dict['taxa'][t] = {}
        afd_dict['taxa'][t]['abundance'] = []


    #taxa_to_keep = list(afd_dict['taxa'].keys())
    for i in range(len(SADs)):
        SAD_i = SADs[i]
        taxonomy_names_i = numpy.asarray(taxonomy_names[i])
        for t in taxonomy_names_set:
            if t in taxonomy_names_i:
                afd_dict['taxa'][t]['abundance'].append(SAD_i[numpy.where(taxonomy_names_i==t)[0][0]])
            else:
                afd_dict['taxa'][t]['abundance'].append(float(0))



    sad_path = sad_phylo_no_subsampling_path_template % (environment.replace(' ', '_'))

    sys.stderr.write("Saving SAD dictionary...\n")
    with open(sad_path, 'wb') as handle:
        pickle.dump(afd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)





def load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied=False):

    sad_path = sad_phylo_no_subsampling_path_template % (environment.replace(' ', '_'))

    with open(sad_path, 'rb') as handle:
        presence_absence_dict = pickle.load(handle)
    return presence_absence_dict





def make_sad_annotated_phylo_dict(environment='human gut metagenome', rarefied=False, n_reads=10**3):

    sad_annotated_taxon_dict = load_sad_annotated_taxon_dict(environment, rarefied=rarefied)
    samples = sad_annotated_taxon_dict['samples']

    sys.stderr.write("Getting SADs...\n")
    SADs, taxonomy_names, samples_to_keep = diversity_utils.get_SADs(samples, rarefied=rarefied)

    tree = tree_utils.get_emp_tree()
    taxa_to_keep = tree.get_leaf_names()

    afd_dict = {}
    afd_dict['taxa'] = {}
    taxonomy_names_all = []
    for t in taxonomy_names:
        taxonomy_names_all.extend(t.tolist())
    taxonomy_names_set = list(set(taxonomy_names_all))
    taxonomy_names_set_to_keep = []
    for t in taxonomy_names_set:

        # make sure ASV is in the phylogeny
        if t in taxa_to_keep:
            afd_dict['taxa'][t] = {}
            afd_dict['taxa'][t]['abundance'] = []
            taxonomy_names_set_to_keep.append(t)


    for i in range(len(SADs)):

        SAD_i = SADs[i]
        taxonomy_names_i = taxonomy_names[i]

        for t in taxonomy_names_set_to_keep:
             if t in taxonomy_names_i:
                 afd_dict['taxa'][t]['abundance'].append(SAD_i[numpy.where(taxonomy_names_i==t)[0][0]])
             else:
                 afd_dict['taxa'][t]['abundance'].append(float(0))


    all_taxa = numpy.asarray(list(afd_dict['taxa'].keys()))
    s_by_s = numpy.asarray([afd_dict['taxa'][t]['abundance'] for t in all_taxa])

    # subsample without replacement
    subsample_sad_all = []
    for sad in s_by_s.T:
        subsample_sad = get_subsampled_sad(sad, replace=False, n_subsample = n_reads)
        subsample_sad_all.append(subsample_sad)

    subsample_s_by_s = numpy.array(subsample_sad_all)

    # remove absent species, species are columns
    idx_no_zeros = ~(numpy.all(subsample_s_by_s == 0, axis=0))
    taxa_to_keep = all_taxa[idx_no_zeros]
    subsample_s_by_s = subsample_s_by_s[:,idx_no_zeros]

    subsample_afd_dict = {}
    subsample_afd_dict['samples'] = samples
    subsample_afd_dict['taxa'] = {}

    for t_idx, t in enumerate(taxa_to_keep):
        subsample_afd_dict['taxa'][t] = {}
        subsample_afd_dict['taxa'][t]['abundance'] = subsample_s_by_s[:,t_idx].tolist()

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    sad_path = sad_phylo_path_template % (environment.replace(' ', '_') + rarefied_label)

    sys.stderr.write("Saving SAD dictionary...\n")
    with open(sad_path, 'wb') as handle:
        pickle.dump(subsample_afd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_sad_annotated_phylo_dict(environment, rarefied=False):

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    sad_path = sad_phylo_path_template % (environment.replace(' ', '_') + rarefied_label)

    with open(sad_path, 'rb') as handle:
        presence_absence_dict = pickle.load(handle)
    return presence_absence_dict










def make_diversity_slope_taxon_dict_old(environment, iter=1000, rarefied=False):

    sys.stderr.write("Loading taxon dict...\n")
    sad_annotated_dict = load_sad_annotated_taxon_dict(environment, rarefied=rarefied)

    dbd_slope_dict = {}
    for rank in diversity_utils.taxa_ranks:
        dbd_slope_dict[rank] = {}

        sys.stderr.write("Starting %s level analysis...\n" % rank)

        sys.stderr.write("Estimating dbd slope...\n")
        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
        all_genera_idx = numpy.arange(len(all_genera))

        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

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
        # remove sites where there are no observations
        s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]

        for focal_genus in all_genera:

            #presence_absence_all = []
            sad_all = []
            s_by_s_focal_genus = []
            sad_non_focal_genus = []
            non_focal_genera_all = []

            for t in taxa:
                genus_t = sad_annotated_dict['taxa'][t][rank]

                if genus_t == focal_genus:
                    s_by_s_focal_genus.append(sad_annotated_dict['taxa'][t]['abundance'])

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

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    slope_path = slope_emp_taxa_template % (environment.replace(' ', '_') + rarefied_label)
    sys.stderr.write("Saving slope dictionary...\n")
    with open(slope_path, 'wb') as handle:
        pickle.dump(dbd_slope_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)





def make_diversity_slope_taxon_dict(environment, measure = 'richness', iter = 1000, rarefied=True):

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
            focal_s_by_s_coarse_idx = numpy.asarray(idx_to_sort[focal_coarse_idx])
            focal_fine_rel_s_by_s = fine_rel_s_by_s[focal_s_by_s_coarse_idx,:]
            nonfocal_coarse_rel_s_by_s = numpy.delete(coarse_rel_s_by_s, focal_coarse_idx, axis=0)

            n_fine_focal_coarse = len(focal_s_by_s_coarse_idx)

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

                focal_fine_rel_s_by_s_null = fine_rel_s_by_s_null[focal_s_by_s_coarse_idx,:]
                nonfocal_coarse_rel_s_by_s_null = numpy.delete(coarse_rel_s_by_s_null, focal_coarse_idx, axis=0)
                focal_fine_rel_s_by_s_null = fine_rel_s_by_s_null[focal_s_by_s_coarse_idx,:]

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










def load_diversity_slope_taxon_dict(environment, measure, rarefied=True):

    slope_path = slope_emp_taxa_template % (diversity_utils.get_label(environment, rarefied) + '_' + measure)

    with open(slope_path, 'rb') as handle:
        dbd_slope_dict = pickle.load(handle)
    return dbd_slope_dict






''' dalballo resource analyses '''


def get_dalbello_16_samples_and_metadata():

    #sample_list = list(range(11941650, 11942242 + 1)) + list(range(11942243, 11942410 + 1))
    sample_to_keep = []
    metadata_to_keep = []
    sra_to_sample_dict = {}
    #for s in ['SRR13989233']:
    for s in open('%sdalbello/SraAccList.txt' % config.data_directory, 'r'):

        #s = 'SAMN%d' % s
        s = s.strip()

        proc = subprocess.Popen("/Users/williamrshoemaker/edirect/esearch -db sra -query '%s'  | /Users/williamrshoemaker/edirect/efetch -format runinfo" % s, stdout=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()

        out = out.decode("utf-8")
        out = out.strip()
        out_split = out.split(',')
        if len(out_split) <= 1:
            continue
        sample = out_split[57].split('-')[0]
        sra_to_sample_dict[s] = sample

        metadata_to_keep.append(out)
        sample_to_keep.append(s)


    metadata_file = open('%sdalbello/sra_metadata.txt' % config.data_directory, "w")
    metadata_file.write("\n".join(metadata_to_keep))
    metadata_file.close()

    sample_file = open('%sdalbello/sra_samples.txt' % config.data_directory, "w")
    sample_file.write("\n".join(sample_to_keep))
    sample_file.close()

    dict_path = "%sdalbello/sra_sample.pickle" % config.data_directory

    with open(dict_path, 'wb') as handle:
        pickle.dump(sra_to_sample_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_sra_to_sample_dict():

    dict_path = "%sdalbello/sra_sample.pickle" % config.data_directory

    with open(dict_path, 'rb') as handle:
        sample_dict = pickle.load(handle)
    return sample_dict


def make_fasta_from_dada2_output():

    dada2_path = "%sdalbello/seqtab-nochim.txt" % config.data_directory
    fasta_path = "%sdalbello/seqtab-nochim.fna" % config.data_directory
    fasta_file = open(fasta_path, 'w')

    for line_idx, line in enumerate(open(dada2_path, 'r')):

        if line_idx == 0:
            continue

        line = line.strip()

        sequence = line.split('\t')[0]
        fasta_file.write('>%s\n' % sequence)

        for i in range(0, len(sequence), n_fna_characters):
            sequence_i = sequence[i : i + n_fna_characters]
            fasta_file.write('%s\n' % sequence_i)
        fasta_file.write('\n')

    fasta_file.close()




class classFASTA:

    # class to load FASTA file

    def __init__(self, fileFASTA):
        self.fileFASTA = fileFASTA

    def readFASTA(self):
        '''Checks for fasta by file extension'''
        file_lower = self.fileFASTA.lower()
        '''Check for three most common fasta file extensions'''
        if file_lower.endswith('.txt') or file_lower.endswith('.fa') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.fna') or \
        file_lower.endswith('.fasta') or file_lower.endswith('.frn') or \
        file_lower.endswith('.faa') or file_lower.endswith('.ffn'):
            with open(self.fileFASTA, "r") as f:
                return self.ParseFASTA(f)
        else:
            print("Not in FASTA format.")

    def ParseFASTA(self, fileFASTA):
        '''Gets the sequence name and sequence from a FASTA formatted file'''
        fasta_list=[]
        for line in fileFASTA:
            if line[0] == '>':
                try:
                    fasta_list.append(current_dna)
            	#pass if an error comes up
                except UnboundLocalError:
                    #print "Inproper file format."
                    pass
                current_dna = [line.lstrip('>').rstrip('\n'),'']
            else:
                current_dna[1] += "".join(line.split())
        fasta_list.append(current_dna)
        '''Returns fasa as nested list, containing line identifier \
            and sequence'''
        return fasta_list




def clean_alignment(muscle_path, muscle_clean_path, min_n_sites=100, max_fraction_empty=0.8):

    # removes all sites where the fraction of empty bases across ASVs is greater than max_fraction_empty (putatively uninformative)
    # removes a sequencies with fewer than min_n_sites informative sites

    frn_aligned = classFASTA(muscle_path).readFASTA()

    n = len(frn_aligned)

    frn_aligned_seqs = [x[1] for x in frn_aligned]
    frn_aligned_seqs_names = [x[0] for x in frn_aligned]

    frns = []
    for site in zip(*frn_aligned_seqs):

        fraction_empty = site.count('-')/n

        if fraction_empty > max_fraction_empty:
            continue

        # skip site if it is uninformative
        if len(set([s for s in site if s != '-'])) == 1:
            continue

        frns.append(site)

    if len(frns) < min_n_sites:
        exit()

    clean_sites_list = zip(*frns)

    # skip if there are too few sites
    #if len(list(clean_sites_list)[0]) < min_n_sites:
    #    continue

    frn_aligned_clean = open(muscle_clean_path, 'w')

    for clean_sites_idx, clean_sites in enumerate(clean_sites_list):
        clean_sites_species = frn_aligned_seqs_names[clean_sites_idx]
        clean_sites_seq = "".join(clean_sites)

        frn_aligned_clean.write('>%s\n' % clean_sites_species)

        clean_sites_seq_split = [clean_sites_seq[i:i+n_fna_characters] for i in range(0, len(clean_sites_seq), n_fna_characters)]

        for seq in clean_sites_seq_split:

            frn_aligned_clean.write('%s\n' % seq)

        frn_aligned_clean.write('\n')


    frn_aligned_clean.close()





def get_sample_dict_debello():

    metadata_path = '%sdalbello/metadata.txt' % config.data_directory

    file = open(metadata_path, 'r')
    header = file.readline()
    header = header.strip()
    header = header.split('\t')

    carbon_array = numpy.asarray(header[19:])

    sample_dict = {}
    for line in file:
        line = line.strip()
        line = line.split('\t')

        if (line[9] == 'NA') or (line[9] == 'H2O') or (line[9] == 'M9'):
            continue

        plate_and_well = line[0] + line[1]

        n_resournces = int(line[10])

        if n_resournces not in sample_dict:
            sample_dict[n_resournces] = {}

        carbon_idx = (numpy.asarray(line[19:]) == 'Y')

        carbon_line = carbon_array[carbon_idx]
        carbon_set = tuple(carbon_line.tolist())

        if carbon_set not in sample_dict[n_resournces]:
            sample_dict[n_resournces][carbon_set] = []



        sample_dict[n_resournces][carbon_set].append(plate_and_well)

    file.close()

    return sample_dict



def get_taxon_dict_dalballo():

    taxon_dict = {}

    metadata_path = '%sdalbello/seqtab-nochim-taxa.txt' % config.data_directory

    file = open(metadata_path, 'r')
    header = file.readline()
    header = header.strip()
    header = header.split('\t')

    for line in file:
        line = line.strip()
        line = line.split('\t')

        asv = line[0]

        taxon_dict[asv] = {}

        taxon_dict[asv]['kingdom'] = line[1]
        taxon_dict[asv]['phylum'] = line[2]
        taxon_dict[asv]['class'] = line[3]
        taxon_dict[asv]['order'] = line[4]
        taxon_dict[asv]['family'] = line[5]
        taxon_dict[asv]['genus'] = line[6]
        taxon_dict[asv]['species'] = line[7]

    file.close()

    return taxon_dict





def get_s_by_s_dalballo():

    sample_dict = get_sample_dict_debello()

    sra_to_sample_dict = load_sra_to_sample_dict()

    data_path = '%sdalbello/seqtab-nochim.txt' % config.data_directory

    file = open(data_path, 'r')
    header = file.readline()
    header = header.strip()
    header = header.split('\t')

    sample_names = numpy.asarray([sra_to_sample_dict[h] for h in header])
    species = []
    species_abundances_all = []

    for line in file:
        line = line.strip()
        line = line.split('\t')

        species.append(line[0])
        species_abundances = [int(s) for s in line[1:]]
        species_abundances_all.append(species_abundances)

    s_by_s = numpy.asarray(species_abundances_all)

    return s_by_s, species, sample_names





def make_diversity_slope_taxon_dalballo_dict(iter=1000):

    sample_dict = get_sample_dict_debello()
    taxon_dict = get_taxon_dict_dalballo()
    sra_to_sample_dict = load_sra_to_sample_dict()

    n_carbons = list(sample_dict.keys())
    n_carbons.sort()

    s_by_s, species, sample_names = get_s_by_s_dalballo()
    species = numpy.asarray(species)

    dbd_slope_dict = {}

    n_replicates = 9

    for n in n_carbons:

        dbd_slope_dict[n] = {}

        # get all samples with n resources
        flat_array = numpy.asarray(numpy.concatenate(list(sample_dict[n].values())).flat)
        idx_sample_range = numpy.asarray([numpy.where(sample_names==s)[0][0] for s in flat_array])

        # sub-sample communities
        idx_sample_n = numpy.random.choice(idx_sample_range, size=n_replicates, replace=False)
        s_by_s_subsample = s_by_s[:,idx_sample_n]
        samples_names_subsample = sample_names[idx_sample_n].tolist()
        dbd_slope_dict[n]['samples'] = samples_names_subsample
        dbd_slope_dict[n]['ranks'] = {}

        for rank in diversity_utils.taxa_ranks:

            # make copy for each rank
            s_by_s_n = numpy.copy(s_by_s_subsample)

            dbd_slope_dict[n]['ranks'][rank] = {}

            taxon_all = numpy.asarray([taxon_dict[s][rank] for s in species])

            # remove species with no annotation
            idx_to_keep = (taxon_all!='NA')
            s_by_s_n = s_by_s_n[idx_to_keep,:]
            species_n = species[idx_to_keep]

            # remove sites with all zeros
            s_by_s_n = s_by_s_n[:,~(numpy.all(s_by_s_n == 0, axis=0))]
            # remove species with all zeros
            idx_species_to_keep = ~(numpy.all(s_by_s_n == 0, axis=1))
            s_by_s_n = s_by_s_n[idx_species_to_keep,:]

            species_n = species_n[idx_species_to_keep]
            taxon_n_all = numpy.asarray([taxon_dict[s][rank] for s in species_n])
            focal_all = list(set(taxon_n_all.tolist()))

            # make coarse-grained s-by-by
            sad_coarse_all = []
            for focal in focal_all:
                focal_idx = numpy.where(taxon_n_all == focal)[0]
                sad_coarse_all.append(s_by_s_n[focal_idx,:].sum(axis=0))

            s_by_s_coarse = numpy.stack(sad_coarse_all, axis=0)
            idx_s_by_s_coarse = numpy.arange(s_by_s_coarse.shape[0])

            idx_taxon_n_all = numpy.arange(len(taxon_n_all))
            for focal_idx, focal in enumerate(focal_all):

                idx_focal = numpy.where(taxon_n_all==focal)[0]

                n_focal = len(idx_focal)

                if n_focal < 10:
                    continue

                s_by_s_n_focal = s_by_s_n[idx_focal,: ]
                s_by_s_n_non_focal = s_by_s_coarse[idx_s_by_s_coarse != focal_idx,: ]

                diversity_focal = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_n_focal)
                diversity_non_focal = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_n_non_focal)

                idx_diversity_to_keep = (diversity_focal>0) & (diversity_non_focal>0)
                diversity_focal = diversity_focal[idx_diversity_to_keep]
                diversity_non_focal = diversity_non_focal[idx_diversity_to_keep]

                if len(diversity_focal) < 5:
                    continue

                slope, intercept, r_value, p_value, std_err = stats.linregress(diversity_non_focal, diversity_focal)


                # get indices for null genera coarse graining
                # exclude focal genus!
                taxon_n_all_exclude_focal = taxon_n_all[taxon_n_all != focal]
                unique_non_focal, counts_non_focal = numpy.unique(taxon_n_all_exclude_focal, return_counts=True)
                coarse_grain_idx = numpy.append([0], numpy.cumsum(counts_non_focal))[:-1]

                #s_by_s_copy = numpy.copy(s_by_s)
                s_by_s_copy = numpy.copy(s_by_s_n)
                slope_null_all = []

                while len(slope_null_all) < iter:

                    numpy.random.shuffle(s_by_s_copy)
                    # each element is the number of ASVS belonging to the non-focal genus in a given site
                    s_by_s_focal_null = s_by_s_copy[:n_focal,:]
                    s_by_s_non_focal_null = s_by_s_copy[n_focal:,:]

                    s_by_s_non_focal_coarse_null = numpy.add.reduceat(s_by_s_non_focal_null, coarse_grain_idx, axis=0)

                    # remove sites where there are no observations in either focal or non-focal
                    idx_to_remove_null = numpy.all(s_by_s_focal_null == 0, axis=0) | numpy.all(s_by_s_non_focal_coarse_null == 0, axis=0)
                    s_by_s_focal_null = s_by_s_focal_null[:,~idx_to_remove_null]
                    s_by_s_non_focal_coarse_null = s_by_s_non_focal_coarse_null[:,~idx_to_remove_null]

                    diversity_focal_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_focal_null)
                    diversity_non_focal_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_non_focal_coarse_null)

                    idx_diversity_to_keep_null = (diversity_focal_null>0) & (diversity_non_focal_null>0)
                    diversity_focal_null = diversity_focal_null[idx_diversity_to_keep_null]
                    diversity_non_focal_null = diversity_non_focal_null[idx_diversity_to_keep_null]

                    # not enough
                    if len(diversity_non_focal_null) < len(diversity_non_focal):
                        continue

                    # randomly sample indexes
                    idx_subsample = numpy.random.choice(numpy.arange(len(diversity_focal_null)), size=len(diversity_non_focal), replace=False)

                    diversity_focal_null_subsample = diversity_focal_null[idx_subsample]
                    diversity_non_focal_null_subsample = diversity_non_focal_null[idx_subsample]

                    slope_null, intercept_null, r_value_null, p_value_null, std_err_null = stats.linregress(diversity_non_focal_null_subsample, diversity_focal_null_subsample)
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

                dbd_slope_dict[n]['ranks'][rank][focal] = {}
                dbd_slope_dict[n]['ranks'][rank][focal]['slope'] = slope
                dbd_slope_dict[n]['ranks'][rank][focal]['null_dist'] = slope_null_all.tolist()
                dbd_slope_dict[n]['ranks'][rank][focal]['p_value'] = p_value
                dbd_slope_dict[n]['ranks'][rank][focal]['percentile'] = percentile
                dbd_slope_dict[n]['ranks'][rank][focal]['lower_ci'] = lower_ci
                dbd_slope_dict[n]['ranks'][rank][focal]['upper_ci'] = upper_ci
                dbd_slope_dict[n]['ranks'][rank][focal]['n_asvs_focal_genus'] = n_focal
                dbd_slope_dict[n]['ranks'][rank][focal]['standard_score'] =  (slope - numpy.mean(slope_null_all)) / numpy.std(slope_null_all)



    sys.stderr.write("Saving slope dictionary...\n")
    with open(slope_resource_taxa_path, 'wb') as handle:
        pickle.dump(dbd_slope_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_diversity_slope_taxon_dalballo_dict():

    with open(slope_resource_taxa_path, 'rb') as handle:
        dbd_slope_dict = pickle.load(handle)
    return dbd_slope_dict



def make_diversity_slope_phylo_dalballo_dict(iter=1000):

    sample_dict = get_sample_dict_debello()
    taxon_dict = get_taxon_dict_dalballo()
    sra_to_sample_dict = load_sra_to_sample_dict()
    dbd_slope_taxon_dict = load_diversity_slope_taxon_dalballo_dict()

    n_carbons = list(sample_dict.keys())
    n_carbons.sort()

    sys.stderr.write("Loading site-by-species matrix...\n")
    s_by_s, species, sample_names = get_s_by_s_dalballo()
    species = numpy.asarray(species)

    sys.stderr.write("Subsetting samples...\n")
    samples_names_n_all = list(itertools.chain(*[dbd_slope_taxon_dict[n]['samples'] for n in n_carbons]))
    idx_samples_names_n_all = numpy.asarray([numpy.where(sample_names==s)[0][0] for s in samples_names_n_all])
    s_by_s_samples = s_by_s[:,idx_samples_names_n_all]

    species_samples = species[~(numpy.all(s_by_s_samples == 0, axis=1))].tolist()

    sys.stderr.write("Loading tree...\n")
    tree = tree_utils.get_dalbello_tree()
    tree_copy = tree.copy()

    sys.stderr.write("Collapsing tree by distance...\n")
    #taxonomy_names_union = list(set(itertools.chain(*taxonomy_names)))
    tree_subset = tree_utils.subset_tree(species_samples, tree_copy)

    distances = numpy.logspace(-3, 0.1, num=100)
    #distance_collapsed_dict = tree_utils.make_distance_collapsed_dict(tree_subset, distances)
    dbd_slope_dict = {}

    for n in n_carbons:

        sys.stderr.write("# resources = %d \n" % n)

        dbd_slope_dict[n] = {}

        sample_names_n = dbd_slope_taxon_dict[n]['samples']
        idx_samples_names_n = numpy.asarray([numpy.where(sample_names==s)[0][0] for s in sample_names_n])
        s_by_s_n = s_by_s[:,idx_samples_names_n]

        idx_non_zero_species = ~(numpy.all(s_by_s_n == 0, axis=1))
        s_by_s_n = s_by_s_n[idx_non_zero_species,:]
        species_n = species[idx_non_zero_species]

        tree_n = tree_utils.subset_tree(species_n.tolist(), tree_copy)
        distance_collapsed_n_dict = tree_utils.make_distance_collapsed_dict(tree_n, distances)

        #min_distance = min([d for d in distance_collapsed_n_dict.keys() if len(distance_collapsed_n_dict[d]) != len(species_n)])
        distances_n = list(distance_collapsed_n_dict.keys())
        distances_n.sort()

        for d in distances_n:

            dbd_slope_dict[n][d] = {}
            dbd_slope_dict[n][d]['slope'] = []
            dbd_slope_dict[n][d]['collapsed_nodes'] = []
            dbd_slope_dict[n][d]['null_dist'] = []
            dbd_slope_dict[n][d]['p_value'] = []
            dbd_slope_dict[n][d]['percentile'] = []
            dbd_slope_dict[n][d]['lower_ci'] = []
            dbd_slope_dict[n][d]['upper_ci'] = []
            dbd_slope_dict[n][d]['n_asvs_focal_genus'] = []
            dbd_slope_dict[n][d]['standard_score'] = []
            dbd_slope_dict[n][d]['asv'] = []

            coarse_grained_list = distance_collapsed_n_dict[d]
            coarse_grained_n = [len(i) for i in coarse_grained_list]

            # get indexes for each clade
            coarse_grained_idx_all = [numpy.asarray([numpy.where(species_n==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            # coarse grain s-by-s for all clades
            s_by_s_n_all_clades = numpy.stack([numpy.sum(s_by_s_n[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
            idx_all = numpy.arange(len(coarse_grained_idx_all))
            for i in idx_all:

                if len(coarse_grained_idx_all[i]) < 10:
                    continue

                s_by_s_n_focal = s_by_s_n[coarse_grained_idx_all[i],:]
                s_by_s_n_non_focal = s_by_s_n_all_clades[idx_all!=i,:]

                diversity_focal = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_n_focal)
                diversity_non_focal = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_n_non_focal)

                idx_to_keep = (diversity_focal>0) & (diversity_non_focal>0) & (numpy.isfinite(diversity_focal)) & (numpy.isfinite(diversity_non_focal))
                diversity_focal = diversity_focal[idx_to_keep]
                diversity_non_focal = diversity_non_focal[idx_to_keep]

                # ignore cut-offs with low richness
                if len(diversity_focal) < 5:
                    continue

                # ignore cut-offs where all diversity estimates are the same
                if (sum(diversity_non_focal==1) == len(diversity_non_focal)) or (sum(diversity_focal==1) == len(diversity_focal)):
                    continue

                slope, intercept, r_value, p_value, std_err = stats.linregress(diversity_non_focal, diversity_focal)
                s_by_s_n_copy = numpy.copy(s_by_s_n)

                slope_null_all = []
                while len(slope_null_all) < iter:
                    numpy.random.shuffle(s_by_s_n_copy)

                    s_by_s_n_focal_null = s_by_s_n_copy[:coarse_grained_n[i],:]
                    s_by_s_n_non_focal_null = s_by_s_n_copy[coarse_grained_n[i]:,:]

                    diversity_focal_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_n_focal_null)
                    diversity_non_focal_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_n_non_focal_null)

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

                    #if (len(slope_null_all) % 100 == 0):
                    #    sys.stderr.write("Iteration %d...\n" % len(slope_null_all))


                slope_null_all = numpy.asarray(slope_null_all)
                slope_null_all = numpy.sort(slope_null_all)

                slope_null_median = numpy.median(slope_null_all)
                if slope >= slope_null_median:
                    p_value = sum(slope_null_all > slope)/iter
                else:
                    p_value = sum(slope_null_all < slope)/iter

                percentile = sum(slope_null_all < slope)/iter

                lower_ci = slope_null_all[int(iter*0.025)]
                upper_ci = slope_null_all[int(iter*0.975)]

                n_asvs_focal_genus = s_by_s_n_focal.shape[0]

                dbd_slope_dict[n][d]['slope'].append(slope)
                dbd_slope_dict[n][d]['collapsed_nodes'].append(coarse_grained_list[i])
                dbd_slope_dict[n][d]['null_dist'].append(slope_null_all.tolist())
                dbd_slope_dict[n][d]['p_value'].append(p_value)
                dbd_slope_dict[n][d]['percentile'].append(percentile)
                dbd_slope_dict[n][d]['lower_ci'].append(lower_ci)
                dbd_slope_dict[n][d]['upper_ci'].append(upper_ci)
                dbd_slope_dict[n][d]['n_asvs_focal_genus'].append(n_asvs_focal_genus)
                dbd_slope_dict[n][d]['standard_score'].append((slope - numpy.mean(slope_null_all)) / numpy.std(slope_null_all))
                dbd_slope_dict[n][d]['asv'].append(coarse_grained_list[i])


    sys.stderr.write("Saving slope dictionary...\n")
    with open(slope_resource_phylo_path, 'wb') as handle:
        pickle.dump(dbd_slope_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_diversity_slope_phylo_dalballo_dict():

    with open(slope_resource_phylo_path, 'rb') as handle:
        dbd_slope_dict = pickle.load(handle)
    return dbd_slope_dict





def make_diversity_vs_diversity_taxon_dict(rarefied = True, iter = 1000):

    diversity_vs_diversity_taxon_dict = {}
    for environment in diversity_utils.environments_to_keep:

        sys.stderr.write("Loading taxon dict for %s...\n" % environment)
        sad_annotated_dict = load_sad_annotated_taxon_dict(environment, rarefied=rarefied)
        #diversity_slope_dict = dbd_utils.load_diversity_slope_taxon_dict(environment)

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        rel_s_by_s = s_by_s/numpy.sum(s_by_s, axis=0)
        mean_asv = numpy.mean(rel_s_by_s, axis=1)
        mean_asv_copy = numpy.copy(mean_asv)
        y_rel_s_by_s = (rel_s_by_s.T/mean_asv).T
        y_rel_s_by_s_copy = numpy.copy(y_rel_s_by_s)

        diversity_species = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s)

        diversity_vs_diversity_taxon_dict[environment] = {}
        diversity_vs_diversity_taxon_dict[environment]['diversity_species'] = diversity_species.tolist()

        for rank in diversity_utils.taxa_ranks:

            sys.stderr.write("Starting %s level analysis...\n" % rank)

            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa])))
            all_genera_idx = numpy.arange(len(all_genera))

            genus_to_taxa_dict = {}
            g_taxa_idx_all = []
            sad_genera_all = []

            for genus in all_genera:
                genus_to_taxa_dict[genus] = []
                for t in taxa:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
                g_taxa_idx_all.append(g_taxa_idx)
                sad_genera_all.append(s_by_s[g_taxa_idx,:].sum(axis=0))

            s_by_s_genera = numpy.stack(sad_genera_all, axis=0)
            # remove sites where there are no observations
            s_by_s_genera = s_by_s_genera[:,~(numpy.all(s_by_s_genera == 0, axis=0))]
            diversity_genera = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_genera)

            slope, intercept, r_value, p_value, std_err = stats.linregress(diversity_species, diversity_genera)

            # get indices for null genera coarse graining
            unique_genera, counts_genera = numpy.unique([sad_annotated_dict['taxa'][t][rank] for t in taxa], return_counts=True)
            coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_genera))[:-1]

            s_by_s_copy = numpy.copy(s_by_s)
            slope_null_all = []
            mean_diversity_coarse_null_all = []
            mean_diversity_coarse_constrain_mad_all = []
            mean_diversity_coarse_constrain_y_all = []
            for i in range(iter):

                numpy.random.shuffle(s_by_s_copy)
                s_by_s_coarse_null = numpy.add.reduceat(s_by_s_copy, coarse_grain_by_genus_idx, axis=0)
                s_by_s_coarse_null = s_by_s_coarse_null[:,~(numpy.all(s_by_s_coarse_null == 0, axis=0))]
                diversity_genera_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_coarse_null)
                slope_null, intercept_null, r_value_null, p_value_null, std_err_null = stats.linregress(diversity_species, diversity_genera_null)
                slope_null_all.append(slope_null)
                mean_diversity_coarse_null_all.append(numpy.mean(diversity_genera_null))


                # constrain on mean, shuffle y
                numpy.random.shuffle(y_rel_s_by_s_copy)
                # multiply by MAD
                rel_s_by_s_constrain_mad = (y_rel_s_by_s_copy.T*mean_asv).T
                # renormalize
                rel_s_by_s_constrain_mad = rel_s_by_s_constrain_mad/numpy.sum(rel_s_by_s_constrain_mad, axis=0)
                # add
                #rel_s_by_s_constrain_mad_coarse = numpy.add.reduceat(rel_s_by_s_constrain_mad, coarse_grain_by_genus_idx, axis=0)
                rel_s_by_s_constrain_mad_coarse = [rel_s_by_s_constrain_mad[g_taxa_idx_i,:].sum(axis=0) for g_taxa_idx_i in g_taxa_idx_all]
                rel_s_by_s_constrain_mad_coarse = numpy.asarray(rel_s_by_s_constrain_mad_coarse)
                diversity_genera_constrain_mad = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_s_by_s_constrain_mad_coarse)
                mean_diversity_coarse_constrain_mad_all.append(numpy.mean(diversity_genera_constrain_mad))


                # constrain on y, shuffle mean
                numpy.random.shuffle(mean_asv_copy)
                rel_s_by_s_constrain_y = (y_rel_s_by_s.T*mean_asv_copy).T
                rel_s_by_s_constrain_y = rel_s_by_s_constrain_y/numpy.sum(rel_s_by_s_constrain_y, axis=0)
                # coarse-grain
                rel_s_by_s_constrain_y_coarse = [rel_s_by_s_constrain_y[g_taxa_idx_i,:].sum(axis=0) for g_taxa_idx_i in g_taxa_idx_all]
                rel_s_by_s_constrain_y_coarse = numpy.asarray(rel_s_by_s_constrain_y_coarse)

                diversity_genera_constrain_y = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_s_by_s_constrain_y_coarse)
                mean_diversity_coarse_constrain_y_all.append(numpy.mean(diversity_genera_constrain_y))



            slope_null_all = numpy.asarray(slope_null_all)
            slope_null_all = numpy.sort(slope_null_all)

            percentile = sum(slope_null_all < slope)/iter

            lower_ci = slope_null_all[int(iter*0.025)]
            upper_ci = slope_null_all[int(iter*0.975)]


            mean_diversity_coarse_null_all = numpy.asarray(mean_diversity_coarse_null_all)
            mean_diversity_coarse_null_all = numpy.sort(mean_diversity_coarse_null_all)

            mean_diversity_coarse_null_lower_ci = mean_diversity_coarse_null_all[int(iter*0.025)]
            mean_diversity_coarse_null_upper_ci = mean_diversity_coarse_null_all[int(iter*0.975)]

            # constrain mad
            mean_diversity_coarse_constrain_mad_all = numpy.asarray(mean_diversity_coarse_constrain_mad_all)
            mean_diversity_coarse_constrain_mad_all = numpy.sort(mean_diversity_coarse_constrain_mad_all)

            mean_diversity_coarse_constrain_mad_lower_ci = mean_diversity_coarse_constrain_mad_all[int(iter*0.025)]
            mean_diversity_coarse_constrain_mad_upper_ci = mean_diversity_coarse_constrain_mad_all[int(iter*0.975)]


            # constrain y
            mean_diversity_coarse_constrain_y_all = numpy.asarray(mean_diversity_coarse_constrain_y_all)
            mean_diversity_coarse_constrain_y_all = numpy.sort(mean_diversity_coarse_constrain_y_all)

            mean_diversity_coarse_constrain_y_lower_ci = mean_diversity_coarse_constrain_y_all[int(iter*0.025)]
            mean_diversity_coarse_constrain_y_upper_ci = mean_diversity_coarse_constrain_y_all[int(iter*0.975)]



            diversity_vs_diversity_taxon_dict[environment][rank] = {}
            diversity_vs_diversity_taxon_dict[environment][rank]['diversity'] = diversity_genera.tolist()
            diversity_vs_diversity_taxon_dict[environment][rank]['slope'] = slope
            diversity_vs_diversity_taxon_dict[environment][rank]['intercept'] = intercept
            diversity_vs_diversity_taxon_dict[environment][rank]['slope_percentile'] = percentile
            diversity_vs_diversity_taxon_dict[environment][rank]['slope_null_lower_ci'] = lower_ci
            diversity_vs_diversity_taxon_dict[environment][rank]['slope_null_upper_ci'] = upper_ci

            diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_observed'] = numpy.mean(diversity_genera)

            diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_null'] = numpy.mean(mean_diversity_coarse_null_all)
            diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_null_lower_ci'] = mean_diversity_coarse_null_lower_ci
            diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_null_upper_ci'] = mean_diversity_coarse_null_upper_ci

            diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_null_constrain_mad'] = numpy.mean(mean_diversity_coarse_constrain_mad_all)
            diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_null_constrain_mad_lower_ci'] = mean_diversity_coarse_constrain_mad_lower_ci
            diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_null_constrain_mad_upper_ci'] = mean_diversity_coarse_constrain_mad_upper_ci


            diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_null_constrain_y'] = numpy.mean(mean_diversity_coarse_constrain_y_all)
            diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_null_constrain_y_lower_ci'] = mean_diversity_coarse_constrain_y_lower_ci
            diversity_vs_diversity_taxon_dict[environment][rank]['mean_diversity_null_constrain_y_upper_ci'] = mean_diversity_coarse_constrain_y_upper_ci


    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    diversity_vs_diversity_taxon_path = diversity_vs_diversity_taxon_template % rarefied_label

    sys.stderr.write("Saving diversity vs. diversity dictionary...\n")
    with open(diversity_vs_diversity_taxon_path, 'wb') as handle:
        pickle.dump(diversity_vs_diversity_taxon_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_diversity_vs_diversity_taxon_dict(rarefied=True):

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    diversity_vs_diversity_taxon_path = diversity_vs_diversity_taxon_template % rarefied_label

    with open(diversity_vs_diversity_taxon_path, 'rb') as handle:
        diversity_vs_diversity_taxon_dict = pickle.load(handle)
    return diversity_vs_diversity_taxon_dict





def make_diversity_vs_diversity_phylo_dict(rarefied = False, iter = 100):

    diversity_vs_diversity_phylo_dict = {}
    for environment in diversity_utils.environments_to_keep:

        sys.stderr.write("Loading tree dict for %s...\n" % environment)
        coarse_grained_tree_dict = load_coarse_grained_tree_dict(environment=environment, rarefied=rarefied)

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        pres_abs_dict = load_sad_annotated_phylo_dict(environment, rarefied = rarefied)

        samples = pres_abs_dict['samples']
        taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))

        sys.stderr.write("Getting site-by-species matrix...\n")
        s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

        rel_s_by_s = s_by_s/numpy.sum(s_by_s, axis=0)
        mean_asv = numpy.mean(rel_s_by_s, axis=1)
        #mean_asv_copy = numpy.copy(mean_asv)
        y_rel_s_by_s = (rel_s_by_s.T/mean_asv).T
        #y_rel_s_by_s_copy = numpy.copy(y_rel_s_by_s)

        distances = list(coarse_grained_tree_dict.keys())
        distances.sort()

        diversity_species = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s)

        diversity_vs_diversity_phylo_dict[environment] = {}
        diversity_vs_diversity_phylo_dict[environment]['diversity_species'] = diversity_species.tolist()

        for distance in distances:

            sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
            coarse_grained_list = coarse_grained_tree_dict[distance]
            coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

            #if len(coarse_grained_n) == len(taxonomy_names):
            #    continue

            # get indexes for each clade
            coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            # re-sort the data
            coarse_grained_idx_all_flat = []
            for c in coarse_grained_idx_all:
                coarse_grained_idx_all_flat.extend(c.tolist())

            coarse_grained_idx_all_flat = numpy.asarray(coarse_grained_idx_all_flat)
            n_coarse_grained_all = numpy.asarray([len(c) for c in coarse_grained_idx_all])
            coarse_grain_idx = numpy.append([0], numpy.cumsum(n_coarse_grained_all))[:-1]

            s_by_s_distance = s_by_s[coarse_grained_idx_all_flat,:]
            rel_s_by_s_distance = rel_s_by_s[coarse_grained_idx_all_flat,:]
            y_rel_s_by_s_distance = y_rel_s_by_s[coarse_grained_idx_all_flat,:]
            mean_asv_distance = mean_asv[coarse_grained_idx_all_flat]

            ########################
            # constrain on occupancy
            ########################
            # figure out number of blocks
            #block_richnes = 100
            #n_blocks = math.ceil(len(mean_asv_distance) / block_richnes)
            #block_richness_all = [block_richnes]*(n_blocks-1)
            # add remaining species
            #block_richness_all.append(len(mean_asv_distance) - block_richnes*(n_blocks-1))
            #block_richness_all = numpy.asarray(block_richness_all)

            # create array of block numbers
            # sort rel_s_by_s by mean abundance
            #s_idx = list(range(s_by_s.shape[0]))
            #mean_abundance_s_idx_tuple_all = list(zip(mean_abundance.tolist(), s_idx))
            # sort tuples by mean abundance, order of increasing abundane
            #mean_abundance_s_idx_tuple_all.sort(key=lambda tup: tup[0])
            #s_sort_idx = numpy.asarray([s[1] for s in mean_abundance_s_idx_tuple_all])
            #mean_abundance_sort = numpy.asarray([s[0] for s in mean_abundance_s_idx_tuple_all])

            #s_by_s_distance = s_by_s_distance[s_sort_idx,:]
            #rel_s_by_s_distance = rel_s_by_s_distance[s_sort_idx,:]
            #y_rel_s_by_s_distance = y_rel_s_by_s_distance[s_sort_idx,:]
            #mean_asv_distance = mean_asv_distance[s_sort_idx]




            mean_asv_distance_copy = numpy.copy(mean_asv_distance)
            y_rel_s_by_s_distance_copy = numpy.copy(y_rel_s_by_s_distance)
            s_by_s_distance_copy = numpy.copy(s_by_s_distance)

            s_by_s_distance_coarse = numpy.add.reduceat(s_by_s_distance_copy, coarse_grain_idx, axis=0)
            diversity_coarse = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_distance_coarse)


            # coarse grain s-by-s for all clades
            #s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
            slope, intercept, r_value, p_value, std_err = stats.linregress(diversity_species, diversity_coarse)

            idx_all = numpy.arange(len(coarse_grained_idx_all))
            #coarse_grain_idx = numpy.append([0], numpy.cumsum(coarse_grained_n))[:-1]
            #s_by_s_copy = numpy.copy(s_by_s)
            slope_null_all = []
            mean_diversity_coarse_null_all = []
            mean_diversity_coarse_constrain_mad_all = []
            mean_diversity_coarse_constrain_y_all = []
            for i in range(iter):

                numpy.random.shuffle(s_by_s_distance_copy)
                s_by_s_coarse_null = numpy.add.reduceat(s_by_s_distance_copy, coarse_grain_idx, axis=0)
                s_by_s_coarse_null = s_by_s_coarse_null[:,~(numpy.all(s_by_s_coarse_null == 0, axis=0))]
                diversity_coarse_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_coarse_null)
                mean_diversity_coarse_null_all.append(numpy.mean(diversity_coarse_null))

                slope_null, intercept_null, r_value_null, p_value_null, std_err_null = stats.linregress(diversity_species, diversity_coarse_null)
                slope_null_all.append(slope_null)


                # constrain on mean, shuffle y
                numpy.random.shuffle(y_rel_s_by_s_distance_copy)
                # multiply by MAD
                rel_s_by_s_constrain_mad = (y_rel_s_by_s_distance_copy.T*mean_asv_distance).T
                # renormalize
                rel_s_by_s_constrain_mad = rel_s_by_s_constrain_mad/numpy.sum(rel_s_by_s_constrain_mad, axis=0)
                # add
                #rel_s_by_s_constrain_mad_coarse = [rel_s_by_s_constrain_mad[g_taxa_idx_i,:].sum(axis=0) for g_taxa_idx_i in coarse_grained_idx_all]
                #rel_s_by_s_constrain_mad_coarse = numpy.asarray(rel_s_by_s_constrain_mad_coarse)
                rel_s_by_s_constrain_mad_coarse = numpy.add.reduceat(rel_s_by_s_constrain_mad, coarse_grain_idx, axis=0)
                diversity_genera_constrain_mad = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_s_by_s_constrain_mad_coarse)
                mean_diversity_coarse_constrain_mad_all.append(numpy.mean(diversity_genera_constrain_mad))


                # constrain on mean and occupancy, shuffle y



                # constrain on y, shuffle mean
                numpy.random.shuffle(mean_asv_distance_copy)
                rel_s_by_s_constrain_y = (y_rel_s_by_s_distance.T*mean_asv_distance_copy).T
                rel_s_by_s_constrain_y = rel_s_by_s_constrain_y/numpy.sum(rel_s_by_s_constrain_y, axis=0)
                # coarse-grain
                #rel_s_by_s_constrain_y_coarse = [rel_s_by_s_constrain_y[g_taxa_idx_i,:].sum(axis=0) for g_taxa_idx_i in coarse_grained_idx_all]
                #rel_s_by_s_constrain_y_coarse = numpy.asarray(rel_s_by_s_constrain_y_coarse)
                rel_s_by_s_constrain_y_coarse = numpy.add.reduceat(rel_s_by_s_constrain_y, coarse_grain_idx, axis=0)

                diversity_genera_constrain_y = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_s_by_s_constrain_y_coarse)
                mean_diversity_coarse_constrain_y_all.append(numpy.mean(diversity_genera_constrain_y))



            slope_null_all = numpy.asarray(slope_null_all)
            slope_null_all = numpy.sort(slope_null_all)

            percentile = sum(slope_null_all < slope)/iter

            lower_ci = slope_null_all[int(iter*0.025)]
            upper_ci = slope_null_all[int(iter*0.975)]

            mean_diversity_coarse_null_all = numpy.asarray(mean_diversity_coarse_null_all)
            mean_diversity_coarse_null_all = numpy.sort(mean_diversity_coarse_null_all)

            mean_diversity_coarse_null_lower_ci = mean_diversity_coarse_null_all[int(iter*0.025)]
            mean_diversity_coarse_null_upper_ci = mean_diversity_coarse_null_all[int(iter*0.975)]


            # constrain mad
            mean_diversity_coarse_constrain_mad_all = numpy.asarray(mean_diversity_coarse_constrain_mad_all)
            mean_diversity_coarse_constrain_mad_all = numpy.sort(mean_diversity_coarse_constrain_mad_all)

            mean_diversity_coarse_constrain_mad_lower_ci = mean_diversity_coarse_constrain_mad_all[int(iter*0.025)]
            mean_diversity_coarse_constrain_mad_upper_ci = mean_diversity_coarse_constrain_mad_all[int(iter*0.975)]


            # constrain y
            mean_diversity_coarse_constrain_y_all = numpy.asarray(mean_diversity_coarse_constrain_y_all)
            mean_diversity_coarse_constrain_y_all = numpy.sort(mean_diversity_coarse_constrain_y_all)

            mean_diversity_coarse_constrain_y_lower_ci = mean_diversity_coarse_constrain_y_all[int(iter*0.025)]
            mean_diversity_coarse_constrain_y_upper_ci = mean_diversity_coarse_constrain_y_all[int(iter*0.975)]


            diversity_vs_diversity_phylo_dict[environment][distance] = {}
            diversity_vs_diversity_phylo_dict[environment][distance]['diversity'] = diversity_coarse.tolist()
            diversity_vs_diversity_phylo_dict[environment][distance]['slope'] = slope
            diversity_vs_diversity_phylo_dict[environment][distance]['slope_percentile'] = percentile
            diversity_vs_diversity_phylo_dict[environment][distance]['slope_null_lower_ci'] = lower_ci
            diversity_vs_diversity_phylo_dict[environment][distance]['slope_null_upper_ci'] = upper_ci


            diversity_vs_diversity_phylo_dict[environment][distance]['mean_diversity_observed'] = numpy.mean(diversity_coarse)


            diversity_vs_diversity_phylo_dict[environment][distance]['mean_diversity_null'] = numpy.mean(mean_diversity_coarse_null_all)
            diversity_vs_diversity_phylo_dict[environment][distance]['mean_diversity_null_lower_ci'] = mean_diversity_coarse_null_lower_ci
            diversity_vs_diversity_phylo_dict[environment][distance]['mean_diversity_null_upper_ci'] = mean_diversity_coarse_null_upper_ci


            diversity_vs_diversity_phylo_dict[environment][distance]['mean_diversity_null_constrain_mad'] = numpy.mean(mean_diversity_coarse_constrain_mad_all)
            diversity_vs_diversity_phylo_dict[environment][distance]['mean_diversity_null_constrain_mad_lower_ci'] = mean_diversity_coarse_constrain_mad_lower_ci
            diversity_vs_diversity_phylo_dict[environment][distance]['mean_diversity_null_constrain_mad_upper_ci'] = mean_diversity_coarse_constrain_mad_upper_ci

            diversity_vs_diversity_phylo_dict[environment][distance]['mean_diversity_null_constrain_y'] = numpy.mean(mean_diversity_coarse_constrain_y_all)
            diversity_vs_diversity_phylo_dict[environment][distance]['mean_diversity_null_constrain_y_lower_ci'] = mean_diversity_coarse_constrain_y_lower_ci
            diversity_vs_diversity_phylo_dict[environment][distance]['mean_diversity_null_constrain_y_upper_ci'] = mean_diversity_coarse_constrain_y_upper_ci



    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    diversity_vs_diversity_phylo_path = diversity_vs_diversity_phylo_emp_template % rarefied_label

    sys.stderr.write("Saving diversity vs. diversity dictionary...\n")
    with open(diversity_vs_diversity_phylo_path, 'wb') as handle:
        pickle.dump(diversity_vs_diversity_phylo_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)





def load_diversity_vs_diversity_phylo_dict(rarefied=True):

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    diversity_vs_diversity_phylo_path = diversity_vs_diversity_phylo_emp_template % rarefied_label

    with open(diversity_vs_diversity_phylo_path, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_






def make_diversity_vs_diversity_taxon_dalbello_dict(iter = 1000):

    sample_dict = get_sample_dict_debello()
    taxon_dict = get_taxon_dict_dalballo()
    sra_to_sample_dict = load_sra_to_sample_dict()

    n_carbons = list(sample_dict.keys())
    n_carbons.sort()

    s_by_s, species, sample_names = get_s_by_s_dalballo()
    species = numpy.asarray(species)

    diversity_vs_diversity_taxon_dict = {}

    for n in n_carbons:

        # use this to fix order of resource combinations
        resource_combo_all = list(sample_dict[n].keys())

        diversity_vs_diversity_taxon_dict[n] = {}

        # get all samples with n resources

        samples_n = numpy.asarray(numpy.concatenate([sample_dict[n][resource_combo_] for resource_combo_ in resource_combo_all]).flat)
        idx_n = numpy.asarray([numpy.where(sample_names==s)[0][0] for s in samples_n])
        # sort sample names by order in s_by_s
        samples_n = sample_names[idx_n]
        s_by_s_n = s_by_s[:,idx_n]

        for rank in diversity_utils.taxa_ranks:

            # make copy for each rank
            s_by_s_n_copy = numpy.copy(s_by_s_n)

            taxon_all = numpy.asarray([taxon_dict[s][rank] for s in species])

            # remove species with no annotation
            idx_to_keep = (taxon_all!='NA')
            s_by_s_n_copy = s_by_s_n_copy[idx_to_keep,:]
            species_n = species[idx_to_keep]

            # remove sites with all zeros
            s_by_s_n_copy = s_by_s_n_copy[:,~(numpy.all(s_by_s_n_copy == 0, axis=0))]
            # remove species with all zeros
            idx_species_to_keep = ~(numpy.all(s_by_s_n_copy == 0, axis=1))
            s_by_s_n_copy = s_by_s_n_copy[idx_species_to_keep,:]

            species_n = species_n[idx_species_to_keep]
            taxon_n_all = numpy.asarray([taxon_dict[s][rank] for s in species_n])
            #focal_all = list(set(taxon_n_all.tolist()))

            unique_taxon_n, counts_taxon_n = numpy.unique(taxon_n_all, return_counts=True)
            coarse_grain_idx = numpy.append([0], numpy.cumsum(counts_taxon_n))[:-1]
            s_by_s_n_coarse = numpy.add.reduceat(s_by_s_n_copy, coarse_grain_idx, axis=0)

            # make coarse-grained s-by-by
            #coarse_grained_n = []
            #sad_coarse_all = []
            #for focal in focal_all:
            #    focal_idx = numpy.where(taxon_n_all == focal)[0]
            #    coarse_grained_n.append(len(focal_idx))
            #    sad_coarse_all.append(s_by_s_n_copy[focal_idx,:].sum(axis=0))

            #coarse_grained_n = numpy.asarray(coarse_grained_n)
            #s_by_s_n_coarse = numpy.stack(sad_coarse_all, axis=0)
            #idx_s_by_s_coarse = numpy.arange(s_by_s_coarse.shape[0])
            #idx_taxon_n_all = numpy.arange(len(taxon_n_all))


            # now generate null coarse-grained counts within each resource_combo block
            diversity_species = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_n_copy)
            diversity_coarse = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_n_coarse)
            slope, intercept, r_value, p_value, std_err = stats.linregress(diversity_species, diversity_coarse)

            # coarse_grained idx for each resource combination in each n carbon
            coarse_grain_n_resource_combo_idx_dict = {}
            s_by_s_n_resource_combo_dict = {}
            for resource_combo in resource_combo_all:

                resource_combo_samples = sample_dict[n][resource_combo]
                idx_n_resource_combo = numpy.asarray([numpy.where(samples_n==s)[0][0] for s in resource_combo_samples])
                s_by_s_n_resource_combo = s_by_s_n_copy[:,idx_n_resource_combo]
                # remove species with all zeros
                idx_n_resource_combo_species_to_keep = ~(numpy.all(s_by_s_n_resource_combo == 0, axis=1))
                s_by_s_n_resource_combo = s_by_s_n_resource_combo[idx_n_resource_combo_species_to_keep,:]
                species_n_resource_combo = species_n[idx_n_resource_combo_species_to_keep]
                taxon_n_resource_combo_all = numpy.asarray([taxon_dict[s][rank] for s in species_n_resource_combo])
                unique_taxon_n_resource_combo, counts_taxon_n_resource_combo = numpy.unique(taxon_n_resource_combo_all, return_counts=True)
                coarse_grain_by_taxon_n_resource_combo_idx = numpy.append([0], numpy.cumsum(counts_taxon_n_resource_combo))[:-1]

                s_by_s_n_resource_combo_dict[resource_combo] = s_by_s_n_resource_combo
                coarse_grain_n_resource_combo_idx_dict[resource_combo] = coarse_grain_by_taxon_n_resource_combo_idx


            slope_null_all = []
            for i in range(iter):

                #diversity_coarse_null_all = []

                #for resource_combo in resource_combo_all:

                #    s_by_s_n_resource_combo = s_by_s_n_resource_combo_dict[resource_combo]
                #    coarse_grain_n_resource_combo_idx = coarse_grain_n_resource_combo_idx_dict[resource_combo]
                #    numpy.random.shuffle(s_by_s_n_resource_combo)
                #    # coarse-grain null
                #    s_by_s_n_resource_combo_coarse_null = numpy.add.reduceat(s_by_s_n_resource_combo, coarse_grain_n_resource_combo_idx, axis=0)
                #    diversity_coarse_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_n_resource_combo_coarse_null)
                #    diversity_coarse_null_all.extend(diversity_coarse_null)

                numpy.random.shuffle(s_by_s_n_copy)
                s_by_s_n_coarse_null = numpy.add.reduceat(s_by_s_n_copy, coarse_grain_idx, axis=0)
                diversity_coarse_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_n_coarse_null)

                slope_null, intercept_null, r_value_null, p_value_null, std_err_null = stats.linregress(diversity_species, diversity_coarse_null)
                slope_null_all.append(slope_null)


            slope_null_all = numpy.asarray(slope_null_all)
            slope_null_all = numpy.sort(slope_null_all)

            percentile = sum(slope_null_all < slope)/iter

            lower_ci = slope_null_all[int(iter*0.025)]
            upper_ci = slope_null_all[int(iter*0.975)]

            diversity_vs_diversity_taxon_dict[n][rank] = {}
            diversity_vs_diversity_taxon_dict[n][rank]['diversity_all'] = diversity_species.tolist()
            diversity_vs_diversity_taxon_dict[n][rank]['diversity_coarse'] = diversity_coarse.tolist()
            diversity_vs_diversity_taxon_dict[n][rank]['slope'] = slope
            diversity_vs_diversity_taxon_dict[n][rank]['percentile'] = percentile
            diversity_vs_diversity_taxon_dict[n][rank]['lower_ci'] = lower_ci
            diversity_vs_diversity_taxon_dict[n][rank]['upper_ci'] = upper_ci


    diversity_vs_diversity_taxon_path = diversity_vs_diversity_taxon_template % ('_resource')

    sys.stderr.write("Saving diversity vs. diversity dictionary...\n")
    with open(diversity_vs_diversity_taxon_path, 'wb') as handle:
        pickle.dump(diversity_vs_diversity_taxon_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)



def load_diversity_vs_diversity_taxon_resource_dict():

    diversity_vs_diversity_taxon_path = diversity_vs_diversity_taxon_template % ('_resource')

    with open(diversity_vs_diversity_taxon_path, 'rb') as handle:
        diversity_vs_diversity_taxon_dict = pickle.load(handle)
    return diversity_vs_diversity_taxon_dict






def make_diversity_vs_diversity_phylo_dalbello_dict(iter = 1000):

    sample_dict = get_sample_dict_debello()
    taxon_dict = get_taxon_dict_dalballo()
    sra_to_sample_dict = load_sra_to_sample_dict()

    n_carbons = list(sample_dict.keys())
    n_carbons.sort()

    s_by_s, species, sample_names = get_s_by_s_dalballo()
    species = numpy.asarray(species)

    sys.stderr.write("Loading tree...\n")
    tree = tree_utils.get_dalbello_tree()
    tree_copy = tree.copy()

    distances = numpy.logspace(-3, 1, num=40)

    diversity_vs_diversity_phylo_dict = {}
    for n in n_carbons:

        diversity_vs_diversity_phylo_dict[n] = {}

        resource_combo_all = list(sample_dict[n].keys())
        # get all samples with n resources
        samples_n = numpy.asarray(numpy.concatenate([sample_dict[n][resource_combo_] for resource_combo_ in resource_combo_all]).flat)
        idx_n = numpy.asarray([numpy.where(sample_names==s)[0][0] for s in samples_n])
        s_by_s_n = s_by_s[:,idx_n]

        idx_non_zero_species = ~(numpy.all(s_by_s_n == 0, axis=1))
        s_by_s_n = s_by_s_n[idx_non_zero_species,:]
        species_n = species[idx_non_zero_species]

        tree_n = tree_utils.subset_tree(species_n.tolist(), tree_copy)
        # only keeps a distance if there at least ten taxa
        distance_collapsed_n_dict = tree_utils.make_distance_collapsed_dict(tree_n, distances)
        distances_n = list(distance_collapsed_n_dict.keys())
        distances_n.sort()


        diversity_species = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_n)

        for d in distances_n:

            coarse_grained_list = distance_collapsed_n_dict[d]
            coarse_grained_n = [len(i) for i in coarse_grained_list]

            # get indexes for each clade
            coarse_grained_idx_all = [numpy.asarray([numpy.where(species_n==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            # coarse grain s-by-s for all clades
            s_by_s_n_coarse = numpy.stack([numpy.sum(s_by_s_n[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
            idx_all = numpy.arange(len(coarse_grained_idx_all))

            diversity_coarse = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_n_coarse)
            slope, intercept, r_value, p_value, std_err = stats.linregress(diversity_species, diversity_coarse)


            coarse_grained_n = numpy.asarray(coarse_grained_n)
            coarse_grain_idx = numpy.append([0], numpy.cumsum(coarse_grained_n))[:-1]

            # coarse_grained idx for each resource combination in each n carbon
            #coarse_grain_n_resource_combo_idx_dict = {}
            #s_by_s_n_resource_combo_dict = {}
            #for resource_combo in resource_combo_all:

            #    resource_combo_samples = sample_dict[n][resource_combo]
            #    idx_n_resource_combo = numpy.asarray([numpy.where(samples_n==s)[0][0] for s in resource_combo_samples])
            #    s_by_s_n_resource_combo = s_by_s_n[:,idx_n_resource_combo]
            #    # remove species with all zeros
            #    idx_n_resource_combo_species_to_keep = ~(numpy.all(s_by_s_n_resource_combo == 0, axis=1))
            #    s_by_s_n_resource_combo = s_by_s_n_resource_combo[idx_n_resource_combo_species_to_keep,:]
            #    species_n_resource_combo = species_n[idx_n_resource_combo_species_to_keep]

            #    taxon_n_resource_combo_all = numpy.asarray([taxon_dict[s][rank] for s in species_n_resource_combo])
            #    unique_taxon_n_resource_combo, counts_taxon_n_resource_combo = numpy.unique(taxon_n_resource_combo_all, return_counts=True)
            #    coarse_grain_by_taxon_n_resource_combo_idx = numpy.append([0], numpy.cumsum(counts_taxon_n_resource_combo))[:-1]

            #    s_by_s_n_resource_combo_dict[resource_combo] = s_by_s_n_resource_combo
            #    coarse_grain_n_resource_combo_idx_dict[resource_combo] = coarse_grain_by_taxon_n_resource_combo_idx


            #
            slope_null_all = []
            s_by_s_n_copy = numpy.copy(s_by_s_n)
            for i in range(iter):

                numpy.random.shuffle(s_by_s_n_copy)
                s_by_s_n_coarse_null = numpy.add.reduceat(s_by_s_n_copy, coarse_grain_idx, axis=0)
                diversity_coarse_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_n_coarse_null)

                slope_null, intercept_null, r_value_null, p_value_null, std_err_null = stats.linregress(diversity_species, diversity_coarse_null)
                slope_null_all.append(slope_null)

            slope_null_all = numpy.asarray(slope_null_all)
            slope_null_all = numpy.sort(slope_null_all)

            percentile = sum(slope_null_all < slope)/iter

            lower_ci = slope_null_all[int(iter*0.025)]
            upper_ci = slope_null_all[int(iter*0.975)]

            diversity_vs_diversity_phylo_dict[n][d] = {}
            diversity_vs_diversity_phylo_dict[n][d]['n_taxa'] = len(coarse_grained_n)
            diversity_vs_diversity_phylo_dict[n][d]['fraction_of_taxa_coarse_grained'] = len(coarse_grained_n)/len(species_n)
            diversity_vs_diversity_phylo_dict[n][d]['diversity_all'] = diversity_species.tolist()
            diversity_vs_diversity_phylo_dict[n][d]['diversity_coarse'] = diversity_coarse.tolist()
            diversity_vs_diversity_phylo_dict[n][d]['slope'] = slope
            diversity_vs_diversity_phylo_dict[n][d]['percentile'] = percentile
            diversity_vs_diversity_phylo_dict[n][d]['lower_ci'] = lower_ci
            diversity_vs_diversity_phylo_dict[n][d]['upper_ci'] = upper_ci




    diversity_vs_diversity_phylo_path = diversity_vs_diversity_phylo_resource_template % ('_resource')

    sys.stderr.write("Saving diversity vs. diversity dictionary...\n")
    with open(diversity_vs_diversity_phylo_path, 'wb') as handle:
        pickle.dump(diversity_vs_diversity_phylo_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)






def load_diversity_vs_diversity_phylo_resource_dict():

    diversity_vs_diversity_phylo_path = diversity_vs_diversity_phylo_resource_template % ('_resource')

    with open(diversity_vs_diversity_phylo_path, 'rb') as handle:
        diversity_vs_diversity_phylo_dict = pickle.load(handle)
    return diversity_vs_diversity_phylo_dict



def load_diversity_vs_diversity_phylo_emp_dict(environment, rarefied):

    diversity_vs_diversity_phylo_path = diversity_vs_diversity_phylo_emp_template % ('_' + diversity_utils.get_environment_label(environment) + diversity_utils.get_rarefied_label(rarefied))
    with open(diversity_vs_diversity_phylo_path, 'rb') as handle:
        diversity_vs_diversity_phylo_dict = pickle.load(handle)
    return diversity_vs_diversity_phylo_dict










def make_permutation_rescaled_abundance_taxon_dict(block_richnes = 100, iter = 1000):

    taxa_ranks = diversity_utils.taxa_ranks
    idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
    taxa_ranks_label = diversity_utils.taxa_ranks_label

    diversity_dict = {}
    for environment in diversity_utils.environments_to_keep:

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        sad_annotated_dict = load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied=False)
        samples = sad_annotated_dict['samples']

        sys.stderr.write("Building site-by-species matrix...\n")
        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
        richness = len(taxa)

        # figure out number of blocks
        n_blocks = math.ceil(richness / block_richnes)
        block_richness_all = [block_richnes]*(n_blocks-1)
        # add remaining species
        block_richness_all.append(richness - block_richnes*(n_blocks-1))
        block_richness_all = numpy.asarray(block_richness_all)

        rel_s_by_s = s_by_s/numpy.sum(s_by_s,axis=0)
        mean_abundance = numpy.mean(rel_s_by_s, axis=1)

        # create array of block numbers
        # sort rel_s_by_s by mean abundance
        s_idx = list(range(s_by_s.shape[0]))
        mean_abundance_s_idx_tuple_all = list(zip(mean_abundance.tolist(), s_idx))

        # sort tuples by mean abundance, order of increasing abundane
        mean_abundance_s_idx_tuple_all.sort(key=lambda tup: tup[0])
        s_sort_idx = numpy.asarray([s[1] for s in mean_abundance_s_idx_tuple_all])
        mean_abundance_sort = numpy.asarray([s[0] for s in mean_abundance_s_idx_tuple_all])

        s_by_s_sort = s_by_s[s_sort_idx,:]
        rel_s_by_s_sort = rel_s_by_s[s_sort_idx,:]
        taxa_sort = taxa[s_sort_idx]
        scaled_rel_s_by_s_sort = (rel_s_by_s_sort.T/mean_abundance_sort).T

        # get coarse_grained
        rank_coarse_grain_dict = {}
        diversity_dict[environment] = {}
        for rank_idx, rank in enumerate(taxa_ranks):

            sys.stderr.write("Starting %s level analysis...\n" % rank)
            all_genera = numpy.asarray(list(set([sad_annotated_dict['taxa'][t][rank] for t in taxa_sort])))
            all_genera_idx = numpy.arange(len(all_genera))

            genus_to_taxa_dict = {}
            sad_genera_all = []
            #coarse_grained_n = []
            #coarse_grained_list = []
            coarse_grained_idx_all = []
            for genus_idx, genus in enumerate(all_genera):
                genus_to_taxa_dict[genus] = []
                for t in taxa_sort:
                    if sad_annotated_dict['taxa'][t][rank] == genus:
                        genus_to_taxa_dict[genus].append(t)

                g_taxa_idx = numpy.asarray([numpy.where(taxa_sort == g_taxa_i)[0][0] for g_taxa_i in genus_to_taxa_dict[genus]])
                s_by_s_sort_genus = s_by_s_sort[g_taxa_idx,:]
                #coarse_grained_n.append(len(g_taxa_idx))
                sad_genera_all.append(s_by_s_sort_genus.sum(axis=0))
                coarse_grained_idx_all.append(g_taxa_idx)


            #coarse_grained_n = numpy.asarray(coarse_grained_n)
            # get indexes for each clade
            #coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa_sort==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            # re-sort the data
            coarse_grained_idx_all_flat = []
            for c in coarse_grained_idx_all:
                coarse_grained_idx_all_flat.extend(c.tolist())

            coarse_grained_idx_all_flat = numpy.asarray(coarse_grained_idx_all_flat)
            n_coarse_grained_all = numpy.asarray([len(c) for c in coarse_grained_idx_all])
            coarse_grain_idx = numpy.append([0], numpy.cumsum(n_coarse_grained_all))[:-1]

            rank_coarse_grain_dict[rank] = coarse_grain_idx

            # get diversity
            rel_s_by_s_sort_coarse = numpy.add.reduceat(rel_s_by_s_sort, coarse_grain_idx, axis=0)
            diversity = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_s_by_s_sort_coarse)

            diversity_dict[environment][rank] = {}
            diversity_dict[environment][rank]['observed'] = numpy.mean(diversity)
            diversity_dict[environment][rank]['null_constran_mean_and_occupancy_all'] = []
            diversity_dict[environment][rank]['null_constran_mean_all'] = []
            diversity_dict[environment][rank]['null_all'] = []


        for i in range(iter):

            if (i%100 == 0) and (i>0):
                print(i)

            # shuffle unconstrained
            numpy.random.shuffle(rel_s_by_s_sort)

            # shuffle constrained on mean and occupancy
            scaled_rel_s_by_s_sort_copy = numpy.copy(scaled_rel_s_by_s_sort)

            numpy.random.shuffle(block_richness_all)
            permutation_block_idx = numpy.append([0], numpy.cumsum(block_richness_all))[:-1]
            # add final index
            permutation_block_idx = numpy.append(permutation_block_idx, richness)

            # shuffle each block
            for n in range(n_blocks):
                numpy.random.shuffle(scaled_rel_s_by_s_sort_copy[permutation_block_idx[n]:permutation_block_idx[n+1],:])

            unscaled_rel_s_by_s_sort_copy = ((scaled_rel_s_by_s_sort_copy.T)*mean_abundance_sort).T
            # renormalize in the off-change rows sum to more than one
            unscaled_rel_s_by_s_sort_copy = unscaled_rel_s_by_s_sort_copy/numpy.sum(unscaled_rel_s_by_s_sort_copy,axis=0)

            # shuffle constrained on just the mean
            scaled_rel_s_by_s_sort_just_mean_copy = numpy.copy(scaled_rel_s_by_s_sort)
            numpy.random.shuffle(scaled_rel_s_by_s_sort_just_mean_copy)
            unscaled_rel_s_by_s_sort_just_mean_copy = ((scaled_rel_s_by_s_sort_just_mean_copy.T)*mean_abundance_sort).T


            # coarse-grain
            for rank in taxa_ranks:

                unscaled_rel_s_by_s_sort_copy_coarse = numpy.add.reduceat(unscaled_rel_s_by_s_sort_copy, rank_coarse_grain_dict[rank], axis=0)
                diversity_null_constrain_mean_and_occupancy = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, unscaled_rel_s_by_s_sort_copy_coarse)
                diversity_dict[environment][rank]['null_constran_mean_and_occupancy_all'].append(numpy.mean(diversity_null_constrain_mean_and_occupancy))

                rel_s_by_s_sort_coarse = numpy.add.reduceat(rel_s_by_s_sort, rank_coarse_grain_dict[rank], axis=0)
                diversity_null_unconstrained = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, rel_s_by_s_sort_coarse)
                diversity_dict[environment][rank]['null_all'].append(numpy.mean(diversity_null_unconstrained))

                unscaled_rel_s_by_s_sort_copy_constrain_mean_coarse = numpy.add.reduceat(unscaled_rel_s_by_s_sort_just_mean_copy, rank_coarse_grain_dict[rank], axis=0)
                diversity_null_constrain_mean = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, unscaled_rel_s_by_s_sort_copy_constrain_mean_coarse)
                diversity_dict[environment][rank]['null_constran_mean_all'].append(numpy.mean(diversity_null_constrain_mean))


        # get confidence intervals
        for rank in taxa_ranks:

            for measure in ['null_constran_mean_and_occupancy_all', 'null_all', 'null_constran_mean_all']:

                data = numpy.sort(numpy.asarray(diversity_dict[environment][rank][measure]))
                lower_ci_data = data[int((0.025*iter))]
                upper_ci_data = data[int((0.975*iter))]

                del diversity_dict[environment][rank][measure]
                diversity_dict[environment][rank]['%s_lower_ci'%  measure] = lower_ci_data
                diversity_dict[environment][rank]['%s_upper_ci'%  measure] = upper_ci_data


    permutation_rescaled_abundance_taxon_emp_path = permutation_rescaled_abundance_taxon_emp_template % block_richnes

    sys.stderr.write("Saving diversity dictionary...\n")
    with open(permutation_rescaled_abundance_taxon_emp_path, 'wb') as handle:
        pickle.dump(diversity_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_permutation_rescaled_abundance_taxon_dict(block_richnes=100):

    permutation_rescaled_abundance_taxon_emp_path = permutation_rescaled_abundance_taxon_emp_template % block_richnes

    with open(permutation_rescaled_abundance_taxon_emp_path, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_













def load_permutation_rescaled_abundance_phylo_dict(block_richnes=100):

    permutation_rescaled_abundance_phylo_emp_path = permutation_rescaled_abundance_phylo_emp_template % block_richnes

    with open(permutation_rescaled_abundance_phylo_emp_path, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_




def make_coarse_grained_tree_no_subsampling_all_dict():

    coarse_grained_tree_dict_all = {}
    distances_all = []
    for environment in diversity_utils.environments_to_keep:
        sys.stderr.write("Loading tree dict for %s...\n" % environment)
        coarse_grained_tree_dict = load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=False)
        coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict

    return coarse_grained_tree_dict_all









def predict_richness_dbd(focal_mean_fine_rel_s_by_s, focal_var_fine_rel_s_by_s, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, nonfocal_var_coarse_rel_s_by_s, total_reads_coarse):

    focal_fine_beta = (focal_mean_fine_rel_s_by_s**2)/focal_var_fine_rel_s_by_s
    nonfocal_coarse_beta = (nonfocal_mean_coarse_rel_s_by_s**2)/nonfocal_var_coarse_rel_s_by_s
   
    focal_fine_theta = focal_mean_fine_rel_s_by_s/focal_fine_beta
    nonfocal_coarse_theta = nonfocal_mean_coarse_rel_s_by_s/nonfocal_coarse_beta

    focal_fine_richness_predicted = numpy.asarray([sum(1-((1+focal_fine_theta*totreads_i)**(-1*focal_fine_beta))) for totreads_i in total_reads_fine])
    nonfocal_coarse_richness_predicted = numpy.asarray([sum(1-((1+nonfocal_coarse_theta*totreads_i)**(-1*nonfocal_coarse_beta))) for totreads_i in total_reads_fine])

    return focal_fine_richness_predicted, nonfocal_coarse_richness_predicted


def predict_richness_dbd_with_beta(focal_mean_fine_rel_s_by_s, focal_fine_beta, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, nonfocal_coarse_beta, total_reads_coarse):

    focal_fine_theta = focal_mean_fine_rel_s_by_s/focal_fine_beta
    nonfocal_coarse_theta = nonfocal_mean_coarse_rel_s_by_s/nonfocal_coarse_beta

    focal_fine_richness_predicted = numpy.asarray([sum(1-((1+focal_fine_theta*totreads_i)**(-1*focal_fine_beta))) for totreads_i in total_reads_fine])
    nonfocal_coarse_richness_predicted = numpy.asarray([sum(1-((1+nonfocal_coarse_theta*totreads_i)**(-1*nonfocal_coarse_beta))) for totreads_i in total_reads_fine])

    return focal_fine_richness_predicted, nonfocal_coarse_richness_predicted






def predict_diversity_dbd(focal_mean_fine_rel_s_by_s, focal_var_fine_rel_s_by_s, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, nonfocal_var_coarse_rel_s_by_s, total_reads_coarse):

    focal_fine_beta = (focal_mean_fine_rel_s_by_s**2)/focal_var_fine_rel_s_by_s
    nonfocal_coarse_beta = (nonfocal_mean_coarse_rel_s_by_s**2)/nonfocal_var_coarse_rel_s_by_s
   
    focal_fine_theta = focal_mean_fine_rel_s_by_s/focal_fine_beta
    nonfocal_coarse_theta = nonfocal_mean_coarse_rel_s_by_s/nonfocal_coarse_beta

    focal_fine_richness_predicted = numpy.asarray([sum(1-((1+focal_fine_theta*totreads_i)**(-1*focal_fine_beta))) for totreads_i in total_reads_fine])
    nonfocal_coarse_richness_predicted = numpy.asarray([sum(1-((1+nonfocal_coarse_theta*totreads_i)**(-1*nonfocal_coarse_beta))) for totreads_i in total_reads_fine])

    return focal_fine_richness_predicted, nonfocal_coarse_richness_predicted


def predict_diversity_dbd_with_beta(focal_mean_fine_rel_s_by_s, focal_fine_beta, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, nonfocal_coarse_beta, total_reads_coarse):

    focal_fine_theta = focal_mean_fine_rel_s_by_s/focal_fine_beta
    nonfocal_coarse_theta = nonfocal_mean_coarse_rel_s_by_s/nonfocal_coarse_beta

    focal_fine_richness_predicted = numpy.asarray([sum(1-((1+focal_fine_theta*totreads_i)**(-1*focal_fine_beta))) for totreads_i in total_reads_fine])
    nonfocal_coarse_richness_predicted = numpy.asarray([sum(1-((1+nonfocal_coarse_theta*totreads_i)**(-1*nonfocal_coarse_beta))) for totreads_i in total_reads_fine])

    return focal_fine_richness_predicted, nonfocal_coarse_richness_predicted









def make_richness_slm_dbd_dict():

    coarse_grained_tree_dict = make_coarse_grained_tree_no_subsampling_all_dict()

    dbd_dict = {}

    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        sys.stderr.write("Loading taxon dict for %s...\n" % environment)
        sad_annotated_dict = load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied=False)
        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        dbd_dict[environment] = {}
        dbd_dict[environment]['taxon'] = {}
        dbd_dict[environment]['phylo'] = {}

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

            total_reads_fine = numpy.sum(s_by_s_fine, axis=0)
            rel_s_by_s_fine =  s_by_s_fine/numpy.sum(s_by_s_fine, axis=0)
            mean_rel_s_by_s_fine = numpy.mean(rel_s_by_s_fine, axis=1)
            var_rel_s_by_s_fine = numpy.var(rel_s_by_s_fine, axis=1)

            coarse_idx = numpy.append([0], numpy.cumsum(counts_coarse))[:-1]
            s_by_s_coarse = numpy.add.reduceat(s_by_s_fine, coarse_idx, axis=0)       

            total_reads_coarse = numpy.sum(s_by_s_coarse, axis=0)
            rel_s_by_s_coarse = s_by_s_coarse/numpy.sum(s_by_s_coarse, axis=0)
            mean_rel_s_by_s_coarse = numpy.mean(rel_s_by_s_coarse, axis=1)
            var_rel_s_by_s_coarse = numpy.var(rel_s_by_s_coarse, axis=1)


            # mean CV^2
            beta_fine = (mean_rel_s_by_s_fine**2)/var_rel_s_by_s_fine
            mean_mean_rel_s_by_s_fine = numpy.mean(mean_rel_s_by_s_fine)
            mean_beta_fine = numpy.mean(beta_fine)

            beta_coarse = (mean_rel_s_by_s_coarse**2)/var_rel_s_by_s_coarse
            mean_mean_rel_s_by_s_coarse = numpy.mean(mean_rel_s_by_s_coarse)
            mean_beta_coarse = numpy.mean(beta_coarse)


            # get index for fine for null
            #unique_fine, counts_fine = numpy.unique(all_fine, return_counts=True)
            #fine_idx = numpy.append([0], numpy.cumsum(counts_fine))[:-1]

            dbd_dict[environment]['taxon'][coarse_rank] = {}
            dbd_dict[environment]['taxon'][coarse_rank]['focal_coarse'] = {}

            slope_all = []
            slope_slm_all = []
            slope_slm_fix_beta_all = []
            slope_slm_fix_mean_all = []
            slope_slm_fix_beta_across_all = []
            slope_slm_fix_mean_across_all = []
            for focal_coarse_idx, focal_coarse in enumerate(set_coarse):

                # ignore coarse-grained taxa with less than five fine-grained taxa
                if counts_coarse[focal_coarse_idx] < 5:
                    continue

                # all the fine-scale indices for the focal
                focal_s_by_s_coarse_idx = numpy.asarray(idx_to_sort[focal_coarse_idx])
                focal_fine_s_by_s = s_by_s_fine[focal_s_by_s_coarse_idx,:]

                focal_mean_fine_rel_s_by_s = mean_rel_s_by_s_fine[focal_s_by_s_coarse_idx]
                focal_var_fine_rel_s_by_s = var_rel_s_by_s_fine[focal_s_by_s_coarse_idx]

                nonfocal_s_by_s_coarse = numpy.delete(s_by_s_coarse, focal_coarse_idx, axis=0)
                nonfocal_mean_coarse_rel_s_by_s = numpy.delete(mean_rel_s_by_s_coarse, focal_coarse_idx, axis=0)          
                nonfocal_var_coarse_rel_s_by_s = numpy.delete(var_rel_s_by_s_coarse, focal_coarse_idx, axis=0)      

                focal_fine_richness_predicted, nonfocal_coarse_richness_predicted = predict_richness_dbd(focal_mean_fine_rel_s_by_s, focal_var_fine_rel_s_by_s, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, nonfocal_var_coarse_rel_s_by_s, total_reads_coarse)    

                # fix CV
                mean_beta_fine_fix_array = numpy.asarray([mean_beta_fine]*len(focal_mean_fine_rel_s_by_s))
                mean_beta_coarse_fix_array = numpy.asarray([mean_beta_coarse]*len(nonfocal_mean_coarse_rel_s_by_s))
                focal_fine_richness_predicted_fix_beta, nonfocal_coarse_richness_predicted_fix_beta = predict_richness_dbd_with_beta(focal_mean_fine_rel_s_by_s, mean_beta_fine_fix_array, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, mean_beta_coarse_fix_array, total_reads_coarse)

                # fix CV same across scales
                mean_beta_fine_fix_across_array = numpy.asarray([mean_beta_fine]*len(focal_mean_fine_rel_s_by_s))
                mean_beta_coarse_fix_across_array = numpy.asarray([mean_beta_fine]*len(nonfocal_mean_coarse_rel_s_by_s))
                focal_fine_richness_predicted_fix_beta_across, nonfocal_coarse_richness_predicted_fix_beta_across = predict_richness_dbd_with_beta(focal_mean_fine_rel_s_by_s, mean_beta_fine_fix_across_array, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, mean_beta_coarse_fix_across_array, total_reads_coarse)


                # fix mean
                beta_fine_array = (focal_mean_fine_rel_s_by_s**2)/focal_var_fine_rel_s_by_s
                beta_coarse_array = (nonfocal_mean_coarse_rel_s_by_s**2)/nonfocal_var_coarse_rel_s_by_s
                mean_mean_rel_s_by_s_fine_array = numpy.asarray([mean_mean_rel_s_by_s_fine]*len(focal_mean_fine_rel_s_by_s))
                mean_mean_rel_s_by_s_coarse_array = numpy.asarray([mean_mean_rel_s_by_s_coarse]*len(nonfocal_mean_coarse_rel_s_by_s))
                focal_fine_richness_predicted_fix_mean, nonfocal_coarse_richness_predicted_fix_mean = predict_richness_dbd_with_beta(mean_mean_rel_s_by_s_fine_array, beta_fine_array, total_reads_fine, mean_mean_rel_s_by_s_coarse_array, beta_coarse_array, total_reads_coarse)


                # fix mean same across scales
                mean_mean_rel_s_by_s_fine_across_array = numpy.asarray([mean_mean_rel_s_by_s_fine]*len(focal_mean_fine_rel_s_by_s))
                mean_mean_rel_s_by_s_coarse_across_array = numpy.asarray([mean_mean_rel_s_by_s_fine]*len(nonfocal_mean_coarse_rel_s_by_s))
                focal_fine_richness_predicted_fix_mean_across, nonfocal_coarse_richness_predicted_fix_mean_across = predict_richness_dbd_with_beta(mean_mean_rel_s_by_s_fine_across_array, beta_fine_array, total_reads_fine, mean_mean_rel_s_by_s_coarse_across_array, beta_coarse_array, total_reads_coarse)


                #n_fine_focal_coarse = len(focal_coarse_s_by_s_idx)
                measure_coarse = numpy.sum(nonfocal_s_by_s_coarse > 0, axis=0)
                measure_fine = numpy.sum(focal_fine_s_by_s > 0, axis=0)

                # remove sites where there are no observations in either focal or non-focal
                #idx_to_remove = (measure_fine>0) | (measure_coarse>0)

                #measure_coarse = measure_coarse[idx_to_remove]
                #measure_fine = measure_fine[idx_to_remove]

                #nonfocal_coarse_richness_predicted = nonfocal_coarse_richness_predicted[idx_to_remove]
                #focal_fine_richness_predicted = focal_fine_richness_predicted[idx_to_remove]

                #nonfocal_coarse_richness_predicted_fix_beta = nonfocal_coarse_richness_predicted_fix_beta[idx_to_remove]
                #focal_fine_richness_predicted_fix_beta = focal_fine_richness_predicted_fix_beta[idx_to_remove]

                #focal_fine_richness_predicted_fix_beta_across = focal_fine_richness_predicted_fix_beta_across[idx_to_remove]
                #nonfocal_coarse_richness_predicted_fix_beta_across = nonfocal_coarse_richness_predicted_fix_beta_across[idx_to_remove]

                #nonfocal_coarse_richness_predicted_fix_mean = nonfocal_coarse_richness_predicted_fix_mean[idx_to_remove]
                #focal_fine_richness_predicted_fix_mean = focal_fine_richness_predicted_fix_mean[idx_to_remove]

                if len(measure_fine) < 5:
                    continue


                slope, intercept, r_value, p_value, std_err = stats.linregress(measure_coarse, measure_fine)
                slope_slm, intercept_slm, r_value_slm, p_value_slm, std_err_slm = stats.linregress(nonfocal_coarse_richness_predicted, focal_fine_richness_predicted)
                slope_slm_fix_beta, intercept_slm_fix_beta, r_value_slm_fix_beta, p_value_slm_fix_beta, std_err_slm_fix_beta = stats.linregress(nonfocal_coarse_richness_predicted_fix_beta, focal_fine_richness_predicted_fix_beta)
                slope_slm_fix_mean, intercept_slm_fix_mean, r_value_slm_fix_mean, p_value_slm_fix_mean, std_err_slm_fix_mean = stats.linregress(nonfocal_coarse_richness_predicted_fix_mean, focal_fine_richness_predicted_fix_mean)
                slope_slm_fix_beta_across, intercept_slm_fix_beta, r_value_slm_fix_beta, p_value_slm_fix_beta, std_err_slm_fix_beta = stats.linregress(nonfocal_coarse_richness_predicted_fix_beta_across, focal_fine_richness_predicted_fix_beta_across)
                slope_slm_fix_mean_across, intercept_slm_fix_beta, r_value_slm_fix_beta, p_value_slm_fix_beta, std_err_slm_fix_beta = stats.linregress(nonfocal_coarse_richness_predicted_fix_mean_across, focal_fine_richness_predicted_fix_mean_across)


                slope_all.append(slope)
                slope_slm_all.append(slope_slm)
                slope_slm_fix_beta_all.append(slope_slm_fix_beta)
                slope_slm_fix_mean_all.append(slope_slm_fix_mean)
                slope_slm_fix_beta_across_all.append(slope_slm_fix_beta_across)
                slope_slm_fix_mean_across_all.append(slope_slm_fix_mean_across)


                dbd_dict[environment]['taxon'][coarse_rank]['focal_coarse'][focal_coarse] = {}
                dbd_dict[environment]['taxon'][coarse_rank]['focal_coarse'][focal_coarse]['measure_coarse'] = measure_coarse.tolist()
                dbd_dict[environment]['taxon'][coarse_rank]['focal_coarse'][focal_coarse]['measure_fine'] = measure_fine.tolist()
                dbd_dict[environment]['taxon'][coarse_rank]['focal_coarse'][focal_coarse]['measure_coarse_predicted'] = nonfocal_coarse_richness_predicted.tolist()
                dbd_dict[environment]['taxon'][coarse_rank]['focal_coarse'][focal_coarse]['measure_fine_predicted'] = focal_fine_richness_predicted.tolist()



                dbd_dict[environment]['taxon'][coarse_rank]['focal_coarse'][focal_coarse]['slope'] = slope
                dbd_dict[environment]['taxon'][coarse_rank]['focal_coarse'][focal_coarse]['slope_slm'] = slope_slm

        

            slope_all = numpy.asarray(slope_all)
            slope_slm_all = numpy.asarray(slope_slm_all)
            slope_slm_fix_beta_all = numpy.asarray(slope_slm_fix_beta_all)
            slope_slm_fix_mean_all = numpy.asarray(slope_slm_fix_mean_all)
            slope_slm_fix_beta_across_all = numpy.asarray(slope_slm_fix_beta_across_all)
            slope_slm_fix_mean_across_all = numpy.asarray(slope_slm_fix_mean_across_all)

            idx_to_keep = ~(numpy.isnan(slope_all) | numpy.isnan(slope_slm_all) | numpy.isnan(slope_slm_fix_beta_all) | numpy.isnan(slope_slm_fix_mean_all) | numpy.isnan(slope_slm_fix_beta_across_all) | numpy.isnan(slope_slm_fix_mean_across_all) )  

            if sum(idx_to_keep) < 3:
                continue

            slope_all = slope_all[idx_to_keep]
            slope_slm_all = slope_slm_all[idx_to_keep]
            slope_slm_fix_beta_all = slope_slm_fix_beta_all[idx_to_keep]
            slope_slm_fix_mean_all = slope_slm_fix_mean_all[idx_to_keep]
            slope_slm_fix_beta_across_all = slope_slm_fix_beta_across_all[idx_to_keep]
            slope_slm_fix_mean_across_all = slope_slm_fix_mean_across_all[idx_to_keep]

            mean_slope_all = numpy.mean(slope_all)
            mean_slope_slm_all = numpy.mean(slope_slm_all)
            mean_slope_slm_fix_beta_all = numpy.mean(slope_slm_fix_beta_all)
            mean_slope_slm_fix_mean_all = numpy.mean(slope_slm_fix_mean_all)
            mean_slope_slm_fix_beta_across_all = numpy.mean(slope_slm_fix_beta_across_all)
            mean_slope_slm_fix_mean_across_all = numpy.mean(slope_slm_fix_mean_across_all)

            # error
            mean_error_slope_slm = numpy.mean(numpy.absolute( slope_slm_all - slope_all) / numpy.absolute(slope_all))
            mean_error_slope_slm_fix_beta = numpy.mean(numpy.absolute( slope_slm_fix_beta_all - slope_all) / numpy.absolute(slope_all))
            mean_error_slope_slm_fix_mean = numpy.mean(numpy.absolute( slope_slm_fix_mean_all - slope_all) / numpy.absolute(slope_all))
            mean_error_slope_slm_fix_beta_across = numpy.mean(numpy.absolute( slope_slm_fix_beta_across_all - slope_all) / numpy.absolute(slope_all))
            mean_error_slope_slm_fix_mean_across = numpy.mean(numpy.absolute( slope_slm_fix_mean_across_all - slope_all) / numpy.absolute(slope_all))


            #dbd_dict[environment]['taxon'][coarse_rank] = {}
            dbd_dict[environment]['taxon'][coarse_rank]['slope_all'] = slope_all.tolist()
            dbd_dict[environment]['taxon'][coarse_rank]['slope_slm_all'] = slope_slm_all.tolist()
            dbd_dict[environment]['taxon'][coarse_rank]['slope_slm_fix_beta_all'] = slope_slm_fix_beta_all.tolist()
            dbd_dict[environment]['taxon'][coarse_rank]['slope_slm_fix_mean_all'] = slope_slm_fix_mean_all.tolist()
            dbd_dict[environment]['taxon'][coarse_rank]['slope_slm_fix_beta_across_all'] = slope_slm_fix_beta_across_all.tolist()
            dbd_dict[environment]['taxon'][coarse_rank]['slope_slm_fix_mean_across_all'] = slope_slm_fix_mean_across_all.tolist()

            dbd_dict[environment]['taxon'][coarse_rank]['mean_slope'] = mean_slope_all
            dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm'] = mean_slope_slm_all
            dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm_fix_beta'] = mean_slope_slm_fix_beta_all
            dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm_fix_mean'] = mean_slope_slm_fix_mean_all
            dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm_fix_beta_across'] = mean_slope_slm_fix_beta_across_all
            dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm_fix_mean_across'] = mean_slope_slm_fix_mean_across_all

            dbd_dict[environment]['taxon'][coarse_rank]['mean_error_slope_slm'] = mean_error_slope_slm
            dbd_dict[environment]['taxon'][coarse_rank]['mean_error_slope_slm_fix_beta'] = mean_error_slope_slm_fix_beta
            dbd_dict[environment]['taxon'][coarse_rank]['mean_error_slope_slm_fix_mean'] = mean_error_slope_slm_fix_mean
            dbd_dict[environment]['taxon'][coarse_rank]['mean_error_slope_slm_fix_beta_across'] = mean_error_slope_slm_fix_beta_across
            dbd_dict[environment]['taxon'][coarse_rank]['mean_error_slope_slm_fix_mean_across'] = mean_error_slope_slm_fix_mean_across


        # phylogenetic coarse-graining
        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        coarse_grained_tree_dict_env = coarse_grained_tree_dict[environment]
        sad_annotated_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        distances = list(coarse_grained_tree_dict_env.keys())
        distances.sort()

        distances = numpy.asarray(distances)

        for distance_idx in range(len(distances) - 1):
            
            fine_distance = distances[distance_idx]
            coarse_distance = distances[distance_idx + 1]

            sys.stderr.write("Phylo distance = %s \n" % round(coarse_distance, 7))

            fine_grained_list = coarse_grained_tree_dict_env[fine_distance]
            coarse_grained_list = coarse_grained_tree_dict_env[coarse_distance]

            coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
            fine_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in fine_grained_list_i]) for fine_grained_list_i in fine_grained_list]

            # coarse grain s-by-s for all clades
            s_by_s_coarse = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
            s_by_s_fine = numpy.stack([numpy.sum(s_by_s[fine_grained_idx,:], axis=0) for fine_grained_idx in fine_grained_idx_all], axis=0)

            total_reads_fine = numpy.sum(s_by_s_fine, axis=0)
            rel_s_by_s_fine =  s_by_s_fine/numpy.sum(s_by_s_fine, axis=0)
            mean_rel_s_by_s_fine = numpy.mean(rel_s_by_s_fine, axis=1)
            var_rel_s_by_s_fine = numpy.var(rel_s_by_s_fine, axis=1)

            total_reads_coarse = numpy.sum(s_by_s_coarse, axis=0)
            rel_s_by_s_coarse = s_by_s_coarse/numpy.sum(s_by_s_coarse, axis=0)
            mean_rel_s_by_s_coarse = numpy.mean(rel_s_by_s_coarse, axis=1)
            var_rel_s_by_s_coarse = numpy.var(rel_s_by_s_coarse, axis=1)

            # mean CV^2
            beta_fine = (mean_rel_s_by_s_fine**2)/var_rel_s_by_s_fine
            mean_mean_rel_s_by_s_fine = numpy.mean(mean_rel_s_by_s_fine)
            mean_beta_fine = numpy.mean(beta_fine)

            beta_coarse = (mean_rel_s_by_s_coarse**2)/var_rel_s_by_s_coarse
            mean_mean_rel_s_by_s_coarse = numpy.mean(mean_rel_s_by_s_coarse)
            mean_beta_coarse = numpy.mean(beta_coarse)


            slope_all = []
            slope_slm_all = []
            slope_slm_fix_beta_all = []
            slope_slm_fix_mean_all = []
            slope_slm_fix_beta_across_all = []
            slope_slm_fix_mean_across_all = []
            for focal_coarse_idx, focal_coarse in enumerate(coarse_grained_list):

                # ignore coarse-grained taxa with less than five fine-grained taxa
                if len(focal_coarse) < 5:
                    continue

                fine_in_coarse_all = []
                # identify fine-grain taxa containing coarse-grained members
                f_idx_all = []
                for f_idx, f in enumerate(fine_grained_list):

                    is_fine_in_coarse = bool(set(focal_coarse) & set(f))

                    if is_fine_in_coarse == True:
                        fine_in_coarse_all.extend(f)
                        f_idx_all.append(f_idx)

                f_idx_all = numpy.asarray(f_idx_all)
                # ignore coarse-grained taxa with less than five fine-grained taxa
                if len(f_idx_all) < 5:
                    continue

                s_by_s_fine_focal = s_by_s_fine[f_idx_all,:]
                s_by_s_coarse_nonfocal = numpy.delete(s_by_s_coarse, focal_coarse_idx, axis=0)

                richness_fine_focal = numpy.sum(s_by_s_fine_focal>0, axis=0)
                richness_coarse_nonfocal = numpy.sum(s_by_s_coarse_nonfocal>0, axis=0)

                # predict using SLM
                focal_mean_fine_rel_s_by_s = mean_rel_s_by_s_fine[f_idx_all]
                focal_var_fine_rel_s_by_s = var_rel_s_by_s_fine[f_idx_all]

                nonfocal_mean_coarse_rel_s_by_s = numpy.delete(mean_rel_s_by_s_coarse, focal_coarse_idx, axis=0)          
                nonfocal_var_coarse_rel_s_by_s = numpy.delete(var_rel_s_by_s_coarse, focal_coarse_idx, axis=0)      

                focal_fine_richness_predicted, nonfocal_coarse_richness_predicted = predict_richness_dbd(focal_mean_fine_rel_s_by_s, focal_var_fine_rel_s_by_s, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, nonfocal_var_coarse_rel_s_by_s, total_reads_coarse)    
                
                # fix CV
                mean_beta_fine_array = numpy.asarray([mean_beta_fine]*len(focal_mean_fine_rel_s_by_s))
                mean_beta_coarse_array = numpy.asarray([mean_beta_coarse]*len(nonfocal_mean_coarse_rel_s_by_s))
                focal_fine_richness_predicted_fix_beta, nonfocal_coarse_richness_predicted_fix_beta = predict_richness_dbd_with_beta(focal_mean_fine_rel_s_by_s, mean_beta_fine_array, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, mean_beta_coarse_array, total_reads_coarse)

                # fix CV same across scales
                mean_beta_fine_fix_across_array = numpy.asarray([mean_beta_fine]*len(focal_mean_fine_rel_s_by_s))
                mean_beta_coarse_fix_across_array = numpy.asarray([mean_beta_fine]*len(nonfocal_mean_coarse_rel_s_by_s))
                focal_fine_richness_predicted_fix_beta_across, nonfocal_coarse_richness_predicted_fix_beta_across = predict_richness_dbd_with_beta(focal_mean_fine_rel_s_by_s, mean_beta_fine_fix_across_array, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, mean_beta_coarse_fix_across_array, total_reads_coarse)

                # fix mean
                beta_fine_array = (focal_mean_fine_rel_s_by_s**2)/focal_var_fine_rel_s_by_s
                beta_coarse_array = (nonfocal_mean_coarse_rel_s_by_s**2)/nonfocal_var_coarse_rel_s_by_s
                mean_mean_rel_s_by_s_fine_array = numpy.asarray([mean_mean_rel_s_by_s_fine]*len(focal_mean_fine_rel_s_by_s))
                mean_mean_rel_s_by_s_coarse_array = numpy.asarray([mean_mean_rel_s_by_s_coarse]*len(nonfocal_mean_coarse_rel_s_by_s))
                focal_fine_richness_predicted_fix_mean, nonfocal_coarse_richness_predicted_fix_mean = predict_richness_dbd_with_beta(mean_mean_rel_s_by_s_fine_array, beta_fine_array, total_reads_fine, mean_mean_rel_s_by_s_coarse_array, beta_coarse_array, total_reads_coarse)

                # fix mean same across scales
                mean_mean_rel_s_by_s_fine_across_array = numpy.asarray([mean_mean_rel_s_by_s_fine]*len(focal_mean_fine_rel_s_by_s))
                mean_mean_rel_s_by_s_coarse_across_array = numpy.asarray([mean_mean_rel_s_by_s_fine]*len(nonfocal_mean_coarse_rel_s_by_s))
                focal_fine_richness_predicted_fix_mean_across, nonfocal_coarse_richness_predicted_fix_mean_across = predict_richness_dbd_with_beta(mean_mean_rel_s_by_s_fine_across_array, beta_fine_array, total_reads_fine, mean_mean_rel_s_by_s_coarse_across_array, beta_coarse_array, total_reads_coarse)

                measure_coarse = numpy.sum(nonfocal_s_by_s_coarse > 0, axis=0)
                measure_fine = numpy.sum(focal_fine_s_by_s > 0, axis=0)

                if len(measure_fine) < 5:
                    continue

                slope, intercept, r_value, p_value, std_err = stats.linregress(richness_coarse_nonfocal, richness_fine_focal)
                slope_slm, intercept_slm, r_value_slm, p_value_slm, std_err_slm = stats.linregress(nonfocal_coarse_richness_predicted, focal_fine_richness_predicted)
                slope_slm_fix_beta, intercept_slm_fix_beta, r_value_slm_fix_beta, p_value_slm_fix_beta, std_err_slm_fix_beta = stats.linregress(nonfocal_coarse_richness_predicted_fix_beta, focal_fine_richness_predicted_fix_beta)
                slope_slm_fix_mean, intercept_slm_fix_mean, r_value_slm_fix_mean, p_value_slm_fix_mean, std_err_slm_fix_mean = stats.linregress(nonfocal_coarse_richness_predicted_fix_mean, focal_fine_richness_predicted_fix_mean)
                slope_slm_fix_beta_across, intercept_slm_fix_beta, r_value_slm_fix_beta, p_value_slm_fix_beta, std_err_slm_fix_beta = stats.linregress(nonfocal_coarse_richness_predicted_fix_beta_across, focal_fine_richness_predicted_fix_beta_across)
                slope_slm_fix_mean_across, intercept_slm_fix_beta, r_value_slm_fix_beta, p_value_slm_fix_beta, std_err_slm_fix_beta = stats.linregress(nonfocal_coarse_richness_predicted_fix_mean_across, focal_fine_richness_predicted_fix_mean_across)

                slope_all.append(slope)
                slope_slm_all.append(slope_slm)
                slope_slm_fix_beta_all.append(slope_slm_fix_beta)
                slope_slm_fix_mean_all.append(slope_slm_fix_mean)
                slope_slm_fix_beta_across_all.append(slope_slm_fix_beta_across)
                slope_slm_fix_mean_across_all.append(slope_slm_fix_mean_across)


            slope_all = numpy.asarray(slope_all)
            slope_slm_all = numpy.asarray(slope_slm_all)
            slope_slm_fix_beta_all = numpy.asarray(slope_slm_fix_beta_all)
            slope_slm_fix_mean_all = numpy.asarray(slope_slm_fix_mean_all)
            slope_slm_fix_beta_across_all = numpy.asarray(slope_slm_fix_beta_across_all)
            slope_slm_fix_mean_across_all = numpy.asarray(slope_slm_fix_mean_across_all)

            idx_to_keep = ~(numpy.isnan(slope_all) | numpy.isnan(slope_slm_all) | numpy.isnan(slope_slm_fix_beta_all) | numpy.isnan(slope_slm_fix_mean_all) | numpy.isnan(slope_slm_fix_beta_across_all) | numpy.isnan(slope_slm_fix_mean_across_all) )

            if sum(idx_to_keep) < 3:
                continue
            
            slope_all = slope_all[idx_to_keep]
            slope_slm_all = slope_slm_all[idx_to_keep]
            slope_slm_fix_beta_all = slope_slm_fix_beta_all[idx_to_keep]
            slope_slm_fix_mean_all = slope_slm_fix_mean_all[idx_to_keep]
            slope_slm_fix_beta_across_all = slope_slm_fix_beta_across_all[idx_to_keep]
            slope_slm_fix_mean_across_all = slope_slm_fix_mean_across_all[idx_to_keep]

            mean_slope_all = numpy.mean(slope_all)
            mean_slope_slm_all = numpy.mean(slope_slm_all)
            mean_slope_slm_fix_beta_all = numpy.mean(slope_slm_fix_beta_all)
            mean_slope_slm_fix_mean_all = numpy.mean(slope_slm_fix_mean_all)
            mean_slope_slm_fix_beta_across_all = numpy.mean(slope_slm_fix_beta_across_all)
            mean_slope_slm_fix_mean_across_all = numpy.mean(slope_slm_fix_mean_across_all)

            # error
            mean_error_slope_slm = numpy.mean(numpy.absolute( slope_slm_all - slope_all) / numpy.absolute(slope_all))
            mean_error_slope_slm_fix_beta = numpy.mean(numpy.absolute( slope_slm_fix_beta_all - slope_all) / numpy.absolute(slope_all))
            mean_error_slope_slm_fix_mean = numpy.mean(numpy.absolute( slope_slm_fix_mean_all - slope_all) / numpy.absolute(slope_all))
            mean_error_slope_slm_fix_beta_across = numpy.mean(numpy.absolute( slope_slm_fix_beta_across_all - slope_all) / numpy.absolute(slope_all))
            mean_error_slope_slm_fix_mean_across = numpy.mean(numpy.absolute( slope_slm_fix_mean_across_all - slope_all) / numpy.absolute(slope_all))
            
            dbd_dict[environment]['phylo'][coarse_distance] = {}
            dbd_dict[environment]['phylo'][coarse_distance]['slope_all'] = slope_all.tolist()
            dbd_dict[environment]['phylo'][coarse_distance]['slope_slm_all'] = slope_slm_all.tolist()
            dbd_dict[environment]['phylo'][coarse_distance]['slope_slm_fix_beta_all'] = slope_slm_fix_beta_all.tolist()
            dbd_dict[environment]['phylo'][coarse_distance]['slope_slm_fix_mean_all'] = slope_slm_fix_mean_all.tolist()
            dbd_dict[environment]['phylo'][coarse_distance]['slope_slm_fix_beta_across_all'] = slope_slm_fix_beta_across_all.tolist()
            dbd_dict[environment]['phylo'][coarse_distance]['slope_slm_fix_mean_across_all'] = slope_slm_fix_mean_across_all.tolist()

            dbd_dict[environment]['phylo'][coarse_distance]['mean_slope'] = mean_slope_all
            dbd_dict[environment]['phylo'][coarse_distance]['mean_slope_slm'] = mean_slope_slm_all
            dbd_dict[environment]['phylo'][coarse_distance]['mean_slope_slm_fix_beta'] = mean_slope_slm_fix_beta_all
            dbd_dict[environment]['phylo'][coarse_distance]['mean_slope_slm_fix_mean'] = mean_slope_slm_fix_mean_all
            dbd_dict[environment]['phylo'][coarse_distance]['mean_slope_slm_fix_beta_across'] = mean_slope_slm_fix_beta_across_all
            dbd_dict[environment]['phylo'][coarse_distance]['mean_slope_slm_fix_mean_across'] = mean_slope_slm_fix_mean_across_all

            dbd_dict[environment]['phylo'][coarse_distance]['mean_error_slope_slm'] = mean_error_slope_slm
            dbd_dict[environment]['phylo'][coarse_distance]['mean_error_slope_slm_fix_beta'] = mean_error_slope_slm_fix_beta
            dbd_dict[environment]['phylo'][coarse_distance]['mean_error_slope_slm_fix_mean'] = mean_error_slope_slm_fix_mean
            dbd_dict[environment]['phylo'][coarse_distance]['mean_error_slope_slm_fix_beta_across'] = mean_error_slope_slm_fix_beta_across
            dbd_dict[environment]['phylo'][coarse_distance]['mean_error_slope_slm_fix_mean_across'] = mean_error_slope_slm_fix_mean_across



    sys.stderr.write("Saving diversity dictionary...\n")
    with open(richness_dbd_dict_path, 'wb') as handle:
        pickle.dump(dbd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    



def load_richness_slm_dbd_dict():

    with open(richness_dbd_dict_path, 'rb') as handle:
        dbd_slope_dict = pickle.load(handle)
    return dbd_slope_dict




def make_diversity_slm_phlyo_integral_dict(environment):

    environments_to_keep = diversity_utils.environments_to_keep

    coarse_grained_tree_dict_all = {}
    distances_all = []
    for e in environments_to_keep:
        sys.stderr.write("Loading tree dict for %s...\n" % e)
        coarse_grained_tree_dict = load_coarse_grained_tree_no_subsampling_dict(environment=e, rarefied=False)
        coarse_grained_tree_dict_all[e] = coarse_grained_tree_dict
        distances = list(coarse_grained_tree_dict.keys())
        distances_all.extend(distances)

    distances_all = list(set(distances_all))
    distances_all.sort()

    coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    pres_abs_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
    samples = pres_abs_dict['samples']
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
    
    sys.stderr.write("Getting site-by-species matrix...\n")
    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

    total_reads = numpy.sum(s_by_s, axis=0)
    
    distances = list(coarse_grained_tree_dict.keys())
    distances.sort()
    distances = numpy.asarray(distances)
    #distances = [0.21199020238496094]

    dict_ = {}    
    dict_['phylo'] = {}
    for distance_idx, distance in enumerate(distances):

        if distance not in coarse_grained_tree_dict:
            continue

        sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
        
        coarse_grained_list = coarse_grained_tree_dict[distance]
        #coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

        # get indexes for each clade
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_coarse = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
        rel_s_by_s_coarse =  s_by_s_coarse/numpy.sum(s_by_s_coarse, axis=0)
        mean_rel_s_by_s_coarse = numpy.mean(rel_s_by_s_coarse, axis=1)
        var_rel_s_by_s_coarse = numpy.var(rel_s_by_s_coarse, axis=1)
        beta_coarse = (mean_rel_s_by_s_coarse**2)/var_rel_s_by_s_coarse

        # first moment
        # second moment
        first_moment_integral_all = []
        second_moment_integral_all = []
        for i in range(len(mean_rel_s_by_s_coarse)):
            
            mean_rel_s_by_s_coarse_i = mean_rel_s_by_s_coarse[i]
            beta_coarse_i = beta_coarse[i]

            first_moment_integral_i_all = []
            second_moment_integral_i_all = []

            for m in range(len(total_reads)):

                total_reads_m = total_reads[m]

                first_moment_integral_i_m = integrate.quad(diversity_utils.integrand_first_moment, 0, total_reads_m, args=(total_reads_m, mean_rel_s_by_s_coarse_i, beta_coarse_i))[0]
                second_moment_integral_i_m = integrate.quad(diversity_utils.integrand_second_moment, 0, total_reads_m, args=(total_reads_m, mean_rel_s_by_s_coarse_i, beta_coarse_i))[0]

                first_moment_integral_i_all.append(first_moment_integral_i_m)
                second_moment_integral_i_all.append(second_moment_integral_i_m)

            first_moment_integral_all.append(first_moment_integral_i_all)
            second_moment_integral_all.append(second_moment_integral_i_all)

        
        dict_['phylo'][distance] = {}
        dict_['phylo'][distance]['first_moment'] = first_moment_integral_all
        dict_['phylo'][distance]['second_moment'] = second_moment_integral_all


    diversity_slm_phylo_integral_dict_path_ = diversity_slm_phylo_integral_dict_path % diversity_utils.get_environment_label(environment)

    sys.stderr.write("Saving integral dictionary...\n")
    with open(diversity_slm_phylo_integral_dict_path_, 'wb') as handle:
        pickle.dump(dict_, handle, protocol=pickle.HIGHEST_PROTOCOL)




def load_diversity_slm_phlyo_integral_dict(environment):

    diversity_slm_phylo_integral_dict_path_ = diversity_slm_phylo_integral_dict_path % diversity_utils.get_environment_label(environment)

    with open(diversity_slm_phylo_integral_dict_path_, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_





def make_diversity_slm_phlyo_integral_otu_level_dict(environment):

    '''
    Solves integrals for each OTU for all OTUs used in phylogenetic coarse-graining
    '''

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    pres_abs_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
    samples = pres_abs_dict['samples']
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
    
    sys.stderr.write("Getting site-by-species matrix...\n")
    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

    rel_s_by_s =  s_by_s/numpy.sum(s_by_s, axis=0)
    mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
    total_reads = numpy.sum(s_by_s, axis=0)
    beta_coarse = (mean_rel_s_by_s**2)/var_rel_s_by_s
    
    dict_ = {}
    sys.stderr.write("Solving integrals...\n")
    for t_idx, t in enumerate(taxa):
        dict_[t] = {}

        mean_rel_s_by_s_t = mean_rel_s_by_s[t_idx]
        beta_t = beta_coarse[t_idx]

        # first moment
        # second moment
        first_moment_integral_t_all = []
        second_moment_integral_t_all = []

        for m in range(len(total_reads)):

            total_reads_m = total_reads[m]

            first_moment_integral_t_m = integrate.quad(diversity_utils.integrand_first_moment, 0, total_reads_m, args=(total_reads_m, mean_rel_s_by_s_t, beta_t), epsabs=1e-20)[0]
            second_moment_integral_t_m = integrate.quad(diversity_utils.integrand_second_moment, 0, total_reads_m, args=(total_reads_m, mean_rel_s_by_s_t, beta_t), epsabs=1e-25)[0]

            first_moment_integral_t_all.append(first_moment_integral_t_m)
            second_moment_integral_t_all.append(second_moment_integral_t_m)

        dict_[t]['first_moment'] = first_moment_integral_t_all
        dict_[t]['second_moment'] = second_moment_integral_t_all


    diversity_slm_phylo_integral_dict_path_ = diversity_slm_phylo_integral_otu_level_dict_path % diversity_utils.get_environment_label(environment)

    sys.stderr.write("Saving integral dictionary...\n")
    with open(diversity_slm_phylo_integral_dict_path_, 'wb') as handle:
        pickle.dump(dict_, handle, protocol=pickle.HIGHEST_PROTOCOL)




def load_diversity_slm_phlyo_integral_otu_level_dict(environment):

    diversity_slm_phylo_integral_dict_path_ = diversity_slm_phylo_integral_otu_level_dict_path % diversity_utils.get_environment_label(environment)

    with open(diversity_slm_phylo_integral_dict_path_, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_







def make_diversity_slm_dbd_integral_dict(environment):

    coarse_grained_tree_dict = make_coarse_grained_tree_no_subsampling_all_dict()

    #for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    sys.stderr.write("Loading taxon dict for %s...\n" % environment)
    sad_annotated_dict = load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied=False)
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    dbd_dict = {}
    dbd_dict['taxon'] = {}
    dbd_dict['phylo'] = {}

    # each coarse-graining threshold refers to a rank that's compared
    # within each rank there's a fine and coarse grained
    for coarse_rank_idx, coarse_rank in enumerate(diversity_utils.taxa_ranks):

        sys.stderr.write("Starting %s level analysis...\n" % coarse_rank)

        if coarse_rank == 'genus':
            fine_rank = 'asv'
            all_fine = numpy.copy(taxa).tolist()
        else:
            fine_rank = diversity_utils.taxa_ranks[coarse_rank_idx-1]
            all_fine = list([sad_annotated_dict['taxa'][t][fine_rank] for t in taxa])

        sys.stderr.write("Solving integrals...\n")

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

        total_reads_fine = numpy.sum(s_by_s_fine, axis=0)
        rel_s_by_s_fine =  s_by_s_fine/numpy.sum(s_by_s_fine, axis=0)
        mean_rel_s_by_s_fine = numpy.mean(rel_s_by_s_fine, axis=1)
        var_rel_s_by_s_fine = numpy.var(rel_s_by_s_fine, axis=1)

        coarse_idx = numpy.append([0], numpy.cumsum(counts_coarse))[:-1]
        s_by_s_coarse = numpy.add.reduceat(s_by_s_fine, coarse_idx, axis=0)       

        total_reads_coarse = numpy.sum(s_by_s_coarse, axis=0)
        rel_s_by_s_coarse = s_by_s_coarse/numpy.sum(s_by_s_coarse, axis=0)
        mean_rel_s_by_s_coarse = numpy.mean(rel_s_by_s_coarse, axis=1)
        var_rel_s_by_s_coarse = numpy.var(rel_s_by_s_coarse, axis=1)


        # mean CV^2
        beta_fine = (mean_rel_s_by_s_fine**2)/var_rel_s_by_s_fine
        beta_coarse = (mean_rel_s_by_s_coarse**2)/var_rel_s_by_s_coarse


        # get index for fine for null
        #unique_fine, counts_fine = numpy.unique(all_fine, return_counts=True)
        #fine_idx = numpy.append([0], numpy.cumsum(counts_fine))[:-1]
        # make matrix of predicted x*logx
        coarse_integral_all = []
        for i in range(len(mean_rel_s_by_s_coarse)):
            
            mean_rel_s_by_s_coarse_i = mean_rel_s_by_s_coarse[i]
            beta_coarse_i = beta_coarse[i]

            integral_coarse_i_all = []

            for m in range(len(total_reads_coarse)):

                total_reads_coarse_m = total_reads_coarse[m]

                integral_coarse_i_m = integrate.quad(diversity_utils.integrand_first_moment, 0, total_reads_coarse_m, args=(total_reads_coarse_m, mean_rel_s_by_s_coarse_i, beta_coarse_i))[0]
                integral_coarse_i_all.append(integral_coarse_i_m)

            coarse_integral_all.append(integral_coarse_i_all)


        fine_integral_all = []
        for i in range(len(mean_rel_s_by_s_fine)):
            
            mean_rel_s_by_s_fine_i = mean_rel_s_by_s_fine[i]
            beta_fine_i = beta_fine[i]

            integral_fine_i_all = []
            for m in range(len(total_reads_fine)):

                total_reads_fine_m = total_reads_fine[m]

                integral_fine_i_m = integrate.quad(diversity_utils.integrand_first_moment, 0, total_reads_fine_m, args=(total_reads_fine_m, mean_rel_s_by_s_fine_i, beta_fine_i))[0]
                integral_fine_i_all.append(integral_fine_i_m)

            fine_integral_all.append(integral_fine_i_all)

        

        dbd_dict['taxon'][coarse_rank] = {}
        dbd_dict['taxon'][coarse_rank]['coarse_integral'] = coarse_integral_all
        dbd_dict['taxon'][coarse_rank]['fine_integral'] = fine_integral_all


    # phylogenetic coarse-graining
    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    coarse_grained_tree_dict_env = coarse_grained_tree_dict[environment]
    sad_annotated_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
    samples = sad_annotated_dict['samples']

    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    distances = list(coarse_grained_tree_dict_env.keys())
    distances.sort()

    distances = numpy.asarray(distances)

    for distance_idx in range(len(distances) - 1):
        
        fine_distance = distances[distance_idx]
        coarse_distance = distances[distance_idx + 1]

        sys.stderr.write("Phylo distance = %s \n" % round(coarse_distance, 7))

        fine_grained_list = coarse_grained_tree_dict_env[fine_distance]
        coarse_grained_list = coarse_grained_tree_dict_env[coarse_distance]

        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
        fine_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in fine_grained_list_i]) for fine_grained_list_i in fine_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_coarse = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
        s_by_s_fine = numpy.stack([numpy.sum(s_by_s[fine_grained_idx,:], axis=0) for fine_grained_idx in fine_grained_idx_all], axis=0)

        total_reads_fine = numpy.sum(s_by_s_fine, axis=0)
        rel_s_by_s_fine =  s_by_s_fine/numpy.sum(s_by_s_fine, axis=0)
        mean_rel_s_by_s_fine = numpy.mean(rel_s_by_s_fine, axis=1)
        var_rel_s_by_s_fine = numpy.var(rel_s_by_s_fine, axis=1)

        total_reads_coarse = numpy.sum(s_by_s_coarse, axis=0)
        rel_s_by_s_coarse = s_by_s_coarse/numpy.sum(s_by_s_coarse, axis=0)
        mean_rel_s_by_s_coarse = numpy.mean(rel_s_by_s_coarse, axis=1)
        var_rel_s_by_s_coarse = numpy.var(rel_s_by_s_coarse, axis=1)

        # mean CV^2
        beta_fine = (mean_rel_s_by_s_fine**2)/var_rel_s_by_s_fine
        beta_coarse = (mean_rel_s_by_s_coarse**2)/var_rel_s_by_s_coarse


        sys.stderr.write("Solving integrals...\n")
        coarse_integral_all = []
        for i in range(len(mean_rel_s_by_s_coarse)):
            
            mean_rel_s_by_s_coarse_i = mean_rel_s_by_s_coarse[i]
            beta_coarse_i = beta_coarse[i]

            integral_coarse_i_all = []

            for m in range(len(total_reads_coarse)):

                total_reads_coarse_m = total_reads_coarse[m]

                integral_coarse_i_m = integrate.quad(diversity_utils.integrand_first_moment, 0, total_reads_coarse_m, args=(total_reads_coarse_m, mean_rel_s_by_s_coarse_i, beta_coarse_i))[0]
                integral_coarse_i_all.append(integral_coarse_i_m)

            coarse_integral_all.append(integral_coarse_i_all)


        fine_integral_all = []
        for i in range(len(mean_rel_s_by_s_fine)):
            
            mean_rel_s_by_s_fine_i = mean_rel_s_by_s_fine[i]
            beta_fine_i = beta_fine[i]

            integral_fine_i_all = []
            for m in range(len(total_reads_fine)):

                total_reads_fine_m = total_reads_fine[m]

                integral_fine_i_m = integrate.quad(diversity_utils.integrand_first_moment, 0, total_reads_fine_m, args=(total_reads_fine_m, mean_rel_s_by_s_fine_i, beta_fine_i))[0]
                integral_fine_i_all.append(integral_fine_i_m)

            fine_integral_all.append(integral_fine_i_all)

        
        dbd_dict['phylo'][coarse_distance] = {}
        dbd_dict['phylo'][coarse_distance]['coarse_integral'] = coarse_integral_all
        dbd_dict['phylo'][coarse_distance]['fine_integral'] = fine_integral_all



    diversity_slm_dbd_integral_dict_path_ = diversity_slm_dbd_integral_dict_path % diversity_utils.get_environment_label(environment)

    sys.stderr.write("Saving diversity dictionary...\n")
    with open(diversity_slm_dbd_integral_dict_path_, 'wb') as handle:
        pickle.dump(dbd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

        

def load_diversity_slm_dbd_integral_dict(environment):

    diversity_slm_dbd_integral_dict_path_ = diversity_slm_dbd_integral_dict_path % diversity_utils.get_environment_label(environment)

    with open(diversity_slm_dbd_integral_dict_path_, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_






def make_diversity_slm_dbd_dict(environment):

    coarse_grained_tree_dict = make_coarse_grained_tree_no_subsampling_all_dict()

    dbd_dict = {}

    sys.stderr.write("Loading taxon dict for %s...\n" % environment)
    sad_annotated_dict = load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied=False)
    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    dbd_dict['taxon'] = {}
    dbd_dict['phylo'] = {}

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

        total_reads_fine = numpy.sum(s_by_s_fine, axis=0)
        rel_s_by_s_fine =  s_by_s_fine/numpy.sum(s_by_s_fine, axis=0)
        mean_rel_s_by_s_fine = numpy.mean(rel_s_by_s_fine, axis=1)
        var_rel_s_by_s_fine = numpy.var(rel_s_by_s_fine, axis=1)

        coarse_idx = numpy.append([0], numpy.cumsum(counts_coarse))[:-1]
        s_by_s_coarse = numpy.add.reduceat(s_by_s_fine, coarse_idx, axis=0)       

        total_reads_coarse = numpy.sum(s_by_s_coarse, axis=0)
        rel_s_by_s_coarse = s_by_s_coarse/numpy.sum(s_by_s_coarse, axis=0)
        mean_rel_s_by_s_coarse = numpy.mean(rel_s_by_s_coarse, axis=1)
        var_rel_s_by_s_coarse = numpy.var(rel_s_by_s_coarse, axis=1)


        # mean CV^2
        beta_fine = (mean_rel_s_by_s_fine**2)/var_rel_s_by_s_fine
        mean_mean_rel_s_by_s_fine = numpy.mean(mean_rel_s_by_s_fine)
        mean_beta_fine = numpy.mean(beta_fine)

        beta_coarse = (mean_rel_s_by_s_coarse**2)/var_rel_s_by_s_coarse
        mean_mean_rel_s_by_s_coarse = numpy.mean(mean_rel_s_by_s_coarse)
        mean_beta_coarse = numpy.mean(beta_coarse)

        # get index for fine for null
        #unique_fine, counts_fine = numpy.unique(all_fine, return_counts=True)
        #fine_idx = numpy.append([0], numpy.cumsum(counts_fine))[:-1]
        # make matrix of predicted x*logx
        coarse_integral_all = []
        for i in range(len(mean_rel_s_by_s_coarse)):
            
            mean_rel_s_by_s_coarse_i = mean_rel_s_by_s_coarse[i]
            beta_coarse_i = beta_coarse[i]

            integral_coarse_i_all = []

            for m in range(len(total_reads_coarse)):

                total_reads_coarse_m = total_reads_coarse[m]

                integral_coarse_i_m = integrate.quad(diversity_utils.integrand_first_moment, 0, total_reads_coarse_m, args=(total_reads_coarse_m, mean_rel_s_by_s_coarse_i, beta_coarse_i))[0]
                integral_coarse_i_all.append(integral_coarse_i_m)

            coarse_integral_all.append(integral_coarse_i_all)

        coarse_integral_all = numpy.array(coarse_integral_all)


        fine_integral_all = []
        for i in range(len(mean_rel_s_by_s_fine)):
            
            mean_rel_s_by_s_fine_i = mean_rel_s_by_s_fine[i]
            beta_fine_i = beta_fine[i]

            integral_fine_i_all = []
            for m in range(len(total_reads_fine)):

                total_reads_fine_m = total_reads_fine[m]

                integral_fine_i_m = integrate.quad(diversity_utils.integrand_first_moment, 0, total_reads_fine_m, args=(total_reads_fine_m, mean_rel_s_by_s_fine_i, beta_fine_i))[0]
                integral_fine_i_all.append(integral_fine_i_m)

            fine_integral_all.append(integral_fine_i_all)


        fine_integral_all = numpy.array(fine_integral_all)

        slope_all = []
        slope_slm_all = []
        
        dbd_dict['taxon'][coarse_rank] = {}
        dbd_dict['taxon'][coarse_rank]['focal_coarse'] = {}
        for focal_coarse_idx, focal_coarse in enumerate(set_coarse):

            # ignore coarse-grained taxa with less than five fine-grained taxa
            if counts_coarse[focal_coarse_idx] < 5:
                continue

            # all the fine-scale indices for the focal
            focal_s_by_s_coarse_idx = numpy.asarray(idx_to_sort[focal_coarse_idx])
            focal_fine_s_by_s = s_by_s_fine[focal_s_by_s_coarse_idx,:]
            nonfocal_s_by_s_coarse = numpy.delete(s_by_s_coarse, focal_coarse_idx, axis=0)

            # predict
            focal_fine_s_by_s_predicted = fine_integral_all[focal_s_by_s_coarse_idx,:]
            nonfocal_s_by_s_coarse_predicted = numpy.delete(coarse_integral_all, focal_coarse_idx, axis=0)

            focal_fine_diversity_predicted = -1*numpy.sum(focal_fine_s_by_s_predicted, axis=0)
            nonfocal_coarse_diversity_predicted = -1*numpy.sum(nonfocal_s_by_s_coarse_predicted, axis=0)

            measure_coarse = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, nonfocal_s_by_s_coarse)
            measure_fine = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, focal_fine_s_by_s)

            if len(measure_fine) < 5:
                continue
            

            slope, intercept, r_value, p_value, std_err = stats.linregress(measure_coarse, measure_fine)
            slope_slm, intercept_slm, r_value_slm, p_value_slm, std_err_slm = stats.linregress(nonfocal_coarse_diversity_predicted, focal_fine_diversity_predicted)

            slope_all.append(slope)
            slope_slm_all.append(slope_slm)

            dbd_dict['taxon'][coarse_rank]['focal_coarse'][focal_coarse] = {}
            dbd_dict['taxon'][coarse_rank]['focal_coarse'][focal_coarse]['measure_coarse'] = measure_coarse.tolist()
            dbd_dict['taxon'][coarse_rank]['focal_coarse'][focal_coarse]['measure_fine'] = measure_fine.tolist()
            dbd_dict['taxon'][coarse_rank]['focal_coarse'][focal_coarse]['measure_coarse_predicted'] = nonfocal_coarse_diversity_predicted.tolist()
            dbd_dict['taxon'][coarse_rank]['focal_coarse'][focal_coarse]['measure_fine_predicted'] = focal_fine_diversity_predicted.tolist()

            dbd_dict['taxon'][coarse_rank]['focal_coarse'][focal_coarse]['slope'] = slope
            dbd_dict['taxon'][coarse_rank]['focal_coarse'][focal_coarse]['slope_slm'] = slope_slm



        slope_all = numpy.asarray(slope_all)
        slope_slm_all = numpy.asarray(slope_slm_all)


        #idx_to_keep = ~(numpy.isnan(slope_all) | numpy.isnan(slope_slm_all) | numpy.isnan(slope_slm_fix_beta_all) | numpy.isnan(slope_slm_fix_mean_all) | numpy.isnan(slope_slm_fix_beta_across_all) | numpy.isnan(slope_slm_fix_mean_across_all) )  
        idx_to_keep = ~(numpy.isnan(slope_all) | numpy.isnan(slope_slm_all) )  

        if sum(idx_to_keep) < 3:
            continue

        slope_all = slope_all[idx_to_keep]
        slope_slm_all = slope_slm_all[idx_to_keep]


        mean_slope_all = numpy.mean(slope_all)
        mean_slope_slm_all = numpy.mean(slope_slm_all)


        # error
        mean_error_slope_slm = numpy.mean(numpy.absolute( slope_slm_all - slope_all) / numpy.absolute(slope_all))

        
        dbd_dict['taxon'][coarse_rank]['slope_all'] = slope_all.tolist()
        dbd_dict['taxon'][coarse_rank]['slope_slm_all'] = slope_slm_all.tolist()

        dbd_dict['taxon'][coarse_rank]['mean_slope'] = mean_slope_all
        dbd_dict['taxon'][coarse_rank]['mean_slope_slm'] = mean_slope_slm_all

        dbd_dict['taxon'][coarse_rank]['mean_error_slope_slm'] = mean_error_slope_slm



    # load integral for phylogeny
    diversity_slm_dbd_integral_dict = load_diversity_slm_dbd_integral_dict(environment)

    # phylogenetic coarse-graining
    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    coarse_grained_tree_dict_env = coarse_grained_tree_dict[environment]
    sad_annotated_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
    samples = sad_annotated_dict['samples']

    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    distances = list(coarse_grained_tree_dict_env.keys())
    distances.sort()

    distances = numpy.asarray(distances)

    for distance_idx in range(len(distances) - 1):
        
        fine_distance = distances[distance_idx]
        coarse_distance = distances[distance_idx + 1]

        sys.stderr.write("Phylo distance = %s \n" % round(coarse_distance, 7))

        coarse_integral = diversity_slm_dbd_integral_dict['phylo'][coarse_distance]['coarse_integral']
        fine_integral = diversity_slm_dbd_integral_dict['phylo'][coarse_distance]['fine_integral']

        coarse_integral = numpy.asarray(coarse_integral)
        fine_integral = numpy.asarray(fine_integral)

        coarse_grained_list = coarse_grained_tree_dict_env[coarse_distance]
        fine_grained_list = coarse_grained_tree_dict_env[fine_distance]

        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
        fine_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in fine_grained_list_i]) for fine_grained_list_i in fine_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_coarse = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
        s_by_s_fine = numpy.stack([numpy.sum(s_by_s[fine_grained_idx,:], axis=0) for fine_grained_idx in fine_grained_idx_all], axis=0)


        slope_all = []
        slope_slm_all = []
        for focal_coarse_idx, focal_coarse in enumerate(coarse_grained_list):

            # ignore coarse-grained taxa with less than five fine-grained taxa
            if len(focal_coarse) < 5:
                continue

            fine_in_coarse_all = []
            # identify fine-grain taxa containing coarse-grained members
            f_idx_all = []
            for f_idx, f in enumerate(fine_grained_list):

                is_fine_in_coarse = bool(set(focal_coarse) & set(f))

                if is_fine_in_coarse == True:
                    fine_in_coarse_all.extend(f)
                    f_idx_all.append(f_idx)

            
            f_idx_all = numpy.asarray(f_idx_all)
            
            # ignore coarse-grained taxa with less than five fine-grained taxa
            if len(f_idx_all) < 5:
                continue

            s_by_s_fine_focal = s_by_s_fine[f_idx_all,:]
            s_by_s_coarse_nonfocal = numpy.delete(s_by_s_coarse, focal_coarse_idx, axis=0)

            measure_coarse = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_coarse_nonfocal)
            measure_fine = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_fine_focal)

                # predict
            focal_fine_s_by_s_predicted = fine_integral[f_idx_all,:]
            nonfocal_s_by_s_coarse_predicted = numpy.delete(coarse_integral, focal_coarse_idx, axis=0)                

            focal_fine_diversity_predicted = -1*numpy.sum(focal_fine_s_by_s_predicted, axis=0)
            nonfocal_coarse_diversity_predicted = -1*numpy.sum(nonfocal_s_by_s_coarse_predicted, axis=0)

            if len(measure_fine) < 5:
                continue

            slope, intercept, r_value, p_value, std_err = stats.linregress(measure_coarse, measure_fine)
            slope_slm, intercept_slm, r_value_slm, p_value_slm, std_err_slm = stats.linregress(nonfocal_coarse_diversity_predicted, focal_fine_diversity_predicted)

            slope_all.append(slope)
            slope_slm_all.append(slope_slm)


        slope_all = numpy.asarray(slope_all)
        slope_slm_all = numpy.asarray(slope_slm_all)

        idx_to_keep = ~(numpy.isnan(slope_all) | numpy.isnan(slope_slm_all))

        if sum(idx_to_keep) < 3:
            continue
        
        slope_all = slope_all[idx_to_keep]
        slope_slm_all = slope_slm_all[idx_to_keep]

        mean_slope_all = numpy.mean(slope_all)
        mean_slope_slm_all = numpy.mean(slope_slm_all)

        # error
        mean_error_slope_slm = numpy.mean(numpy.absolute( slope_slm_all - slope_all) / numpy.absolute(slope_all))

        dbd_dict['phylo'][coarse_distance] = {}
        dbd_dict['phylo'][coarse_distance]['slope_all'] = slope_all.tolist()
        dbd_dict['phylo'][coarse_distance]['slope_slm_all'] = slope_slm_all.tolist()

        dbd_dict['phylo'][coarse_distance]['mean_slope'] = mean_slope_all
        dbd_dict['phylo'][coarse_distance]['mean_slope_slm'] = mean_slope_slm_all

        dbd_dict['phylo'][coarse_distance]['mean_error_slope_slm'] = mean_error_slope_slm



    diversity_dbd_dict_path_ = diversity_dbd_dict_path % environment.replace(' ', '_')

    sys.stderr.write("Saving diversity dictionary...\n")
    with open(diversity_dbd_dict_path_, 'wb') as handle:
        pickle.dump(dbd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    




def load_diversity_slm_dbd_dict():

    dict_final = {}

    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        diversity_dbd_dict_path_ = diversity_dbd_dict_path % environment.replace(' ', '_')

        with open(diversity_dbd_dict_path_, 'rb') as handle:
            dict_ = pickle.load(handle)
        dict_final[environment] = dict_
    
    
    return dict_final




def make_richness_diversity_prediction_taxon_dict(iter_=1000):

    richness_diversity_prediction_dict = {}    

    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        sad_annotated_dict = load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = False)

        samples = sad_annotated_dict['samples']
        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        sys.stderr.write("Getting site-by-species matrix...\n")

        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        richness_diversity_prediction_dict[environment] = {}
        richness_diversity_prediction_dict[environment]['taxon'] = {}

        #for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

        for rank_idx, rank in enumerate(['family']):

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



            mean_richness_observed_gamma_rvs, var_richness_observed_gamma_rvs, mean_richness_predicted_gamma_rvs, var_richness_predicted_gamma_rvs, mean_diversity_observed_gamma_rvs, var_diversity_observed_gamma_rvs, mean_diversity_predicted_gamma_rvs, var_diversity_predicted_gamma_rvs = diversity_utils.predict_mean_and_var_richness_and_diversity_using_gamma_rv(s_by_s_genera, iter_=iter_)


            mean_diversity_observed, mean_diversity_predicted, var_diversity_observed, var_diversity_predicted, var_diversity_predicted_plus_covariance = diversity_utils.predict_mean_and_var_diversity_analytic(s_by_s_genera, s_by_s_genera.shape[0])

            mean_richness_observed, mean_richness_predicted = diversity_utils.predict_mean_richness(s_by_s_genera, s_by_s_genera.shape[0])
            var_richness_observed, var_richness_predicted, var_richness_predicted_plus_covariance = diversity_utils.predict_var_richness(s_by_s_genera, s_by_s_genera.shape[0])
            

            richness_diversity_prediction_dict[environment]['taxon'][rank] = {}

            print(mean_diversity_observed, mean_diversity_predicted, mean_diversity_predicted_gamma_rvs)

            # mean richness
            richness_diversity_prediction_dict[environment]['taxon'][rank]['mean_richness_observed'] = mean_richness_observed
            richness_diversity_prediction_dict[environment]['taxon'][rank]['mean_richness_predicted'] = mean_richness_predicted
            richness_diversity_prediction_dict[environment]['taxon'][rank]['mean_richness_predicted_gamma'] = mean_richness_predicted_gamma_rvs

            # var richness
            richness_diversity_prediction_dict[environment]['taxon'][rank]['var_richness_observed'] = var_richness_observed
            richness_diversity_prediction_dict[environment]['taxon'][rank]['var_richness_predicted'] = var_richness_predicted
            richness_diversity_prediction_dict[environment]['taxon'][rank]['var_richness_predicted_plus_covariance'] = var_richness_predicted_plus_covariance
            richness_diversity_prediction_dict[environment]['taxon'][rank]['var_richness_predicted_gamma'] = var_richness_predicted_gamma_rvs

            # mean diversity
            richness_diversity_prediction_dict[environment]['taxon'][rank]['mean_diversity_observed'] = mean_diversity_observed
            richness_diversity_prediction_dict[environment]['taxon'][rank]['mean_diversity_predicted'] = mean_diversity_predicted
            richness_diversity_prediction_dict[environment]['taxon'][rank]['mean_diversity_predicted_gamma'] = mean_diversity_predicted_gamma_rvs

            # var diversity
            richness_diversity_prediction_dict[environment]['taxon'][rank]['var_diversity_observed'] = var_diversity_observed
            richness_diversity_prediction_dict[environment]['taxon'][rank]['var_diversity_predicted'] = var_diversity_predicted
            richness_diversity_prediction_dict[environment]['taxon'][rank]['var_diversity_predicted_plus_covariance'] = var_diversity_predicted_plus_covariance
            richness_diversity_prediction_dict[environment]['taxon'][rank]['var_diversity_predicted_gamma'] = var_diversity_predicted_gamma_rvs


            #relative_error = numpy.absolute(var_diversity_predicted-var_diversity_predicted_gamma_rvs)/var_diversity_predicted_gamma_rvs


    sys.stderr.write("Saving prediction dictionary...\n")
    with open(richness_diversity_prediction_taxon_dict_path, 'wb') as handle:
        pickle.dump(richness_diversity_prediction_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




def load_richness_diversity_prediction_taxon_dict():

    with open(richness_diversity_prediction_taxon_dict_path, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_




# simulate and return simulated values
# do not take expectation
def run_richness_diversity_prediction_gamma_rv_phylo_dict(environment, prefix, iter_=100):

    environments_to_keep = diversity_utils.environments_to_keep

    coarse_grained_tree_dict_all = {}
    distances_all = []
    for e in environments_to_keep:
        sys.stderr.write("Loading tree dict for %s...\n" % e)
        coarse_grained_tree_dict = load_coarse_grained_tree_no_subsampling_dict(environment=e, rarefied=False)
        coarse_grained_tree_dict_all[e] = coarse_grained_tree_dict
        distances = list(coarse_grained_tree_dict.keys())
        distances_all.extend(distances)

    distances_all = list(set(distances_all))
    distances_all.sort()

    coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    pres_abs_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
    samples = pres_abs_dict['samples']
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
    
    sys.stderr.write("Getting site-by-species matrix...\n")
    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])
    
    richness_diversity_prediction_dict = {}    
    richness_diversity_prediction_dict['phylo'] = {}

    for distance_idx, distance in enumerate(distances):

        if distance not in coarse_grained_tree_dict:
            continue

        sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
        
        coarse_grained_list = coarse_grained_tree_dict[distance]
        coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

        # get indexes for each clade
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
        
        mean_richness_observed_gamma_rvs, var_richness_observed_gamma_rvs, mean_richness_predicted_gamma_rvs, var_richness_predicted_gamma_rvs, mean_diversity_observed_gamma_rvs, var_diversity_observed_gamma_rvs, mean_diversity_predicted_gamma_rvs, var_diversity_predicted_gamma_rvs = diversity_utils.predict_mean_and_var_richness_and_diversity_using_gamma_rv(s_by_s_all_clades, iter_=iter_)
        mean_diversity_observed, mean_diversity_predicted, var_diversity_observed, var_diversity_predicted, var_diversity_predicted_plus_covariance = diversity_utils.predict_mean_and_var_diversity_analytic(s_by_s_all_clades, s_by_s_all_clades.shape[0])
        
        mean_richness_observed, mean_richness_predicted = diversity_utils.predict_mean_richness(s_by_s_all_clades, s_by_s_all_clades.shape[0])
        var_richness_observed, var_richness_predicted, var_richness_predicted_plus_covariance = diversity_utils.predict_var_richness(s_by_s_all_clades, s_by_s_all_clades.shape[0])


        #predict_richness_and_diversity_using_gamma_rv






def make_integral_and_richness_diversity_prediction_phylo_dict(environment, iter_=1000):

    environments_to_keep = diversity_utils.environments_to_keep

    coarse_grained_tree_dict_all = {}
    distances_all = []
    for e in environments_to_keep:
        sys.stderr.write("Loading tree dict for %s...\n" % e)
        coarse_grained_tree_dict = load_coarse_grained_tree_no_subsampling_dict(environment=e, rarefied=False)
        coarse_grained_tree_dict_all[e] = coarse_grained_tree_dict
        distances = list(coarse_grained_tree_dict.keys())
        distances_all.extend(distances)

    distances_all = list(set(distances_all))
    distances_all.sort()

    coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    pres_abs_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
    samples = pres_abs_dict['samples']
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
    
    sys.stderr.write("Getting site-by-species matrix...\n")
    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])
    
    richness_diversity_prediction_dict = {}    
    richness_diversity_prediction_dict['phylo'] = {}

    for distance_idx, distance in enumerate(distances):

        if distance not in coarse_grained_tree_dict:
            continue

        sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
        
        coarse_grained_list = coarse_grained_tree_dict[distance]
        coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

        # get indexes for each clade
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
        
        mean_richness_observed_gamma_rvs, var_richness_observed_gamma_rvs, mean_richness_predicted_gamma_rvs, var_richness_predicted_gamma_rvs, mean_diversity_observed_gamma_rvs, var_diversity_observed_gamma_rvs, mean_diversity_predicted_gamma_rvs, var_diversity_predicted_gamma_rvs = diversity_utils.predict_mean_and_var_richness_and_diversity_using_gamma_rv(s_by_s_all_clades, iter_=iter_)
        mean_diversity_observed, mean_diversity_predicted, var_diversity_observed, var_diversity_predicted, var_diversity_predicted_plus_covariance = diversity_utils.predict_mean_and_var_diversity_analytic(s_by_s_all_clades, s_by_s_all_clades.shape[0])
        
        mean_richness_observed, mean_richness_predicted = diversity_utils.predict_mean_richness(s_by_s_all_clades, s_by_s_all_clades.shape[0])
        var_richness_observed, var_richness_predicted, var_richness_predicted_plus_covariance = diversity_utils.predict_var_richness(s_by_s_all_clades, s_by_s_all_clades.shape[0])


        richness_diversity_prediction_dict['phylo'][distance] = {}

        # mean richness
        richness_diversity_prediction_dict['phylo'][distance]['mean_richness_observed'] = mean_richness_observed
        richness_diversity_prediction_dict['phylo'][distance]['mean_richness_predicted'] = mean_richness_predicted
        richness_diversity_prediction_dict['phylo'][distance]['mean_richness_predicted_gamma'] = mean_richness_predicted_gamma_rvs

        # var richness
        richness_diversity_prediction_dict['phylo'][distance]['var_richness_observed'] = var_richness_observed
        richness_diversity_prediction_dict['phylo'][distance]['var_richness_predicted'] = var_richness_predicted
        richness_diversity_prediction_dict['phylo'][distance]['var_richness_predicted_plus_covariance'] = var_richness_predicted_plus_covariance
        richness_diversity_prediction_dict['phylo'][distance]['var_richness_predicted_gamma'] = var_richness_predicted_gamma_rvs

        # mean diversity
        richness_diversity_prediction_dict['phylo'][distance]['mean_diversity_observed'] = mean_diversity_observed
        richness_diversity_prediction_dict['phylo'][distance]['mean_diversity_predicted'] = mean_diversity_predicted
        richness_diversity_prediction_dict['phylo'][distance]['mean_diversity_predicted_gamma'] = mean_diversity_predicted_gamma_rvs

        # var diversity
        richness_diversity_prediction_dict['phylo'][distance]['var_diversity_observed'] = var_diversity_observed
        richness_diversity_prediction_dict['phylo'][distance]['var_diversity_predicted'] = var_diversity_predicted
        richness_diversity_prediction_dict['phylo'][distance]['var_diversity_predicted_plus_covariance'] = var_diversity_predicted_plus_covariance
        richness_diversity_prediction_dict['phylo'][distance]['var_diversity_predicted_gamma'] = var_diversity_predicted_gamma_rvs


    richness_diversity_prediction_phylo_dict_path_ = richness_diversity_prediction_phylo_dict_path % environment.replace(' ', '_')

    sys.stderr.write("Saving prediction dictionary...\n")
    with open(richness_diversity_prediction_phylo_dict_path_, 'wb') as handle:
        pickle.dump(richness_diversity_prediction_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




def make_cov_phylo_dict(environment):

    environments_to_keep = diversity_utils.environments_to_keep

    coarse_grained_tree_dict_all = {}
    distances_all = []
    for e in environments_to_keep:
        sys.stderr.write("Loading tree dict for %s...\n" % e)
        coarse_grained_tree_dict = load_coarse_grained_tree_no_subsampling_dict(environment=e, rarefied=False)
        coarse_grained_tree_dict_all[e] = coarse_grained_tree_dict
        distances = list(coarse_grained_tree_dict.keys())
        distances_all.extend(distances)

    distances_all = list(set(distances_all))
    distances_all.sort()

    coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    pres_abs_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
    
    sys.stderr.write("Getting site-by-species matrix...\n")
    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])
    
    cov_dict = {}    
    cov_dict['phylo'] = {}


    for distance_idx, distance in enumerate(distances):

        if distance not in coarse_grained_tree_dict:
            continue

        sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
        
        coarse_grained_list = coarse_grained_tree_dict[distance]

        # get indexes for each clade
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
        
        rel_s_by_s = (s_by_s_all_clades/s_by_s_all_clades.sum(axis=0))

        x_log_x = rel_s_by_s*numpy.log(rel_s_by_s)
        x_log_x[numpy.isnan(x_log_x)] = 0
        cov_x_log_x = numpy.cov(x_log_x)
        cov_sum = 2*sum(cov_x_log_x[numpy.triu_indices(cov_x_log_x.shape[0], k = 1)])

        cov_dict['phylo']['distance'] = cov_sum


    cov_phylo_dict_path_ = cov_phylo_dict_path % environment.replace(' ', '_')

    sys.stderr.write("Saving prediction dictionary...\n")
    with open(cov_phylo_dict_path_, 'wb') as handle:
        pickle.dump(cov_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)







def load_cov_phylo_dict(environment):

    cov_phylo_dict_path_ = cov_phylo_dict_path % environment.replace(' ', '_')

    with open(cov_phylo_dict_path_, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_






def make_richness_diversity_prediction_phylo_dict(environment, iter_=1000):

    environments_to_keep = diversity_utils.environments_to_keep

    coarse_grained_tree_dict_all = {}
    distances_all = []
    for e in environments_to_keep:
        sys.stderr.write("Loading tree dict for %s...\n" % e)
        coarse_grained_tree_dict = load_coarse_grained_tree_no_subsampling_dict(environment=e, rarefied=False)
        coarse_grained_tree_dict_all[e] = coarse_grained_tree_dict
        distances = list(coarse_grained_tree_dict.keys())
        distances_all.extend(distances)

    distances_all = list(set(distances_all))
    distances_all.sort()

    coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]

    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    pres_abs_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
    samples = pres_abs_dict['samples']
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
    
    sys.stderr.write("Getting site-by-species matrix...\n")
    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])
    
    richness_diversity_prediction_dict = {}    
    richness_diversity_prediction_dict['phylo'] = {}

    sys.stderr.write("Loading integral dict for %s...\n" % environment)
    integral_dict = load_diversity_slm_phlyo_integral_dict(environment)

    for distance_idx, distance in enumerate(distances):

        if distance not in coarse_grained_tree_dict:
            continue

        sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
        
        coarse_grained_list = coarse_grained_tree_dict[distance]
        #coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

        # get indexes for each clade
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
        
        # integral dict
        first_moment_integral_dict = numpy.asarray(integral_dict['phylo'][distance]['first_moment'])
        second_moment_integral_dict = numpy.asarray(integral_dict['phylo'][distance]['second_moment'])

        mean_diversity_observed, mean_diversity_predicted, var_diversity_observed, var_diversity_predicted, var_diversity_predicted_plus_covariance = diversity_utils.predict_mean_and_var_diversity_analytic_with_integral(s_by_s_all_clades, s_by_s_all_clades.shape[0], first_moment_integral_dict, second_moment_integral_dict)
        
        mean_richness_observed_gamma_rvs, var_richness_observed_gamma_rvs, mean_richness_predicted_gamma_rvs, var_richness_predicted_gamma_rvs, mean_diversity_observed_gamma_rvs, var_diversity_observed_gamma_rvs, mean_diversity_predicted_gamma_rvs, var_diversity_predicted_gamma_rvs = diversity_utils.predict_mean_and_var_richness_and_diversity_using_gamma_rv(s_by_s_all_clades, iter_=iter_)
       
        
        mean_richness_observed, mean_richness_predicted = diversity_utils.predict_mean_richness(s_by_s_all_clades, s_by_s_all_clades.shape[0])
        var_richness_observed, var_richness_predicted, var_richness_predicted_plus_covariance = diversity_utils.predict_var_richness(s_by_s_all_clades, s_by_s_all_clades.shape[0])


        richness_diversity_prediction_dict['phylo'][distance] = {}

        # mean richness
        richness_diversity_prediction_dict['phylo'][distance]['mean_richness_observed'] = mean_richness_observed
        richness_diversity_prediction_dict['phylo'][distance]['mean_richness_predicted'] = mean_richness_predicted
        richness_diversity_prediction_dict['phylo'][distance]['mean_richness_predicted_gamma'] = mean_richness_predicted_gamma_rvs

        # var richness
        richness_diversity_prediction_dict['phylo'][distance]['var_richness_observed'] = var_richness_observed
        richness_diversity_prediction_dict['phylo'][distance]['var_richness_predicted'] = var_richness_predicted
        richness_diversity_prediction_dict['phylo'][distance]['var_richness_predicted_plus_covariance'] = var_richness_predicted_plus_covariance
        richness_diversity_prediction_dict['phylo'][distance]['var_richness_predicted_gamma'] = var_richness_predicted_gamma_rvs

        # mean diversity
        richness_diversity_prediction_dict['phylo'][distance]['mean_diversity_observed'] = mean_diversity_observed
        richness_diversity_prediction_dict['phylo'][distance]['mean_diversity_predicted'] = mean_diversity_predicted
        richness_diversity_prediction_dict['phylo'][distance]['mean_diversity_predicted_gamma'] = mean_diversity_predicted_gamma_rvs

        # var diversity
        richness_diversity_prediction_dict['phylo'][distance]['var_diversity_observed'] = var_diversity_observed
        richness_diversity_prediction_dict['phylo'][distance]['var_diversity_predicted'] = var_diversity_predicted
        richness_diversity_prediction_dict['phylo'][distance]['var_diversity_predicted_plus_covariance'] = var_diversity_predicted_plus_covariance
        richness_diversity_prediction_dict['phylo'][distance]['var_diversity_predicted_gamma'] = var_diversity_predicted_gamma_rvs


    richness_diversity_prediction_phylo_dict_path_ = richness_diversity_prediction_phylo_dict_path % environment.replace(' ', '_')

    sys.stderr.write("Saving prediction dictionary...\n")
    with open(richness_diversity_prediction_phylo_dict_path_, 'wb') as handle:
        pickle.dump(richness_diversity_prediction_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)







def load_richness_diversity_prediction_phylo_dict():

    richness_diversity_prediction_dict = {}

    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        if environment == 'soil metagenome':
            continue

        richness_diversity_prediction_phylo_dict_path_ = richness_diversity_prediction_phylo_dict_path % environment.replace(' ', '_')

        with open(richness_diversity_prediction_phylo_dict_path_, 'rb') as handle:
            dict_ = pickle.load(handle)
        richness_diversity_prediction_dict[environment] = dict_
    
    
    return richness_diversity_prediction_dict













def make_diversity_slm_dbd_simulation_taxon_dict(iter_=1000):

    dbd_dict = {}
    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        sys.stderr.write("Loading taxon dict for %s...\n" % environment)
        sad_annotated_dict = load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied=False)
        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        rel_s_by_s = s_by_s/numpy.sum(s_by_s, axis=0)
        mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
        var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
        total_reads_coarse = numpy.sum(s_by_s, axis=0)

        # load SVD
        sys.stderr.write("Loading singular value decomposition dictionary ...\n")
        svd_dict = load_svd_taxon_dict(environment)
        svd = (numpy.asarray(svd_dict['u']), numpy.asarray(svd_dict['s']), numpy.asarray(svd_dict['v']))

        sys.stderr.write("Simulating correlated site-by-OTU read count matrices...\n")
        read_counts_gamma_all = []
        for i in range(iter_):    
            rel_abundances_gamma, read_counts_gamma = simulation_utils.genrate_community_from_mean_and_var_svd(mean_rel_s_by_s, var_rel_s_by_s, total_reads_coarse, len(total_reads_coarse), svd)
            read_counts_gamma_all.append(read_counts_gamma)
    
        # simulate null matrices
        dbd_dict[environment] = {}
        dbd_dict[environment]['taxon'] = {}

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
                s_by_s_simulation_fine = read_counts_gamma_all
                
                all_coarse = []
                for fine in set_fine:
                    all_coarse.append(fine_coarse_dict[fine])
                all_coarse = numpy.asarray(all_coarse)


            else:
                sad_fine_all = []
                sad_fine_simulation_all = [[] for i in range(iter_)]
                all_coarse = []
                for fine in set_fine:
                    taxa_in_fine = []
                    for t in taxa:

                        if sad_annotated_dict['taxa'][t][fine_rank] == fine:
                            taxa_in_fine.append(t)

                    # numpy.where() is a major bottleneck
                    g_taxa_idx = numpy.asarray([numpy.where(taxa == taxa_in_fine_i)[0][0] for taxa_in_fine_i in taxa_in_fine])

                    # sum empirical abundances
                    sad_fine_all.append(s_by_s[g_taxa_idx,:].sum(axis=0))

                    # sum abundances of simulation
                    for i in range(iter_):
                        sad_fine_simulation_all[i].append(read_counts_gamma_all[i][g_taxa_idx,:].sum(axis=0))

                    all_coarse.append(fine_coarse_dict[fine])

                all_coarse = numpy.asarray(all_coarse)
                s_by_s_fine = numpy.stack(sad_fine_all, axis=0)
                s_by_s_simulation_fine = [numpy.stack(i, axis=0) for i in sad_fine_simulation_all]


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
            coarse_idx = numpy.append([0], numpy.cumsum(counts_coarse))[:-1]
            s_by_s_coarse = numpy.add.reduceat(s_by_s_fine, coarse_idx, axis=0)      

            # sort simulated date
            s_by_s_simulation_fine_sorted = []
            s_by_s_simulation_coarse_sorted = []
            for i in range(iter_):
                s_by_s_simulation_fine_i = s_by_s_simulation_fine[i][idx_to_sort_flat,:]
                s_by_s_simulation_fine_sorted.append(s_by_s_simulation_fine_i)
                s_by_s_simulation_coarse_sorted.append(numpy.add.reduceat(s_by_s_simulation_fine_i, coarse_idx, axis=0))

            
            slope_all = []
            slope_slm_all = []
            for focal_coarse_idx, focal_coarse in enumerate(set_coarse):

                # ignore coarse-grained taxa with less than five fine-grained taxa
                if counts_coarse[focal_coarse_idx] < 5:
                    continue
               
                # observed
                # all the fine-scale indices for the focal
                focal_s_by_s_coarse_idx = numpy.asarray(idx_to_sort[focal_coarse_idx])
                focal_fine_s_by_s = s_by_s_fine[focal_s_by_s_coarse_idx,:]

                nonfocal_s_by_s_coarse = numpy.delete(s_by_s_coarse, focal_coarse_idx, axis=0)

                measure_coarse = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, nonfocal_s_by_s_coarse)
                measure_fine = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, focal_fine_s_by_s)
                slope, intercept, r_value, p_value, std_err = stats.linregress(measure_coarse, measure_fine)

                # predicted
                focal_slope_slm_predicted_all = []
                for i in range(iter_):

                    rel_read_counts_gamma_coarse_i = s_by_s_simulation_coarse_sorted[i]
                    rel_read_counts_gamma_fine_i = s_by_s_simulation_fine_sorted[i]

                    focal_fine_s_by_s_gamma_i = rel_read_counts_gamma_fine_i[focal_s_by_s_coarse_idx,:]
                    nonfocal_s_by_s_coarse_gamma_i = numpy.delete(rel_read_counts_gamma_coarse_i, focal_coarse_idx, axis=0)

                    measure_coarse_predicted_i = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, nonfocal_s_by_s_coarse_gamma_i)
                    measure_fine_predicted_i = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, focal_fine_s_by_s_gamma_i)

                    slope_slm, intercept_slm, r_value_slm, p_value_slm, std_err_slm = stats.linregress(measure_coarse_predicted_i, measure_fine_predicted_i)
                    focal_slope_slm_predicted_all.append(slope_slm)

                
                mean_focal_slope_slm_predicted = numpy.mean(focal_slope_slm_predicted_all)

                slope_all.append(slope)
                slope_slm_all.append(mean_focal_slope_slm_predicted)

                #print(slope, mean_focal_slope_slm_predicted)

            slope_all = numpy.asarray(slope_all)
            slope_slm_all = numpy.asarray(slope_slm_all)

            mean_slope_all = numpy.mean(slope_all)
            mean_slope_slm_all = numpy.mean(slope_slm_all)

            mean_error_slope_slm = numpy.mean(numpy.absolute( slope_slm_all - slope_all) / numpy.absolute(slope_all))


            dbd_dict[environment]['taxon'][coarse_rank] = {}
            dbd_dict[environment]['taxon'][coarse_rank]['slope_all'] = slope_all.tolist()
            dbd_dict[environment]['taxon'][coarse_rank]['slope_slm_all'] = slope_slm_all.tolist()

            dbd_dict[environment]['taxon'][coarse_rank]['mean_slope'] = mean_slope_all
            dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm'] = mean_slope_slm_all

            dbd_dict[environment]['taxon'][coarse_rank]['mean_error_slope_slm'] = mean_error_slope_slm


    sys.stderr.write("Saving diversity dictionary...\n")
    with open(diversity_dbd_simulation_taxon_dict_path, 'wb') as handle:
        pickle.dump(dbd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    



def load_diversity_slm_dbd_simulation_taxon_dict():

    with open(diversity_dbd_simulation_taxon_dict_path, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_



def make_svd_taxon_dict():

    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        sys.stderr.write("Loading taxon dict for %s...\n" % environment)
        sad_annotated_dict = load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied=False)
        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        rel_s_by_s = s_by_s/numpy.sum(s_by_s, axis=0)

        corr_rel_s_by_s = numpy.corrcoef(rel_s_by_s)

        u, s, v = numpy.linalg.svd(corr_rel_s_by_s)

        svd_dict = {}
        svd_dict['u'] = u.tolist()
        svd_dict['s'] = s.tolist()
        svd_dict['v'] = v.tolist()

        svd_dict_path_ = svd_taxon_dict_path % environment.replace(' ', '_')

        sys.stderr.write("Saving diversity dictionary...\n")
        with open(svd_dict_path_, 'wb') as handle:
            pickle.dump(svd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        




def load_svd_taxon_dict(environment):

    svd_dict_path_ = svd_taxon_dict_path % environment.replace(' ', '_')

    with open(svd_dict_path_, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_





def make_svd_phylo_dict():

    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        pres_abs_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
        samples = pres_abs_dict['samples']
        taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
        sys.stderr.write("Getting site-by-species matrix...\n")
        s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])

        rel_s_by_s = s_by_s/numpy.sum(s_by_s, axis=0)
        corr_rel_s_by_s = numpy.corrcoef(rel_s_by_s)

        u, s, v = numpy.linalg.svd(corr_rel_s_by_s)

        svd_dict = {}
        svd_dict['u'] = u.tolist()
        svd_dict['s'] = s.tolist()
        svd_dict['v'] = v.tolist()

        svd_dict_path_ = svd_phylo_dict_path % environment.replace(' ', '_')

        sys.stderr.write("Saving diversity dictionary...\n")
        with open(svd_dict_path_, 'wb') as handle:
            pickle.dump(svd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        


def load_svd_phylo_dict(environment):

    svd_dict_path_ = svd_phylo_dict_path % environment.replace(' ', '_')

    with open(svd_dict_path_, 'rb') as handle:
        dict_ = pickle.load(handle)
    return dict_






def run_svd_and_make_diversity_slm_dbd_simulation_phylo_dict(environment, iter_=1000):

    environments_to_keep = diversity_utils.environments_to_keep

    coarse_grained_tree_dict_all = {}
    distances_all = []
    for e in environments_to_keep:
        sys.stderr.write("Loading tree dict for %s...\n" % e)
        coarse_grained_tree_dict = load_coarse_grained_tree_no_subsampling_dict(environment=e, rarefied=False)
        coarse_grained_tree_dict_all[e] = coarse_grained_tree_dict
        distances = list(coarse_grained_tree_dict.keys())
        distances_all.extend(distances)

    distances_all = list(set(distances_all))
    distances_all.sort()

    # phylogenetic coarse-graining
    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    coarse_grained_tree_dict_env = coarse_grained_tree_dict_all[environment]
    sad_annotated_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
    samples = sad_annotated_dict['samples']

    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    # get SVD for simulations
    rel_s_by_s = s_by_s/numpy.sum(s_by_s, axis=0)
    total_reads = numpy.sum(s_by_s, axis=0)
    mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
    corr_rel_s_by_s = numpy.corrcoef(rel_s_by_s)
    u, s, v = numpy.linalg.svd(corr_rel_s_by_s)
    svd = (u, s, v)

    sys.stderr.write("Simulating correlated site-by-OTU read count matrices...\n")
    read_counts_gamma_all = []
    for i in range(iter_):    
        rel_abundances_gamma, read_counts_gamma = simulation_utils.genrate_community_from_mean_and_var_svd(mean_rel_s_by_s, var_rel_s_by_s, total_reads, len(total_reads), svd)
        #read_counts_gamma = numpy.random.rand(s_by_s.shape[0], s_by_s.shape[1])
        read_counts_gamma_all.append(read_counts_gamma)


    distances = list(coarse_grained_tree_dict_env.keys())
    distances.sort()

    distances = numpy.asarray(distances)

    dbd_dict = {}
    dbd_dict['phylo'] = {}
    for distance_idx in range(len(distances) - 1):
        
        fine_distance = distances[distance_idx]
        coarse_distance = distances[distance_idx + 1]

        sys.stderr.write("Phylo distance = %s \n" % round(coarse_distance, 7))

        coarse_grained_list = coarse_grained_tree_dict_env[coarse_distance]
        fine_grained_list = coarse_grained_tree_dict_env[fine_distance]

        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
        fine_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in fine_grained_list_i]) for fine_grained_list_i in fine_grained_list]

        # coarse grain s-by-s for all clades
        s_by_s_coarse = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
        s_by_s_fine = numpy.stack([numpy.sum(s_by_s[fine_grained_idx,:], axis=0) for fine_grained_idx in fine_grained_idx_all], axis=0)

        # generate fine and coarse matrices for null matrices
        s_by_s_coarse_null_all = []
        s_by_s_fine_null_all = []
        for i in range(iter_):
            read_counts_gamma_all_i = read_counts_gamma_all[i]
            s_by_s_coarse_null = numpy.stack([numpy.sum(read_counts_gamma_all_i[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
            s_by_s_fine_null = numpy.stack([numpy.sum(read_counts_gamma_all_i[fine_grained_idx,:], axis=0) for fine_grained_idx in fine_grained_idx_all], axis=0)

            s_by_s_coarse_null_all.append(s_by_s_coarse_null)
            s_by_s_fine_null_all.append(s_by_s_fine_null)


        slope_all = []
        slope_slm_all = []
        for focal_coarse_idx, focal_coarse in enumerate(coarse_grained_list):

            # ignore coarse-grained taxa with less than five fine-grained taxa
            if len(focal_coarse) < 5:
                continue

            fine_in_coarse_all = []
            # identify fine-grain taxa containing coarse-grained members
            f_idx_all = []
            for f_idx, f in enumerate(fine_grained_list):

                is_fine_in_coarse = bool(set(focal_coarse) & set(f))

                if is_fine_in_coarse == True:
                    fine_in_coarse_all.extend(f)
                    f_idx_all.append(f_idx)

            
            f_idx_all = numpy.asarray(f_idx_all)
            
            # ignore coarse-grained taxa with less than five fine-grained taxa
            if len(f_idx_all) < 5:
                continue
            
            # observed
            s_by_s_fine_focal = s_by_s_fine[f_idx_all,:]
            s_by_s_coarse_nonfocal = numpy.delete(s_by_s_coarse, focal_coarse_idx, axis=0)

            measure_coarse = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_coarse_nonfocal)
            measure_fine = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_fine_focal)

            if len(measure_fine) < 5:
                continue

            slope, intercept, r_value, p_value, std_err = stats.linregress(measure_coarse, measure_fine)

            # predicted
            slope_slm_iter_all = []
            for i in range(iter_):

                s_by_s_fine_null_i = s_by_s_fine_null_all[i]
                s_by_s_coarse_null_i = s_by_s_coarse_null_all[i]

                focal_fine_s_by_s_predicted = s_by_s_fine_null_i[f_idx_all,:]
                nonfocal_s_by_s_coarse_predicted = numpy.delete(s_by_s_coarse_null_i, focal_coarse_idx, axis=0)                

                focal_fine_diversity_predicted = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, focal_fine_s_by_s_predicted)
                nonfocal_coarse_diversity_predicted = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, nonfocal_s_by_s_coarse_predicted)

                slope_slm, intercept_slm, r_value_slm, p_value_slm, std_err_slm = stats.linregress(nonfocal_coarse_diversity_predicted, focal_fine_diversity_predicted)

                slope_slm_iter_all.append(slope_slm)


            slope_all.append(slope)
            slope_slm_all.append(numpy.mean(slope_slm_iter_all))


        slope_all = numpy.asarray(slope_all)
        slope_slm_all = numpy.asarray(slope_slm_all)

        idx_to_keep = ~(numpy.isnan(slope_all) | numpy.isnan(slope_slm_all))

        if sum(idx_to_keep) < 3:
            continue
        
        slope_all = slope_all[idx_to_keep]
        slope_slm_all = slope_slm_all[idx_to_keep]

        mean_slope_all = numpy.mean(slope_all)
        mean_slope_slm_all = numpy.mean(slope_slm_all)

        # error
        mean_error_slope_slm = numpy.mean(numpy.absolute( slope_slm_all - slope_all) / numpy.absolute(slope_all))

        dbd_dict['phylo'][coarse_distance] = {}
        dbd_dict['phylo'][coarse_distance]['slope_all'] = slope_all.tolist()
        dbd_dict['phylo'][coarse_distance]['slope_slm_all'] = slope_slm_all.tolist()

        dbd_dict['phylo'][coarse_distance]['mean_slope'] = mean_slope_all
        dbd_dict['phylo'][coarse_distance]['mean_slope_slm'] = mean_slope_slm_all

        dbd_dict['phylo'][coarse_distance]['mean_error_slope_slm'] = mean_error_slope_slm



    diversity_dbd_simulation_phylo_dict_path_ = diversity_dbd_simulation_phylo_dict_path % diversity_utils.get_environment_label(environment)

    sys.stderr.write("Saving diversity dictionary...\n")
    with open(diversity_dbd_simulation_phylo_dict_path_, 'wb') as handle:
        pickle.dump(dbd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)





def load_svd_and_make_diversity_slm_dbd_simulation_phylo_dict():

    dict_all = {}

    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        diversity_dbd_simulation_phylo_dict_path_ = diversity_dbd_simulation_phylo_dict_path % diversity_utils.get_environment_label(environment)

        if os.path.exists(diversity_dbd_simulation_phylo_dict_path_) == True:

            with open(diversity_dbd_simulation_phylo_dict_path_, 'rb') as handle:
                dict_ = pickle.load(handle)
            dict_all[environment] = dict_
    
    return dict_all




def richness_diversity_phylo_mean_simulation(environment, n_iter=100):

    '''
    Predict the mean and variance of diversity and richness using 
    the empirical mean and variance of realtive abundances at the 
    OTU level. Higher-level predictions generated by coarse-graining
    simulated data.

    contrasts with make_richness_diversity_prediction_phylo_dict() because
    here we coarse-grain fine-grain simulated data, instead of coarse-graining
    observed data and proceeding with simulations 
    '''

    environments_to_keep = diversity_utils.environments_to_keep

    coarse_grained_tree_dict_all = {}
    distances_all = []
    for e in environments_to_keep:
        sys.stderr.write("Loading tree dict for %s...\n" % e)
        coarse_grained_tree_dict = load_coarse_grained_tree_no_subsampling_dict(environment=e, rarefied=False)
        coarse_grained_tree_dict_all[e] = coarse_grained_tree_dict
        distances = list(coarse_grained_tree_dict.keys())
        distances_all.extend(distances)

    distances_all = list(set(distances_all))
    distances_all.sort()

    # phylogenetic coarse-graining
    sys.stderr.write("Subsetting samples for %s...\n" % environment)
    coarse_grained_tree_dict_env = coarse_grained_tree_dict_all[environment]
    sad_annotated_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
    samples = sad_annotated_dict['samples']

    taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

    # get SVD for simulations
    rel_s_by_s = s_by_s/numpy.sum(s_by_s, axis=0)
    n_reads = numpy.sum(s_by_s, axis=0)
    mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
    #beta_rel_s_by_s = (mean_rel_s_by_s**2)/var_rel_s_by_s

    slope, intercept, r_value, p_value, std_err = stats.linregress(numpy.log10(mean_rel_s_by_s), numpy.log10(var_rel_s_by_s))

    # taylors law
    fixed_null_mean = numpy.mean(mean_rel_s_by_s)
    fixed_null_var = (10**( numpy.log10(fixed_null_mean)*slope + intercept))
    var_taylors = (10**(numpy.log10(mean_rel_s_by_s)*slope + intercept))
    mean_rel_s_by_s_fixed = numpy.asarray([fixed_null_mean]*len(mean_rel_s_by_s))
    var_rel_s_by_s_fixed = numpy.asarray([fixed_null_var]*len(mean_rel_s_by_s))

    #svd = load_svd_phylo_dict(environment)

    u_null = numpy.identity(len(mean_rel_s_by_s))
    s_null = numpy.asarray([1]*len(mean_rel_s_by_s))
    v_null = numpy.identity(len(mean_rel_s_by_s))
    svd_null = (u_null, s_null, v_null)

    measure_dict = {}
    coarse_grain_index_dict = {}
    
    distances_all = [distances_all[10]]
    #distances_all = distances_all[:20]

    for distance_idx, distance in enumerate(distances_all):

        if distance not in coarse_grained_tree_dict:
            continue

        measure_dict[distance] = {}
        measure_dict[distance]['fixed_mean_phylo'] = {}
        measure_dict[distance]['fixed_mean_phylo']['mean_richness'] = []
        measure_dict[distance]['fixed_mean_phylo']['mean_diversity'] = []
        measure_dict[distance]['fixed_mean_phylo']['var_richness'] = []
        measure_dict[distance]['fixed_mean_phylo']['var_diversity'] = []

        # observed MAD, variance from taylor's law
        measure_dict[distance]['observed_mean_phylo'] = {}
        measure_dict[distance]['observed_mean_phylo']['mean_richness'] = []
        measure_dict[distance]['observed_mean_phylo']['mean_diversity'] = []
        measure_dict[distance]['observed_mean_phylo']['var_richness'] = []
        measure_dict[distance]['observed_mean_phylo']['var_diversity'] = []

        # observed MAD, variance from taylor's law, null coarse_graining
        measure_dict[distance]['observed_mean_null'] = {}
        measure_dict[distance]['observed_mean_null']['mean_richness'] = []
        measure_dict[distance]['observed_mean_null']['mean_diversity'] = []
        measure_dict[distance]['observed_mean_null']['var_richness'] = []
        measure_dict[distance]['observed_mean_null']['var_diversity'] = []

        # index dict
        #coarse_grained_list = coarse_grained_tree_dict[distance]
        coarse_grained_list = coarse_grained_tree_dict_env[distance]
        coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])
        #coarse_grain_idx = numpy.append([0], numpy.cumsum(coarse_grained_n))[:-1]

        # get indexes for each clade
        coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

        coarse_grain_index_dict[distance] = {}
        coarse_grain_index_dict[distance]['coarse_grained_n'] = coarse_grained_n
        coarse_grain_index_dict[distance]['coarse_grained_idx_all'] = coarse_grained_idx_all

        print(distance)


    sys.stderr.write("Simulating abundance matrices...\n")

    #mean_richness_observed_gamma_rvs, var_richness_observed_gamma_rvs, mean_richness_predicted_gamma_rvs, var_richness_predicted_gamma_rvs, mean_diversity_observed_gamma_rvs, var_diversity_observed_gamma_rvs, mean_diversity_predicted_gamma_rvs, var_diversity_predicted_gamma_rvs = diversity_utils.predict_mean_and_var_richness_and_diversity_using_gamma_rv(s_by_s, iter_=100)
    #print(mean_diversity_observed_gamma_rvs, mean_diversity_predicted_gamma_rvs)

    distances_to_keep = list(measure_dict.keys())
    distances_to_keep.sort()
    for i in range(n_iter):

        print(i)

        #rel_abundances_null_fixed_mean_i, read_counts_multinomial_null_fixed_mean_i = simulation_utils.genrate_community_from_mean_and_var_svd(mean_rel_s_by_s_fixed, var_rel_s_by_s_fixed, n_reads, len(n_reads), svd_null)
        rel_abundances_null_fixed_mean_i, read_counts_multinomial_null_fixed_mean_i = simulation_utils.genrate_community_from_mean_and_var_using_rvs(mean_rel_s_by_s_fixed, var_rel_s_by_s_fixed, n_reads, len(n_reads))

        #rel_abundances_null_i, read_counts_multinomial_null_i = simulation_utils.genrate_community_from_mean_and_var_svd(mean_rel_s_by_s, var_taylors, n_reads, len(n_reads), svd_null)
        #rel_abundances_null_i, read_counts_multinomial_null_i = simulation_utils.genrate_community_from_mean_and_var_using_rvs(mean_rel_s_by_s, var_taylors, n_reads, len(n_reads))
        rel_abundances_null_i, read_counts_multinomial_null_i = simulation_utils.genrate_community_from_mean_and_var_using_rvs(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads))

        for distance in distances_to_keep:

            coarse_grained_n = coarse_grain_index_dict[distance]['coarse_grained_n']
            coarse_grained_idx_all = coarse_grain_index_dict[distance]['coarse_grained_idx_all']

            s_by_s_null_fixed_mean_coarse = numpy.stack([numpy.sum(read_counts_multinomial_null_fixed_mean_i[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
            
            # phylogenetic coarse-graining
            s_by_s_null_coarse_phylo = numpy.stack([numpy.sum(read_counts_multinomial_null_i[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
            
            # random coarse-graining
            numpy.random.shuffle(coarse_grained_n)
            coarse_grain_null_idx = numpy.append([0], numpy.cumsum(coarse_grained_n))[:-1]
            s_by_s_null_coarse_null = numpy.add.reduceat(read_counts_multinomial_null_i, coarse_grain_null_idx, axis=0)

            # get measures
            richness_null_fixed_mean_i = numpy.apply_along_axis(diversity_utils.calculate_richness, 0, s_by_s_null_fixed_mean_coarse)
            diversity_null_fixed_mean_i = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_null_fixed_mean_coarse)

            richness_null_phylo_i = numpy.apply_along_axis(diversity_utils.calculate_richness, 0, s_by_s_null_coarse_phylo)
            diversity_null_phylo_i = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_null_coarse_phylo)

            richness_null_null_i = numpy.apply_along_axis(diversity_utils.calculate_richness, 0, s_by_s_null_coarse_null)
            diversity_null_null_i = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_null_coarse_null)


            measure_dict[distance]['fixed_mean_phylo']['mean_richness'].append(numpy.mean(richness_null_fixed_mean_i))
            measure_dict[distance]['fixed_mean_phylo']['mean_diversity'].append(numpy.mean(diversity_null_fixed_mean_i))
            measure_dict[distance]['fixed_mean_phylo']['var_richness'].append(numpy.var(richness_null_fixed_mean_i))
            measure_dict[distance]['fixed_mean_phylo']['var_diversity'].append(numpy.var(diversity_null_fixed_mean_i))

            # observed MAD, variance from taylor's law
            measure_dict[distance]['observed_mean_phylo']['mean_richness'].append(numpy.mean(richness_null_phylo_i))
            measure_dict[distance]['observed_mean_phylo']['mean_diversity'].append(numpy.mean(diversity_null_phylo_i))
            measure_dict[distance]['observed_mean_phylo']['var_richness'].append(numpy.var(richness_null_phylo_i))
            measure_dict[distance]['observed_mean_phylo']['var_diversity'].append(numpy.var(diversity_null_phylo_i))

            # observed MAD, variance from taylor's law, null coarse_graining
            measure_dict[distance]['observed_mean_null']['mean_richness'].append(numpy.mean(richness_null_null_i))
            measure_dict[distance]['observed_mean_null']['mean_diversity'].append(numpy.mean(diversity_null_null_i))
            measure_dict[distance]['observed_mean_null']['var_richness'].append(numpy.var(richness_null_null_i))
            measure_dict[distance]['observed_mean_null']['var_diversity'].append(numpy.var(diversity_null_null_i))



    dbd_dict = {}
    dbd_dict['phylo'] = {}

    for distance in distances_to_keep:

        dbd_dict['phylo'][distance] = {}
        dbd_dict['phylo'][distance]['observed'] = {}
        dbd_dict['phylo'][distance]['simulated'] = {}

        coarse_grained_idx_all = coarse_grain_index_dict[distance]['coarse_grained_idx_all']

        # observed data
        # coarse grain s-by-s for all clades
        s_by_s_coarse = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
        richness_observed = numpy.apply_along_axis(diversity_utils.calculate_richness, 0, s_by_s_coarse)
        diversity_observed = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_coarse)

        dbd_dict['phylo'][distance]['observed']['mean_richness'] = numpy.mean(richness_observed)
        dbd_dict['phylo'][distance]['observed']['mean_diversity'] = numpy.mean(diversity_observed)
        dbd_dict['phylo'][distance]['observed']['var_richness'] = numpy.var(richness_observed)
        dbd_dict['phylo'][distance]['observed']['var_diversity'] = numpy.var(diversity_observed)

        # fixed mean
        dbd_dict['phylo'][distance]['simulated']['fixed_mean_phylo'] = {}
        dbd_dict['phylo'][distance]['simulated']['fixed_mean_phylo']['mean_richness'] = numpy.mean(measure_dict[distance]['fixed_mean_phylo']['mean_richness'])
        dbd_dict['phylo'][distance]['simulated']['fixed_mean_phylo']['mean_diversity'] = numpy.mean(measure_dict[distance]['fixed_mean_phylo']['mean_diversity'])
        dbd_dict['phylo'][distance]['simulated']['fixed_mean_phylo']['var_richness'] = numpy.var(measure_dict[distance]['fixed_mean_phylo']['var_richness'])
        dbd_dict['phylo'][distance]['simulated']['fixed_mean_phylo']['var_diversity'] = numpy.var(measure_dict[distance]['fixed_mean_phylo']['var_diversity'])

        # observed MAD, variance from taylor's law
        dbd_dict['phylo'][distance]['simulated']['observed_mean_phylo'] = {}
        dbd_dict['phylo'][distance]['simulated']['observed_mean_phylo']['mean_richness'] = numpy.mean(measure_dict[distance]['observed_mean_phylo']['mean_richness'])
        dbd_dict['phylo'][distance]['simulated']['observed_mean_phylo']['mean_diversity'] = numpy.mean(measure_dict[distance]['observed_mean_phylo']['mean_diversity'])
        dbd_dict['phylo'][distance]['simulated']['observed_mean_phylo']['var_richness'] = numpy.var(measure_dict[distance]['observed_mean_phylo']['var_richness'])
        dbd_dict['phylo'][distance]['simulated']['observed_mean_phylo']['var_diversity'] = numpy.var(measure_dict[distance]['observed_mean_phylo']['var_diversity'])

        # observed MAD, variance from taylor's law, null coarse_graining
        dbd_dict['phylo'][distance]['simulated']['observed_mean_null'] = {}
        dbd_dict['phylo'][distance]['simulated']['observed_mean_null']['mean_richness'] = numpy.mean(measure_dict[distance]['observed_mean_null']['mean_richness'])
        dbd_dict['phylo'][distance]['simulated']['observed_mean_null']['mean_diversity'] = numpy.mean(measure_dict[distance]['observed_mean_null']['mean_diversity'])
        dbd_dict['phylo'][distance]['simulated']['observed_mean_null']['var_richness'] = numpy.var(measure_dict[distance]['observed_mean_null']['var_richness'])
        dbd_dict['phylo'][distance]['simulated']['observed_mean_null']['var_diversity'] = numpy.var(measure_dict[distance]['observed_mean_null']['var_diversity'])


        print(distance, numpy.mean(measure_dict[distance]['observed_mean_phylo']['mean_diversity']),  numpy.mean(diversity_observed))


        # later, add correlation...


    diversity_dbd_simulation_phylo_dict_path_ = richness_diversity_phylo_mean_simulation_path % diversity_utils.get_environment_label(environment)
    sys.stderr.write("Saving diversity dictionary...\n")
    with open(diversity_dbd_simulation_phylo_dict_path_, 'wb') as handle:
        pickle.dump(dbd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)





    
