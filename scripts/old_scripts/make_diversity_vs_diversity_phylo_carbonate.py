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

rarefied = True
iter = 1000


coarse_grained_tree_path_template = config.data_directory + "coarse_grained_tree/%s.pickle"
sad_path_template = config.data_directory + "sad_annotated/sad_annotated_%s.pickle"
#diversity_vs_diversity_phylo_template = config.data_directory + "diversity_vs_diversity_phylo/diversity_vs_diversity_phylo%s.pickle"

diversity_vs_diversity_phylo_emp_template = config.data_directory + "diversity_vs_diversity_phylo_emp/diversity_vs_diversity_phylo%s.pickle"



def load_coarse_grained_tree_dict(environment='human gut metagenome', rarefied=False):

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    coarse_grained_tree_path = coarse_grained_tree_path_template % (environment.replace(' ', '_') + rarefied_label)

    with open(coarse_grained_tree_path, 'rb') as handle:
        distance_collapsed_dict = pickle.load(handle)
    return distance_collapsed_dict



def load_sad_annotated_taxon_dict(environment, rarefied=False):

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    sad_path = sad_path_template % (environment.replace(' ', '_') + rarefied_label)

    with open(sad_path, 'rb') as handle:
        presence_absence_dict = pickle.load(handle)
    return presence_absence_dict




diversity_vs_diversity_phylo_dict = {}
sys.stderr.write("Loading tree dict for %s...\n" % environment)

coarse_grained_tree_dict = load_coarse_grained_tree_dict(environment=environment, rarefied=rarefied)

distances = list(coarse_grained_tree_dict.keys())
distances.sort()

sys.stderr.write("Subsetting samples...\n")
#samples = diversity_utils.subset_observations(environment=environment)
pres_abs_dict = load_sad_annotated_taxon_dict(environment, rarefied = rarefied)
samples = pres_abs_dict['samples']

sys.stderr.write("Getting site-by-species matrix...\n")
s_by_s, taxonomy_names, samples_keep = diversity_utils.get_s_by_s(samples)
rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

diversity_species = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s)

diversity_vs_diversity_phylo_dict['diversity_species'] = diversity_species.tolist()

for distance in distances:

    sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
    coarse_grained_list = coarse_grained_tree_dict[distance]
    coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

    #print(len(coarse_grained_n), len(taxonomy_names))

    if len(coarse_grained_n) == len(taxonomy_names):
        continue

    # get indexes for each clade
    coarse_grained_idx_all = [numpy.asarray([numpy.where(taxonomy_names==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]
    # coarse grain s-by-s for all clades
    s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)

    diversity_coarse = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_all_clades)
    slope, intercept, r_value, p_value, std_err = stats.linregress(diversity_species, diversity_coarse)


    # AFD
    rel_s_by_s_all_clades = (s_by_s_all_clades/s_by_s_all_clades.sum(axis=0))
    clade_log10_rescaled_all = []
    for clade in rel_s_by_s_all_clades:

        clade = clade[clade>0]
        clade_log10 = numpy.log10(clade)

        if len(clade_log10) < 4:
            continue

        clade_log10_rescaled = (clade_log10 - numpy.mean(clade_log10))/numpy.std(clade_log10)
        clade_log10_rescaled_all.append(clade_log10_rescaled)


    # taylors law
    mean_rel_s_by_s_all_clades = numpy.mean(rel_s_by_s_all_clades, axis=1)
    var_rel_s_by_s_all_clades = numpy.var(rel_s_by_s_all_clades, axis=1)
    taylors_law_slope, taylors_law_intercept, taylors_law_r_value, taylors_law_p_value, taylors_law_std_err = stats.linregress(numpy.log10(mean_rel_s_by_s_all_clades), numpy.log10(var_rel_s_by_s_all_clades))


    # MAD
    mean_rel_s_by_s_all_clades_log10 = numpy.log10(mean_rel_s_by_s_all_clades)
    mad_log10_rescaled = (mean_rel_s_by_s_all_clades_log10 - numpy.mean(mean_rel_s_by_s_all_clades_log10)) / numpy.std(mean_rel_s_by_s_all_clades_log10)


    # max rel abundancs fine grain vs rel abundance coarse-grained
    #max_mean_fine_all = numpy.asarray([max(numpy.mean(rel_s_by_s[coarse_grained_idx,:], axis=1)) for coarse_grained_idx in coarse_grained_idx_all])
    #mean_coarse = numpy.mean([numpy.sum(rel_s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=1)
    #mean_ratio_all = max_mean_fine_all/mean_coarse


    # gamma error
    #occupancies, predicted_occupancies, mad, beta, species = diversity_utils.predict_occupancy(s_by_s_all_clades, list(range(len(coarse_grained_idx_all))))
    #error = numpy.absolute(occupancies - predicted_occupancies)/occupancies
    #dominance = numpy.asarray([numpy.mean(numpy.amax(rel_s_by_s[coarse_grained_idx,:], axis=1)) for coarse_grained_idx in coarse_grained_idx_all])
    #error_idx = (error>0)
    #error = error[error_idx]
    #dominance = dominance[error_idx]
    #dominance_vs_error_rho = numpy.corrcoef(numpy.log10(dominance), numpy.log10(error))[1,0]

    #print(dominance_vs_error_rho)

    coarse_grained_n = numpy.asarray([len(coarse_grained_idx) for coarse_grained_idx in coarse_grained_idx_all])
    idx_all = numpy.arange(len(coarse_grained_idx_all))
    coarse_grain_idx = numpy.append([0], numpy.cumsum(coarse_grained_n))[:-1]
    s_by_s_copy = numpy.copy(s_by_s)
    slope_null_all = []
    clade_log10_rescaled_null_all = []
    taylors_law_slope_null_all = []
    taylors_law_intercept_null_all = []
    mad_log10_rescaled_null_all = []
    mean_ratio_all_null_all = []
    dominance_vs_error_rho_null_all = []
    for i in range(iter):

        # diversity vs diversity null
        numpy.random.shuffle(s_by_s_copy)
        s_by_s_coarse_null = numpy.add.reduceat(s_by_s_copy, coarse_grain_idx, axis=0)
        s_by_s_coarse_null = s_by_s_coarse_null[:,~(numpy.all(s_by_s_coarse_null == 0, axis=0))]
        diversity_coarse_null = numpy.apply_along_axis(diversity_utils.calculate_shannon_diversity, 0, s_by_s_coarse_null)
        slope_null, intercept_null, r_value_null, p_value_null, std_err_null = stats.linregress(diversity_species, diversity_coarse_null)
        slope_null_all.append(slope_null)

        # AFD null
        rel_s_by_s_coarse_null = (s_by_s_coarse_null/s_by_s_coarse_null.sum(axis=0))
        for clade in rel_s_by_s_coarse_null:

            clade = clade[clade>0]
            clade_log10 = numpy.log10(clade)

            if len(clade_log10) < 4:
                continue

            clade_log10_rescaled = (clade_log10 - numpy.mean(clade_log10))/numpy.std(clade_log10)
            clade_log10_rescaled_null_all.append(clade_log10_rescaled)


        # taylors law null
        mean_rel_s_by_s_coarse_null = numpy.mean(rel_s_by_s_coarse_null, axis=1)
        var_rel_s_by_s_coarse_null = numpy.var(rel_s_by_s_coarse_null, axis=1)
        taylors_law_slope_null, taylors_law_intercept_null, taylors_law_r_value_null, taylors_law_p_value_null, taylors_law_std_err_null = stats.linregress(numpy.log10(mean_rel_s_by_s_coarse_null), numpy.log10(var_rel_s_by_s_coarse_null))
        taylors_law_slope_null_all.append(taylors_law_slope_null)
        taylors_law_intercept_null_all.append(taylors_law_intercept_null)


        # MAD null
        mean_rel_s_by_s_coarse_null_log10 = numpy.log10(mean_rel_s_by_s_coarse_null)
        mad_log10_rescaled_null = (mean_rel_s_by_s_coarse_null_log10 - numpy.mean(mean_rel_s_by_s_coarse_null_log10)) / numpy.std(mean_rel_s_by_s_coarse_null_log10)
        mad_log10_rescaled_null_all.append(mad_log10_rescaled_null)


        # max rel abundancs fine grain vs rel abundance coarse-grained null
        #rel_s_by_s_copy = (s_by_s_copy/s_by_s_copy.sum(axis=0))
        #max_mean_fine_all_null = numpy.asarray([max(numpy.mean(rel_s_by_s_copy[coarse_grained_idx,:], axis=1)) for coarse_grained_idx in coarse_grained_idx_all])
        #mean_coarse_null = numpy.mean([numpy.sum(rel_s_by_s_copy[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=1)
        #mean_ratio_all_null = max_mean_fine_all_null/mean_coarse_null
        #mean_ratio_all_null_all.append(mean_ratio_all_null)


        #numpy.add.reduceat(s_by_s_copy, coarse_grain_idx, axis=0)

        #print('start gamma error')

        # gamma error null
        #occupancies_null, predicted_occupancies_null, mad_null, beta_null, species_null = diversity_utils.predict_occupancy(s_by_s_coarse_null, list(range(len(coarse_grained_idx_all))))
        #error_null = numpy.absolute(occupancies_null - predicted_occupancies_null)/occupancies_null
        #dominance_null = numpy.asarray([numpy.mean(numpy.amax(rel_s_by_s_copy[coarse_grained_idx,:], axis=1)) for coarse_grained_idx in coarse_grained_idx_all])
        #error_null_idx = (error_null>0)
        #error_null = error_null[error_null_idx]
        #dominance_null = dominance_null[error_null_idx]
        #dominance_vs_error_rho_null = numpy.corrcoef(numpy.log10(dominance_null), numpy.log10(error_null))[1,0]
        #dominance_vs_error_rho_null_all.append(dominance_vs_error_rho_null)



    #diversity_vs_diversity
    slope_null_all = numpy.asarray(slope_null_all)
    slope_null_all = numpy.sort(slope_null_all)

    percentile = sum(slope_null_all < slope)/iter

    lower_ci = slope_null_all[int(iter*0.025)]
    upper_ci = slope_null_all[int(iter*0.975)]

    diversity_vs_diversity_phylo_dict[distance] = {}
    diversity_vs_diversity_phylo_dict[distance]['diversity_vs_diversity'] = {}
    diversity_vs_diversity_phylo_dict[distance]['diversity_vs_diversity']['diversity'] = diversity_coarse.tolist()
    diversity_vs_diversity_phylo_dict[distance]['diversity_vs_diversity']['slope'] = slope
    diversity_vs_diversity_phylo_dict[distance]['diversity_vs_diversity']['percentile'] = percentile
    diversity_vs_diversity_phylo_dict[distance]['diversity_vs_diversity']['lower_ci'] = lower_ci
    diversity_vs_diversity_phylo_dict[distance]['diversity_vs_diversity']['upper_ci'] = upper_ci


    # AFD
    def get_hist_and_bins(flat_array):

        # null is too large, so we are binning it for the plot in this script
        hist_, bin_edges_ = numpy.histogram(flat_array, density=True, bins=20)
        bins_mean_ = numpy.asarray([0.5 * (bin_edges_[i] + bin_edges_[i+1]) for i in range(0, len(bin_edges_)-1 )])
        hist_to_plot = hist_[hist_>0]
        bins_mean_to_plot = bins_mean_[hist_>0]

        return hist_to_plot, bins_mean_to_plot


    clade_log10_rescaled_all_flat = numpy.concatenate(clade_log10_rescaled_all).ravel()
    clade_log10_rescaled_null_all_flat = numpy.concatenate(clade_log10_rescaled_null_all).ravel()
    hist_to_plot, bins_mean_to_plot = get_hist_and_bins(clade_log10_rescaled_all_flat)
    hist_to_plot_null, bins_mean_to_plot_null = get_hist_and_bins(clade_log10_rescaled_null_all_flat)

    diversity_vs_diversity_phylo_dict[distance]['afd'] = {}
    diversity_vs_diversity_phylo_dict[distance]['afd']['observed_hist_to_plot'] = hist_to_plot.tolist()
    diversity_vs_diversity_phylo_dict[distance]['afd']['observed_bins_mean_to_plot'] = bins_mean_to_plot.tolist()
    diversity_vs_diversity_phylo_dict[distance]['afd']['null_hist_to_plot'] = hist_to_plot_null.tolist()
    diversity_vs_diversity_phylo_dict[distance]['afd']['null_bins_mean_to_plot'] = bins_mean_to_plot_null.tolist()


    # taylors law
    taylors_law_slope_null_all = numpy.asarray(taylors_law_slope_null_all)
    taylors_law_slope_null_all = numpy.sort(taylors_law_slope_null_all)
    taylors_law_intercept_null_all = numpy.asarray(taylors_law_intercept_null_all)
    taylors_law_intercept_null_all = numpy.sort(taylors_law_intercept_null_all)

    #percentile = sum(slope_null_all < slope)/iter
    taylors_law_slope_null_lower_ci = taylors_law_slope_null_all[int(iter*0.025)]
    taylors_law_slope_null_upper_ci = taylors_law_slope_null_all[int(iter*0.975)]
    taylors_law_intercept_null_lower_ci = taylors_law_intercept_null_all[int(iter*0.025)]
    taylors_law_intercept_null_upper_ci = taylors_law_intercept_null_all[int(iter*0.975)]

    diversity_vs_diversity_phylo_dict[distance]['taylors_law'] = {}
    diversity_vs_diversity_phylo_dict[distance]['taylors_law']['slope'] = taylors_law_slope
    diversity_vs_diversity_phylo_dict[distance]['taylors_law']['slope_lower_ci'] = taylors_law_slope_null_lower_ci
    diversity_vs_diversity_phylo_dict[distance]['taylors_law']['slope_upper_ci'] = taylors_law_slope_null_upper_ci

    diversity_vs_diversity_phylo_dict[distance]['taylors_law']['intercept'] = taylors_law_intercept
    diversity_vs_diversity_phylo_dict[distance]['taylors_law']['intercept_lower_ci'] = taylors_law_intercept_null_lower_ci
    diversity_vs_diversity_phylo_dict[distance]['taylors_law']['intercept_upper_ci'] = taylors_law_intercept_null_upper_ci


    # MAD
    mad_log10_rescaled_null_all_flat = numpy.concatenate(mad_log10_rescaled_null_all).ravel()
    hist_to_plot, bins_mean_to_plot = get_hist_and_bins(mad_log10_rescaled)
    hist_to_plot_null, bins_mean_to_plot_null = get_hist_and_bins(mad_log10_rescaled_null_all_flat)

    diversity_vs_diversity_phylo_dict[distance]['mad'] = {}
    diversity_vs_diversity_phylo_dict[distance]['mad']['observed_hist_to_plot'] = hist_to_plot.tolist()
    diversity_vs_diversity_phylo_dict[distance]['mad']['observed_bins_mean_to_plot'] = bins_mean_to_plot.tolist()
    diversity_vs_diversity_phylo_dict[distance]['mad']['null_hist_to_plot'] = hist_to_plot_null.tolist()
    diversity_vs_diversity_phylo_dict[distance]['mad']['null_bins_mean_to_plot'] = bins_mean_to_plot_null.tolist()


    # max relative abundance vs. coarse grained abundance
    #mean_ratio_all_null_all = numpy.asarray(mean_ratio_all_null_all)
    #mean_mean_ratio_all_null = []
    #lower_ci_mean_ratio_all_null = []
    #upper_ci_mean_ratio_all_null = []
    #for null_array in mean_ratio_all_null_all.T:
    #    null_array = numpy.sort(null_array)
    #    mean_ratio_null = numpy.mean(null_array)
    #    lower_ci_null_array = null_array[int(iter*0.025)]
    #    upper_ci_null_array = null_array[int(iter*0.975)]

    #    mean_mean_ratio_all_null.append(mean_ratio_null)
    #    lower_ci_mean_ratio_all_null.append(lower_ci_null_array)
    #    upper_ci_mean_ratio_all_null.append(upper_ci_null_array)

    #diversity_vs_diversity_phylo_dict[distance]['max_mean_fine_vs_mean_coarse'] = {}
    #diversity_vs_diversity_phylo_dict[distance]['max_mean_fine_vs_mean_coarse']['observed'] = mean_ratio_all.tolist()
    #diversity_vs_diversity_phylo_dict[distance]['max_mean_fine_vs_mean_coarse']['null_mean'] = mean_mean_ratio_all_null
    #diversity_vs_diversity_phylo_dict[distance]['max_mean_fine_vs_mean_coarse']['null_lower_ci'] = lower_ci_mean_ratio_all_null
    #diversity_vs_diversity_phylo_dict[distance]['max_mean_fine_vs_mean_coarse']['null_upper_ci'] = upper_ci_mean_ratio_all_null


    # gamma error
    #dominance_vs_error_rho_null_all = numpy.asarray(dominance_vs_error_rho_null_all)
    #dominance_vs_error_rho_null_all = numpy.sort(dominance_vs_error_rho_null_all)

    #lower_ci_dominance_vs_error_rho_null = dominance_vs_error_rho_null_all[int(iter*0.025)]
    #upper_ci_dominance_vs_error_rho_null = dominance_vs_error_rho_null_all[int(iter*0.975)]

    #diversity_vs_diversity_phylo_dict[distance]['dominance_vs_error'] = {}
    #diversity_vs_diversity_phylo_dict[distance]['dominance_vs_error']['dominance'] = dominance.tolist()
    #diversity_vs_diversity_phylo_dict[distance]['dominance_vs_error']['error'] = error.tolist()
    #diversity_vs_diversity_phylo_dict[distance]['dominance_vs_error']['rho'] = dominance_vs_error_rho
    #diversity_vs_diversity_phylo_dict[distance]['dominance_vs_error']['rho_null_lower_ci'] = lower_ci_dominance_vs_error_rho_null
    #diversity_vs_diversity_phylo_dict[distance]['dominance_vs_error']['rho_null_upper_ci'] = upper_ci_dominance_vs_error_rho_null




if rarefied == True:
    rarefied_label = '_rarefied'
else:
    rarefied_label = ''

environment_label = environment.replace(' ', '_')

diversity_vs_diversity_phylo_path = diversity_vs_diversity_phylo_emp_template % ('_' + environment_label + rarefied_label)

sys.stderr.write("Saving diversity vs. diversity dictionary...\n")
with open(diversity_vs_diversity_phylo_path, 'wb') as handle:
    pickle.dump(diversity_vs_diversity_phylo_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
