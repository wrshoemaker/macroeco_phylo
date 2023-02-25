
import os
from Bio import Phylo
import random
import copy
import sys

import treeswift

random.seed(123456789)

directory = os.path.expanduser("~/GitHub/macroeco_phylo/")



#tree = Phylo.read(directory +'97_otus.tre', "newick")

tree = Phylo.read(directory +'test.tre', "newick")


term_names = [term.name for term in tree.get_terminals()]

term_names_sample = random.sample(term_names, 10)
term_names_sample = ['L131', 'L169', 'L170', 'L167', 'L168']

print(term_names_sample)
treecopy = tree

#treecopy = copy.deepcopy(tree)



count = 0

for leaf in treecopy.get_terminals():

    count +=1

    if count % 10000 == 0:
        sys.stderr.write("%dk tips completed...\n" % (count/1000))

    if leaf not in term_names_sample:

        try:
            # preserve_branch_length=True
            treecopy.prune(leaf)
        except ValueError as err:
            continue
    else:
        continue
#



term_names = [term.name for term in treecopy.get_terminals()]
print(term_names)




def make_richness_diversity_prediction_phylo_dict_old():

    environments_to_keep = diversity_utils.environments_to_keep

    coarse_grained_tree_dict_all = {}
    distances_all = []
    for environment in environments_to_keep:
        sys.stderr.write("Loading tree dict for %s...\n" % environment)
        coarse_grained_tree_dict = load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=False)
        coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict
        distances = list(coarse_grained_tree_dict.keys())
        distances_all.extend(distances)

    distances_all = list(set(distances_all))
    distances_all.sort()

    richness_diversity_prediction_dict = {}    

    for environment_idx, environment in enumerate(environments_to_keep):

        coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]

        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        pres_abs_dict = load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
        samples = pres_abs_dict['samples']
        taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
       
        sys.stderr.write("Getting site-by-species matrix...\n")
        s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])


        richness_diversity_prediction_dict[environment] = {}
        richness_diversity_prediction_dict[environment]['phylo'] = {}

        for distance_idx, distance in enumerate(distances):

            sys.stderr.write("Phylo distance = %s \n" % round(distance, 7))
            coarse_grained_list = coarse_grained_tree_dict[distance]
            coarse_grained_n = numpy.asarray([len(i) for i in coarse_grained_list])

            # get indexes for each clade
            coarse_grained_idx_all = [numpy.asarray([numpy.where(taxa==i)[0][0] for i in coarse_grained_list_i]) for coarse_grained_list_i in coarse_grained_list]

            # coarse grain s-by-s for all clades
            s_by_s_all_clades = numpy.stack([numpy.sum(s_by_s[coarse_grained_idx,:], axis=0) for coarse_grained_idx in coarse_grained_idx_all], axis=0)
            
            mean_richness_observed_gamma_rvs, var_richness_observed_gamma_rvs, mean_richness_predicted_gamma_rvs, var_richness_predicted_gamma_rvs, mean_richness_predicted_gamma_rvs_corr, var_richness_predicted_gamma_rvs_corr, mean_diversity_observed_gamma_rvs, var_diversity_observed_gamma_rvs, mean_diversity_predicted_gamma_rvs, var_diversity_predicted_gamma_rvs, mean_diversity_predicted_gamma_rvs_corr, var_diversity_predicted_gamma_rvs_corr = diversity_utils.predict_mean_and_var_richness_and_diversity_using_gamma_rv(s_by_s_all_clades, iter=10)
            mean_diversity_observed, mean_diversity_predicted, var_diversity_observed, var_diversity_predicted = diversity_utils.predict_mean_and_var_diversity_analytic(s_by_s_all_clades, s_by_s_all_clades.shape[0])
            
            mean_richness_observed, mean_richness_predicted = diversity_utils.predict_mean_richness(s_by_s_all_clades, s_by_s_all_clades.shape[0])
            var_richness_observed, var_richness_predicted = diversity_utils.predict_var_richness(s_by_s_all_clades, s_by_s_all_clades.shape[0])


            richness_diversity_prediction_dict[environment]['phylo'][distance] = {}

            # mean richness
            richness_diversity_prediction_dict[environment]['phylo'][distance]['mean_richness_observed'] = mean_richness_observed
            richness_diversity_prediction_dict[environment]['phylo'][distance]['mean_richness_predicted'] = mean_richness_predicted
            richness_diversity_prediction_dict[environment]['phylo'][distance]['mean_richness_predicted_gamma'] = mean_richness_predicted_gamma_rvs
            richness_diversity_prediction_dict[environment]['phylo'][distance]['mean_richness_predicted_gamma_corr'] = mean_richness_predicted_gamma_rvs_corr

            # var richness
            richness_diversity_prediction_dict[environment]['phylo'][distance]['var_richness_observed'] = var_richness_observed
            richness_diversity_prediction_dict[environment]['phylo'][distance]['var_richness_predicted'] = var_richness_predicted
            richness_diversity_prediction_dict[environment]['phylo'][distance]['var_richness_predicted_gamma'] = var_richness_predicted_gamma_rvs
            richness_diversity_prediction_dict[environment]['phylo'][distance]['var_richness_predicted_gamma_corr'] = var_richness_predicted_gamma_rvs_corr

            # mean diversity
            richness_diversity_prediction_dict[environment]['phylo'][distance]['mean_diversity_observed'] = mean_diversity_observed
            richness_diversity_prediction_dict[environment]['phylo'][distance]['mean_diversity_predicted'] = mean_diversity_predicted
            richness_diversity_prediction_dict[environment]['phylo'][distance]['mean_diversity_predicted_gamma'] = mean_diversity_predicted_gamma_rvs
            richness_diversity_prediction_dict[environment]['phylo'][distance]['mean_diversity_predicted_gamma_corr'] = mean_diversity_predicted_gamma_rvs_corr

            # var diversity
            richness_diversity_prediction_dict[environment]['phylo'][distance]['var_diversity_observed'] = var_diversity_observed
            richness_diversity_prediction_dict[environment]['phylo'][distance]['var_diversity_predicted'] = var_diversity_predicted
            richness_diversity_prediction_dict[environment]['phylo'][distance]['var_diversity_predicted_gamma'] = var_diversity_predicted_gamma_rvs
            richness_diversity_prediction_dict[environment]['phylo'][distance]['var_diversity_predicted_gamma_corr'] = var_diversity_predicted_gamma_rvs_corr



    sys.stderr.write("Saving prediction dictionary...\n")
    with open(richness_diversity_prediction_phylo_dict_path, 'wb') as handle:
        pickle.dump(richness_diversity_prediction_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)





#subset_tree = tree.prune(term_names_sample, preserve_branch_length=True)
