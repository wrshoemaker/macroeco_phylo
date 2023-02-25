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
import scipy.stats as stats

import ete3
import tree_utils
import random
import dbd_utils
import plot_utils

import matplotlib.pyplot as plt

rarefied = False



coarse_grained_tree_dict = dbd_utils.get_coarse_grained_tree_dict()

richness_dbd_dict_path = config.data_directory + "richness_dbd_dict.pickle"



def predict_richness_dbd(focal_mean_fine_rel_s_by_s, focal_var_fine_rel_s_by_s, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, nonfocal_var_coarse_rel_s_by_s, total_reads_coarse):

    focal_fine_beta = (focal_mean_fine_rel_s_by_s**2)/focal_var_fine_rel_s_by_s
    nonfocal_coarse_beta = (nonfocal_mean_coarse_rel_s_by_s**2)/nonfocal_var_coarse_rel_s_by_s
   
    focal_fine_theta = focal_mean_fine_rel_s_by_s/focal_fine_beta
    nonfocal_coarse_theta = nonfocal_mean_coarse_rel_s_by_s/nonfocal_coarse_beta

    focal_fine_richness_predicted = numpy.asarray([sum(1-((1+focal_fine_theta*totreads_i)**(-1*focal_fine_beta))) for totreads_i in total_reads_fine])
    nonfocal_coarse_richness_predicted = numpy.asarray([sum(1-((1+nonfocal_coarse_theta*totreads_i)**(-1*nonfocal_coarse_beta))) for totreads_i in total_reads_fine])

    return focal_fine_richness_predicted, nonfocal_coarse_richness_predicted





def make_dbd_slm_dict():

    dbd_dict = {}

    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        sys.stderr.write("Loading taxon dict for %s...\n" % environment)
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied=rarefied)
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
            coarse_s_by_s = numpy.add.reduceat(s_by_s_fine, coarse_idx, axis=0)       

            total_reads_coarse = numpy.sum(coarse_s_by_s, axis=0)
            rel_coarse_s_by_s = coarse_s_by_s/numpy.sum(coarse_s_by_s, axis=0)
            mean_rel_s_by_s_coarse = numpy.mean(rel_coarse_s_by_s, axis=1)
            var_rel_s_by_s_coarse = numpy.var(rel_coarse_s_by_s, axis=1)

            # get index for fine for null
            #unique_fine, counts_fine = numpy.unique(all_fine, return_counts=True)
            #fine_idx = numpy.append([0], numpy.cumsum(counts_fine))[:-1]

            slope_all = []
            slope_slm_all = []

            for focal_coarse_idx, focal_coarse in enumerate(set_coarse):

                # ignore coarse-grained taxa with less than five fine-grained taxa
                if counts_coarse[focal_coarse_idx] < 5:
                    continue

                # all the fine-scale indices for the focal
                focal_coarse_s_by_s_idx = numpy.asarray(idx_to_sort[focal_coarse_idx])
                focal_fine_s_by_s = s_by_s_fine[focal_coarse_s_by_s_idx,:]

                focal_mean_fine_rel_s_by_s = mean_rel_s_by_s_fine[focal_coarse_s_by_s_idx]
                focal_var_fine_rel_s_by_s = var_rel_s_by_s_fine[focal_coarse_s_by_s_idx]

                nonfocal_coarse_s_by_s = numpy.delete(coarse_s_by_s, focal_coarse_idx, axis=0)
                
                nonfocal_mean_coarse_rel_s_by_s = numpy.delete(mean_rel_s_by_s_coarse, focal_coarse_idx, axis=0)          
                nonfocal_var_coarse_rel_s_by_s = numpy.delete(var_rel_s_by_s_coarse, focal_coarse_idx, axis=0)      

                focal_fine_richness_predicted, nonfocal_coarse_richness_predicted = predict_richness_dbd(focal_mean_fine_rel_s_by_s, focal_var_fine_rel_s_by_s, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, nonfocal_var_coarse_rel_s_by_s, total_reads_coarse)    

                #n_fine_focal_coarse = len(focal_coarse_s_by_s_idx)
                measure_coarse = numpy.sum(nonfocal_coarse_s_by_s > 0, axis=0)
                measure_fine = numpy.sum(focal_fine_s_by_s > 0, axis=0)

                # remove sites where there are no observations in either focal or non-focal
                idx_to_remove = (measure_fine>0) | (measure_coarse>0)

                measure_coarse = measure_coarse[idx_to_remove]
                measure_fine = measure_fine[idx_to_remove]

                nonfocal_coarse_richness_predicted = nonfocal_coarse_richness_predicted[idx_to_remove]
                focal_fine_richness_predicted = focal_fine_richness_predicted[idx_to_remove]

                if len(measure_fine) < 5:
                    continue

                slope, intercept, r_value, p_value, std_err = stats.linregress(measure_coarse, measure_fine)
                slope_slm, intercept_slm, r_value_slm, p_value_slm, std_err_slm = stats.linregress(nonfocal_coarse_richness_predicted, focal_fine_richness_predicted)

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

            dbd_dict[environment]['taxon'][coarse_rank] = {}
            dbd_dict[environment]['taxon'][coarse_rank]['mean_slope'] = mean_slope_all
            dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm'] = mean_slope_slm_all
            
            dbd_dict[environment]['taxon'][coarse_rank]['slope_all'] = slope_all.tolist()
            dbd_dict[environment]['taxon'][coarse_rank]['slope_slm_all'] = slope_slm_all.tolist()


        # phylogenetic coarse-graining
        sys.stderr.write("Subsetting samples for %s...\n" % environment)
        coarse_grained_tree_dict_env = coarse_grained_tree_dict[environment]
        sad_annotated_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
        samples = sad_annotated_dict['samples']

        taxa = numpy.asarray(list(sad_annotated_dict['taxa'].keys()))
        s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])

        distances = list(coarse_grained_tree_dict_env.keys())
        distances.sort()

        distances = numpy.asarray(distances)

        for distance_idx in range(len(distances) - 1):
            
            fine_distance = distances[distance_idx]
            coarse_distance = distances[distance_idx + 1]

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

                richness_fine_focal = numpy.sum(s_by_s_fine_focal>0, axis=0)
                richness_coarse_focal = numpy.sum(s_by_s_coarse_nonfocal>0, axis=0)

                # predict using SLM
                focal_mean_fine_rel_s_by_s = mean_rel_s_by_s_fine[f_idx_all]
                focal_var_fine_rel_s_by_s = var_rel_s_by_s_fine[f_idx_all]

                nonfocal_mean_coarse_rel_s_by_s = numpy.delete(mean_rel_s_by_s_coarse, focal_coarse_idx, axis=0)          
                nonfocal_var_coarse_rel_s_by_s = numpy.delete(var_rel_s_by_s_coarse, focal_coarse_idx, axis=0)      

                focal_fine_richness_predicted, nonfocal_coarse_richness_predicted = predict_richness_dbd(focal_mean_fine_rel_s_by_s, focal_var_fine_rel_s_by_s, total_reads_fine, nonfocal_mean_coarse_rel_s_by_s, nonfocal_var_coarse_rel_s_by_s, total_reads_coarse)    

                slope, intercept, r_value, p_value, std_err = stats.linregress(richness_coarse_focal, richness_fine_focal)
                slope_slm, intercept_slm, r_value_slm, p_value_slm, std_err_slm = stats.linregress(nonfocal_coarse_richness_predicted, focal_fine_richness_predicted)

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

            dbd_dict[environment]['phylo'][coarse_distance] = {}
            dbd_dict[environment]['phylo'][coarse_distance]['mean_slope'] = mean_slope_all
            dbd_dict[environment]['phylo'][coarse_distance]['mean_slope_slm'] = mean_slope_slm_all

            dbd_dict[environment]['phylo'][coarse_distance]['slope_all'] = slope_all.tolist()
            dbd_dict[environment]['phylo'][coarse_distance]['slope_slm_all'] = slope_slm_all.tolist()

            #print(len(fine_grained_list), len(coarse_grained_list), s_by_s_coarse.shape, s_by_s_fine.shape)

    sys.stderr.write("Saving diversity dictionary...\n")
    with open(richness_dbd_dict_path, 'wb') as handle:
        pickle.dump(dbd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    







def plot_dbd_slm_dict():

    with open(richness_dbd_dict_path, 'rb') as handle:
        dbd_dict = pickle.load(handle)

    fig = plt.figure(figsize = (8, 4)) #
    fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

    ax_taxon = plt.subplot2grid((1, 2), (0,0))
    ax_phylo = plt.subplot2grid((1, 2), (0,1))

    ax_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_taxon.transAxes)
    ax_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_phylo.transAxes)
    
    for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

        distances = list(dbd_dict[environment]['phylo'].keys())
        distances.sort()
        taxa_ranks = list(dbd_dict[environment]['taxon'].keys())

        color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)
        color_environment_phylo = plot_utils.get_custom_cmap_phylo(environment, len(distances))

        for coarse_rank_idx, coarse_rank in enumerate(taxa_ranks):
            
            observed = dbd_dict[environment]['taxon'][coarse_rank]['mean_slope']
            predicted = dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm']

            color_taxon = color_environment_taxon[coarse_rank_idx+1]

            ax_taxon.scatter(observed, predicted, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)

        for d_idx, d in enumerate(distances):

            observed = dbd_dict[environment]['phylo'][d]['mean_slope']
            predicted = dbd_dict[environment]['phylo'][d]['mean_slope_slm']

            color_phylo = color_environment_phylo[d_idx]
            ax_phylo.scatter(observed, predicted, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)


    ax_taxon.set_xscale('log', base=10)
    ax_taxon.set_yscale('log', base=10)
    #ax.plot([min_,max_], [min_,max_], lw=2,ls='--',c='k',zorder=1, label='1:1')
    #ax.set_xlim([min_,max_])
    #ax.set_ylim([min_,max_])
    ax_taxon.plot([0.004,0.7], [0.004,0.7], lw=2,ls='--',c='k',zorder=1, label='1:1')
    ax_taxon.set_xlabel("Observed mean DBD slope", fontsize = 12)
    ax_taxon.set_ylabel("Predicted mean DBD slope", fontsize = 12)
    ax_taxon.legend(loc="upper left", fontsize=7)
    ax_taxon.set_title("Taxonomic coarse-graining", fontsize=12)



    ax_phylo.set_xscale('log', base=10)
    ax_phylo.set_yscale('log', base=10)
    ax_phylo.plot([0.004,0.7], [0.004,0.7], lw=2,ls='--',c='k',zorder=1, label='1:1')
    ax_phylo.set_xlabel("Observed mean DBD slope", fontsize = 12)
    ax_phylo.set_ylabel("Predicted mean DBD slope", fontsize = 12)
    ax_phylo.set_title("Phylogenetic coarse-graining", fontsize=12)
    #ax_phylo.legend(loc="lower left", fontsize=7)



 
    fig.subplots_adjust(hspace=0.3,wspace=0.35)
    fig_name = "%sdbd_slope_slm.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



#make_dbd_slm_dict()


plot_dbd_slm_dict()

