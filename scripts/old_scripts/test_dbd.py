import os
#from Bio import Phylo
import random
import copy
import sys
import numpy
import random
import pickle
import scipy.stats as stats
import matplotlib.pyplot as plt
from collections import Counter

import functools
import operator

#import dbd_utils

#import diversity_utils
#import config


#import treeswift


random.seed(123456789)
iter = 1000






def plot_slope():


    dbd_slope_dict = dbd_utils.load_slope_dicersity_dict()

    genera_all = list(dbd_slope_dict.keys())

    slopes_all = [dbd_slope_dict[g]['slope'] for g in genera_all]
    standard_score_all = [dbd_slope_dict[g]['standard_score'] for g in genera_all]
    standard_score_all = numpy.asarray(standard_score_all)

    #print(sum(standard_score_all>0) / len(standard_score_all))

    # 0.44 genera have slopes greater than the bulk of the null distribution


    fig = plt.figure(figsize = (5, 8)) #
    fig.subplots_adjust(bottom= 0.15,  wspace=0.25)

    ax_hist = plt.subplot2grid((2, 1), (0, 0), colspan=1)
    ax_bar = plt.subplot2grid((2, 1), (1, 0), colspan=1)



    null_dist_all = [dbd_slope_dict[g]['null_dist'] for g in genera_all]
    null_dist_all = numpy.asarray(functools.reduce(operator.iconcat, null_dist_all, []))

    hist_, bin_edges_ = numpy.histogram(null_dist_all, density=True, bins=15)
    bins_mean_ = [0.5 * (bin_edges_[i] + bin_edges_[i+1]) for i in range(0, len(bin_edges_)-1 )]
    bins_mean_ = numpy.asarray(bins_mean_)
    #hist_to_plot = hist_[hist_>0]
    #bins_mean_to_plot = bins_mean_[hist_>0]
    hist_to_plot = hist_
    hist_to_plot = hist_to_plot/sum(hist_to_plot)
    ax_hist.plot(bins_mean_, hist_to_plot, alpha=1, c='k', zorder=2, label='Null')


    #for g in genera_all:

    #    null_dist = dbd_slope_dict[g]['null_dist']

    #    hist_, bin_edges_ = numpy.histogram(null_dist, density=True, bins=20)
    #    bins_mean_ = [0.5 * (bin_edges_[i] + bin_edges_[i+1]) for i in range(0, len(bin_edges_)-1) ]
    #    bins_mean_ = numpy.asarray(bins_mean_)

    #    #hist_to_plot = hist_[hist_>0]
    #    hist_to_plot = hist_
    #    #ins_mean_to_plot = bins_mean_[hist_>0]

    #    hist_to_plot = hist_to_plot/sum(hist_to_plot)

    #    ax_hist.plot(bins_mean_, hist_to_plot, alpha=0.2, c='k', zorder=1)
    #    #ax_hist.hist(dbd_slope_dict[g]['null_dist'], histtype='step', alpha=0.3, bins= 20, color='k', density=True, zorder=1)

    #hist_, bin_edges_ = numpy.histogram(slopes_all, density=True, bins=30)
    hist_, bin_edges_ = numpy.histogram(slopes_all, density=True, bins=10)
    bins_mean_ = [0.5 * (bin_edges_[i] + bin_edges_[i+1]) for i in range(0, len(bin_edges_)-1 )]
    bins_mean_ = numpy.asarray(bins_mean_)
    #hist_to_plot = hist_[hist_>0]
    #bins_mean_to_plot = bins_mean_[hist_>0]
    hist_to_plot = hist_
    hist_to_plot = hist_to_plot/sum(hist_to_plot)
    ax_hist.plot(bins_mean_, hist_to_plot, alpha=1, c='b', zorder=2, label='Observed')

    ax_hist.legend(loc="upper left", fontsize=8)
    #ax_hist.axvline(x=dbd_slope_dict['Lactobacillus']['slope'], color='k', linestyle=':', lw = 3, zorder=2)



    #ax_hist.hist(slopes_all, alpha=0.8, bins= 20, density=True, zorder=2)
    #fig, ax = plt.subplots(figsize=(4,4))
    #ax.hist(slope_null_all, alpha=0.8, bins=20, density=True)
    #ax_hist.axvline(x=numpy.median(slopes_all), color='dodgerblue', linestyle='--', lw = 3, zorder=2, label='Observed')
    #ax_hist.axvline(x=0, color='k', linestyle=':', lw = 3, zorder=2)
    #ax.legend(loc="upper left", fontsize=8)
    ax_hist.set_xlabel('DBD slope', fontsize=12)
    ax_hist.set_ylabel('Probability', fontsize=12)

    # plot bar plot
    zipped = list(zip(genera_all, slopes_all))
    zipped_sort = sorted(zipped, key = lambda x: x[1])



    for i in range(len(zipped_sort)):
        g_i = zipped_sort[i][0]
        #ax.plot(zipped_sort[1], i, c='k', ls=utils.line_dict[phage_treatment_type], label=label, lw=2, alpha=1)
        p_value_i = dbd_slope_dict[g_i]['p_value']
        if p_value_i < 0.05:
            c = 'r'
        else:
            c='k'

        ax_bar.scatter(i, zipped_sort[i][1], c=c, s=5, zorder=2)
        delta = 0.1
        fill_between_x = numpy.asarray([i-delta, i+delta])

        ax_bar.fill_between(fill_between_x, y1=dbd_slope_dict[g_i]['lower_ci'] , y2=dbd_slope_dict[g_i]['upper_ci'] , color='grey', zorder=1)


    ax_bar.set_ylabel('DBD slope', fontsize=12)
    ax_bar.set_xlabel('Genera', fontsize=12)



    fig.subplots_adjust(wspace=0.3, hspace=0.3)
    fig.savefig("%sdbd_slope.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.5, dpi = 600)
    plt.close()


plot_slope()
