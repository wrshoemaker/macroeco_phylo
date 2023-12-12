import config
import numpy
import sys
import dbd_utils
import dbd_richness_neutral
import matplotlib.pyplot as plt
import plot_utils
import diversity_utils
from matplotlib.lines import Line2D


dbd_slm_dict = dbd_utils.load_richness_slm_dbd_dict()
#dbd_phylo_dict = dbd_utils.load_svd_and_make_diversity_slm_dbd_simulation_phylo_dict()

dbd_dict = dbd_richness_neutral.load_richness_neutral_dbd_dict()


rarefied = False

fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


coarse_grained_tree_dict_all = {}
distances_all = []
for environment in diversity_utils.environments_to_keep:
    sys.stderr.write("Loading tree dict for %s...\n" % environment)
    coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)
    coarse_grained_tree_dict_all[environment] = coarse_grained_tree_dict





environments_to_keep = diversity_utils.environments_to_keep

environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        distances = list(coarse_grained_tree_dict.keys())
        distances.sort()

        color_environment_phylo = plot_utils.get_custom_cmap_phylo(environment, len(distances))

        distance_all = []
        mean_error_slope_all = []
        for distance_idx, distance in enumerate(distances):

            if distance not in dbd_dict[environment]['phylo']:
                continue

            mean_error_slope = dbd_dict[environment]['phylo'][distance]['mean_error_slope_neutral']

            color = color_environment_phylo[distance_idx+1]
            ax.scatter(distance, mean_error_slope, color=color, s=40, linewidth=0.8, edgecolors='k', marker='s', zorder=2)


            distance_all.append(distance)
            mean_error_slope_all.append(mean_error_slope)
        

        ax.plot(distance_all, mean_error_slope_all, lw=1, ls='-', c='k', zorder=1)


        distance_all = []
        mean_error_slope_all = []
        for distance_idx, distance in enumerate(distances):

            if distance not in dbd_slm_dict[environment]['phylo']:
                continue

            mean_error_slope = dbd_slm_dict[environment]['phylo'][distance]['mean_error_slope_slm']

            color = color_environment_phylo[distance_idx+1]
            ax.scatter(distance, mean_error_slope, color=color, s=40, linewidth=0.8, edgecolors='k', zorder=2)


            distance_all.append(distance)
            mean_error_slope_all.append(mean_error_slope)
        

        ax.plot(distance_all, mean_error_slope_all, lw=1, ls='-', c='k', zorder=1)
        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)

        ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)

        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.tick_params(axis='both', which='minor', labelsize=7)

        ax.set_xlabel("Phylogenetic distance", fontsize = 9)
        ax.set_ylabel("Mean relative error of richness slope", fontsize = 9)

        legend_elements = [Line2D([0], [0], marker='o', color='w', label='SLM',
                          markerfacecolor='none',markeredgecolor='k', markersize=10),
                            Line2D([0], [0], marker='s', color='w', label='UNTB',
                          markerfacecolor='none',markeredgecolor='k', markersize=10)]


        if (environment_chunk_idx == 0) and (environment_idx == 0):

            ax.legend(handles=legend_elements, loc='upper left')

            




fig.subplots_adjust(hspace=0.3,wspace=0.3)
fig_name = "%sfigure-5-figure-supplement-9-gamma_error_summary_neutral_phylo.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()

