import config
import numpy
import sys
import dbd_utils
import dbd_richness_neutral
import matplotlib.pyplot as plt
import plot_utils
import diversity_utils
from matplotlib.lines import Line2D


dbd_taxon_dict = dbd_utils.load_richness_slm_dbd_dict()
#dbd_phylo_dict = dbd_utils.load_svd_and_make_diversity_slm_dbd_simulation_phylo_dict()

dbd_dict = dbd_richness_neutral.load_richness_neutral_dbd_dict()



fig = plt.figure(figsize = (10, 10)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)


taxa_ranks = diversity_utils.taxa_ranks_with_asv[1:]
idx_taxa_ranks = numpy.asarray(list(range(len(taxa_ranks))))
taxa_ranks_label = diversity_utils.taxa_ranks_label_with_asv[1:]



environments_to_keep = diversity_utils.environments_to_keep

environment_chunk_all = [environments_to_keep[x:x+3] for x in range(0, len(environments_to_keep), 3)]
for environment_chunk_idx, environment_chunk in enumerate(environment_chunk_all):

    for environment_idx, environment in enumerate(environment_chunk):

        ax = plt.subplot2grid((3, 3), (environment_chunk_idx, environment_idx))

        color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

        rank_idx_all = []
        mean_error_slope_neutral_all = []
        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            if 'mean_error_slope_neutral' not in dbd_dict[environment]['taxon'][rank]:
                continue

            mean_error_slope_neutral = dbd_dict[environment]['taxon'][rank]['mean_error_slope_neutral']

            color_taxon = color_environment_taxon[rank_idx+1]

            ax.scatter(rank_idx, mean_error_slope_neutral, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', marker='s', zorder=2)

            rank_idx_all.append(rank_idx)
            mean_error_slope_neutral_all.append(mean_error_slope_neutral)
        
        ax.plot(rank_idx_all, mean_error_slope_neutral_all, lw=1, ls='-', c='k', zorder=1)


        rank_slm_idx_all = []
        mean_error_slope_slm_all = []
        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            if 'mean_error_slope_slm' not in dbd_taxon_dict[environment]['taxon'][rank]:
                continue

            mean_error_slope_slm = dbd_taxon_dict[environment]['taxon'][rank]['mean_error_slope_slm']
            color_taxon = color_environment_taxon[rank_idx+1]
            ax.scatter(rank_idx, mean_error_slope_slm, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)

            rank_slm_idx_all.append(rank_idx)
            mean_error_slope_slm_all.append(mean_error_slope_slm)

        ax.plot(rank_slm_idx_all, mean_error_slope_slm_all, lw=1, ls='-', c='k', zorder=1)


        ax.set_yscale('log', base=10)

        ax.set_title(' '.join(environment.split(' ')[:-1]).capitalize(), fontsize=11)
        ax.set_xticks(idx_taxa_ranks)
        ax.set_xticklabels(taxa_ranks_label, fontsize=9)
        ax.tick_params(axis='y', which='major', labelsize=7)
        ax.tick_params(axis='y', which='minor', labelsize=7)

        ax.set_ylabel("Mean relative error of richness slope", fontsize = 9)

        legend_elements = [Line2D([0], [0], marker='o', color='w', label='SLM',
                          markerfacecolor='none',markeredgecolor='k', markersize=10),
                            Line2D([0], [0], marker='s', color='w', label='UNTB',
                          markerfacecolor='none',markeredgecolor='k', markersize=10)]


        if (environment_chunk_idx == 0) and (environment_idx == 0):

            ax.legend(handles=legend_elements, loc='upper right')

            




fig.subplots_adjust(hspace=0.2,wspace=0.3)
fig_name = "%sgamma_error_summary_neutral_taxon.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()

