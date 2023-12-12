import config
import numpy
import sys
import dbd_utils
import dbd_richness_neutral
import matplotlib.pyplot as plt
import plot_utils
import diversity_utils


dbd_taxon_dict = dbd_utils.load_diversity_slm_dbd_simulation_taxon_dict()
dbd_phylo_dict = dbd_utils.load_svd_and_make_diversity_slm_dbd_simulation_phylo_dict()

dbd_dict = dbd_richness_neutral.load_richness_neutral_dbd_dict()



fig = plt.figure(figsize = (8.5, 4))
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

# (#rows, # columns)
ax_taxon = plt.subplot2grid((1, 2), (0, 0))
ax_phylo = plt.subplot2grid((1, 2), (0, 1))

ax_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_taxon.transAxes)
ax_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_phylo.transAxes)




for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    if environment not in dbd_dict:
        continue

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

    rho_environment_taxon_all = []
    color_taxon_all = []
    coarse_rank_idx_to_keep = []
    rank_idx_all = []
    mean_error_slope_neutral_all = []
    for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

        if 'mean_error_slope_neutral' not in dbd_dict[environment]['taxon'][rank]:
            continue

        mean_error_slope_neutral = dbd_dict[environment]['taxon'][rank]['mean_error_slope_neutral']

        color_taxon = color_environment_taxon[rank_idx+1]

        ax_taxon.scatter(rank_idx, mean_error_slope_neutral, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)

        rank_idx_all.append(rank_idx)
        mean_error_slope_neutral_all.append(mean_error_slope_neutral)

    ax_taxon.plot(rank_idx_all, mean_error_slope_neutral_all, ls='-', lw=1, c='k', zorder=1)


    # plot for SLM




fig.subplots_adjust(hspace=0.2,wspace=0.25)
fig_name = "%sgamma_error_summary_neutral.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()

