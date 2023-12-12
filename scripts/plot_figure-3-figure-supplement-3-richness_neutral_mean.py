import config
import numpy
import dbd_richness_neutral
import matplotlib.pyplot as plt
import plot_utils
import diversity_utils


dbd_dict = dbd_richness_neutral.load_richness_neutral_dbd_dict()


fig = plt.figure(figsize = (8.5, 4))
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

# (#rows, # columns)
ax_taxon = plt.subplot2grid((1, 2), (0, 0))
ax_phylo = plt.subplot2grid((1, 2), (0, 1))

ax_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_taxon.transAxes)
ax_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_phylo.transAxes)



taxon_mean_slope_all = []
taxon_mean_slope_slm_all = []
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    if environment not in dbd_dict:
        continue

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

    rho_environment_taxon_all = []
    color_taxon_all = []
    coarse_rank_idx_to_keep = []
    for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

        if rank not in dbd_dict[environment]['taxon']:
            continue

        if 'slope_all' not in dbd_dict[environment]['taxon'][rank]:
            continue

        
        mean_slope = dbd_dict[environment]['taxon'][rank]['observed_mean_richness']
        mean_slope_slm = dbd_dict[environment]['taxon'][rank]['predicted_mean_richness']

        color_taxon = color_environment_taxon[rank_idx+1]

        ax_taxon.scatter(mean_slope, mean_slope_slm, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)

        taxon_mean_slope_all.append(mean_slope)
        taxon_mean_slope_slm_all.append(mean_slope_slm)

        color_taxon_all.append(color_taxon)
        coarse_rank_idx_to_keep.append(rank_idx)


taxon_min = min(taxon_mean_slope_all+taxon_mean_slope_slm_all)*0.8
taxon_max = max(taxon_mean_slope_all+taxon_mean_slope_slm_all)*1.2

taxon_rho = numpy.corrcoef(numpy.log10(taxon_mean_slope_all), numpy.log10(taxon_mean_slope_slm_all))[0,1]



ax_taxon.plot([taxon_min,taxon_max], [taxon_min,taxon_max], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax_taxon.set_xlim([taxon_min,taxon_max])
ax_taxon.set_ylim([taxon_min,taxon_max])
ax_taxon.set_xlabel("Observed mean richness", fontsize = 12)
ax_taxon.set_ylabel("Predicted mean richness, UNTB", fontsize = 11)
ax_taxon.legend(loc="upper left", fontsize=7)
ax_taxon.set_title("Taxonomic coarse-graining", fontsize=12, fontweight='bold')
#ax_taxon.text(0.2, 0.85, r'$\rho^{2} =$' + str(round(taxon_rho, 3)), fontsize=10, ha='center', va='center', transform=ax_taxon.transAxes)

ax_taxon.set_xscale('log', base=10)
ax_taxon.set_yscale('log', base=10)


# ploy phylogenetic

# phylogenetic coarse-graining
phylo_mean_slope_all = []
phylo_mean_slope_slm_all = []
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    if environment not in dbd_dict:
        continue

    distances = list(dbd_dict[environment]['phylo'].keys())
    distances.sort()

    color_environment_phylo = plot_utils.get_custom_cmap_phylo(environment, n=len(distances))

    rho_environment_phylo_all = []
    color_diversity_all = []
    diversity_all = []
    for distance_idx, distance in enumerate(distances):

        mean_slope = dbd_dict[environment]['phylo'][distance]['observed_mean_richness']
        mean_slope_slm = dbd_dict[environment]['phylo'][distance]['predicted_mean_richness']

        if (mean_slope<= 0) or (mean_slope_slm <= 0):
            continue
    

        color_phylo = color_environment_phylo[distance_idx+1]

        ax_phylo.scatter(mean_slope, mean_slope_slm, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)

        phylo_mean_slope_all.append(mean_slope)
        phylo_mean_slope_slm_all.append(mean_slope_slm)

        color_diversity_all.append(color_phylo)
        diversity_all.append(distance)





phylo_rho = numpy.corrcoef(numpy.log10(phylo_mean_slope_all), numpy.log10(phylo_mean_slope_slm_all))[0,1]



phylo_min = min(phylo_mean_slope_all+phylo_mean_slope_slm_all)*0.8
phylo_max = max(phylo_mean_slope_all+phylo_mean_slope_slm_all)*1.2

ax_phylo.plot([phylo_min,phylo_max], [phylo_min,phylo_max], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax_phylo.set_xlim([phylo_min,phylo_max])
ax_phylo.set_ylim([phylo_min,phylo_max])
#ax_phylo.legend(loc="upper left", fontsize=9, frameon=False)

#ax_richness_taxon.legend(handles=plot_utils.legend_elements, loc='upper left', fontsize=7)

ax_phylo.set_xlabel("Observed mean richness", fontsize = 12)
ax_phylo.set_ylabel("Predicted mean richness, UNTB", fontsize = 11)
ax_phylo.set_title("Phylogenetic coarse-graining", fontsize=12, fontweight='bold')
#ax_phylo.text(0.17, 0.85, r'$\rho^{2} =$' + str(round(phylo_rho, 3)), fontsize=10, ha='center', va='center', transform=ax_phylo.transAxes)
ax_phylo.set_xscale('log', base=10)
ax_phylo.set_yscale('log', base=10)





fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sfigure-3-figure-supplement-3-richness_neutral_mean.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()