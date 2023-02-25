import numpy
import config
import diversity_utils
import dbd_utils
import matplotlib.pyplot as plt
import plot_utils
import scipy.stats as stats



dbd_taxon_dict = dbd_utils.load_diversity_slm_dbd_simulation_taxon_dict()
dbd_phylo_dict = dbd_utils.load_svd_and_make_diversity_slm_dbd_simulation_phylo_dict()



fig = plt.figure(figsize = (8.5, 4))
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

# (#rows, # columns)
ax_taxon = plt.subplot2grid((1, 2), (0, 0))
ax_phylo = plt.subplot2grid((1, 2), (0, 1))

ax_taxon.text(-0.1, 1.04, plot_utils.sub_plot_labels[0], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_taxon.transAxes)
ax_phylo.text(-0.1, 1.04, plot_utils.sub_plot_labels[1], fontsize=10, fontweight='bold', ha='center', va='center', transform=ax_phylo.transAxes)




taxon_slope_all = []
taxon_slope_slm_all = []
taxon_mean_slope_all = []
taxon_mean_slope_slm_all = []
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    if environment not in dbd_taxon_dict:
        continue

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

    for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

        if rank not in dbd_taxon_dict[environment]['taxon']:
            continue


        mean_slope = dbd_taxon_dict[environment]['taxon'][rank]['mean_slope']
        mean_slope_slm = dbd_taxon_dict[environment]['taxon'][rank]['mean_slope_slm']

        color_taxon = color_environment_taxon[rank_idx+1]

        ax_taxon.scatter(mean_slope, mean_slope_slm, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)

        taxon_slope_all.extend(dbd_taxon_dict[environment]['taxon'][rank]['slope_all'])
        taxon_slope_slm_all.extend(dbd_taxon_dict[environment]['taxon'][rank]['slope_slm_all'])
        taxon_mean_slope_all.append(mean_slope)
        taxon_mean_slope_slm_all.append(mean_slope_slm)






# phylogenetic coarse-graining
phylo_slope_all = []
phylo_slope_slm_all = []
phylo_mean_slope_all = []
phylo_mean_slope_slm_all = []
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    if environment not in dbd_phylo_dict:
        continue

    distances = list(dbd_phylo_dict[environment]['phylo'].keys())
    distances.sort()

    color_environment_phylo = plot_utils.get_custom_cmap_phylo(environment, n=len(distances))

    for distance_idx, distance in enumerate(distances):


        mean_slope = dbd_phylo_dict[environment]['phylo'][distance]['mean_slope']
        mean_slope_slm = dbd_phylo_dict[environment]['phylo'][distance]['mean_slope_slm']

        color_phylo = color_environment_phylo[distance_idx+1]

        ax_phylo.scatter(mean_slope, mean_slope_slm, color=color_phylo, s=40, linewidth=0.8, edgecolors='k', zorder=2)

        phylo_slope_all.extend(dbd_phylo_dict[environment]['phylo'][distance]['slope_all'])
        phylo_slope_slm_all.extend(dbd_phylo_dict[environment]['phylo'][distance]['slope_slm_all'])
        phylo_mean_slope_all.append(mean_slope)
        phylo_mean_slope_slm_all.append(mean_slope_slm)




taxon_rho = numpy.corrcoef(taxon_slope_all, taxon_slope_slm_all)[0,1]

taxon_min = min(taxon_mean_slope_all+taxon_mean_slope_slm_all)*0.8
taxon_max = max(taxon_mean_slope_all+taxon_mean_slope_slm_all)*1.2

ax_taxon.plot([taxon_min,taxon_max], [taxon_min,taxon_max], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax_taxon.set_xlim([taxon_min,taxon_max])
ax_taxon.set_ylim([taxon_min,taxon_max])
ax_taxon.set_xlabel("Observed mean DBD diversity slope", fontsize = 12)
ax_taxon.set_ylabel("Predicted mean DBD diversity slope, simulation", fontsize = 11)
ax_taxon.legend(loc="upper left", fontsize=7)
ax_taxon.set_title("Taxonomic coarse-graining", fontsize=12)
ax_taxon.text(0.2, 0.85, r'$\rho^{2} =$' + str(round(taxon_rho, 3)), fontsize=10, ha='center', va='center', transform=ax_taxon.transAxes)



phylo_rho = numpy.corrcoef(phylo_slope_all, phylo_slope_slm_all)[0,1]

phylo_min = min(phylo_mean_slope_all+phylo_mean_slope_slm_all)*0.8
phylo_max = max(phylo_mean_slope_all+phylo_mean_slope_slm_all)*1.2

ax_phylo.plot([phylo_min,phylo_max], [phylo_min,phylo_max], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax_phylo.set_xlim([phylo_min,phylo_max])
ax_phylo.set_ylim([phylo_min,phylo_max])
ax_phylo.set_xlabel("Observed mean DBD diversity slope", fontsize = 12)
ax_phylo.set_ylabel("Predicted mean DBD diversity slope, simulation", fontsize = 11)
ax_phylo.set_title("Phylogenetic coarse-graining", fontsize=12)
ax_phylo.text(0.2, 0.85, r'$\rho^{2} =$' + str(round(phylo_rho, 3)), fontsize=10, ha='center', va='center', transform=ax_phylo.transAxes)




fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sdbd_diversity_slm_simulation.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
