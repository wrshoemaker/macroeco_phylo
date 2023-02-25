import config
import diversity_utils
import dbd_utils
import matplotlib.pyplot as plt
import plot_utils
import scipy.stats as stats
import numpy



#dbd_utils.make_richness_slAm_dbd_dict()
#dbd_utils.make_diversity_slm_dbd_dict()

#dbd_utils.make_svd_dict()


#dbd_utils.make_diversity_slm_dbd_simulation_dict()

#dbd_utils.make_diversity_slm_dbd_simulation_dict()
dbd_dict = dbd_utils.load_diversity_slm_dbd_simulation_dict()


fig, ax = plt.subplots()

slope_all = []
slope_slm_all = []

mean_slope_all = []
mean_slope_slm_all = []
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

    for coarse_rank_idx, coarse_rank in enumerate(diversity_utils.taxa_ranks):

        if coarse_rank in dbd_dict[environment]['taxon']:

            slope_all.extend(dbd_dict[environment]['taxon'][coarse_rank]['slope_all'])
            slope_slm_all.extend(dbd_dict[environment]['taxon'][coarse_rank]['slope_slm_all'])

            mean_slope = dbd_dict[environment]['taxon'][coarse_rank]['mean_slope']
            mean_slope_slm = dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm']

            color_taxon = color_environment_taxon[coarse_rank_idx+1]

            ax.scatter(mean_slope, mean_slope_slm, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)

            mean_slope_all.append(mean_slope)
            mean_slope_slm_all.append(mean_slope_slm)



#slope, intercept, r_valuer_value, p_value, std_err = stats.linregress(slope_all, slope_slm_all)
rho = numpy.corrcoef(slope_all, slope_slm_all)[0,1]

print('r^2 = ', rho**2)

min_ = min(mean_slope_all+mean_slope_slm_all)
max_ = max(mean_slope_all+mean_slope_slm_all)

ax.plot([min_,max_], [min_,max_], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax.set_xlim([min_,max_])
ax.set_ylim([min_,max_])
ax.set_xlabel("Observed mean DBD diversity slope", fontsize = 12)
ax.set_ylabel("Predicted mean DBD diversity slope, simulation", fontsize = 12)
ax.legend(loc="upper left", fontsize=7)
ax.set_title("Taxonomic coarse-graining", fontsize=12)


fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sdbd_diversity_slope_slm_simulation.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
