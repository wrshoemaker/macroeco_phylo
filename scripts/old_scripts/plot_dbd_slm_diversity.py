import config
import diversity_utils
import dbd_utils
import matplotlib.pyplot as plt
import plot_utils
import scipy.stats as stats

#dbd_utils.make_diversity_slm_dbd_dict(p)


#dbd_utils.make_diversity_slm_dbd_integral_dict('human gut metagenome')



#dbd_dict = dbd_utils.load_diversity_slm_dbd_dict()


fig, ax = plt.subplots()


mean_slope_all = []
mean_slope_slm_all = []
for environment_idx, environment in enumerate(diversity_utils.environments_to_keep):

    color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

    for coarse_rank_idx, coarse_rank in enumerate(diversity_utils.taxa_ranks):

        #predic

        if coarse_rank in dbd_dict[environment]['taxon']:

            mean_slope = dbd_dict[environment]['taxon'][coarse_rank]['mean_slope']
            mean_slope_slm = dbd_dict[environment]['taxon'][coarse_rank]['mean_slope_slm']

            #print(mean_slope, mean_slope_slm)

            color_taxon = color_environment_taxon[coarse_rank_idx+1]

            ax.scatter(mean_slope, mean_slope_slm, color=color_taxon, s=40, linewidth=0.8, edgecolors='k', zorder=2)

            mean_slope_all.append(mean_slope)
            mean_slope_slm_all.append(mean_slope_slm)



slope, intercept, r_valuer_value, p_value, std_err = stats.linregress(mean_slope_all, mean_slope_slm_all)

print(r_valuer_value**2)

min_ = min(mean_slope_all+mean_slope_slm_all)
max_ = max(mean_slope_all+mean_slope_slm_all)

ax.plot([min_,max_], [min_,max_], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax.set_xlim([min_,max_])
ax.set_ylim([min_,max_])
ax.set_xlabel("Observed mean DBD diversity slope", fontsize = 12)
ax.set_ylabel("Predicted mean DBD diversity slope", fontsize = 12)
ax.legend(loc="upper left", fontsize=7)
ax.set_title("Taxonomic coarse-graining", fontsize=12)


fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sdbd_diversity_slope_slm.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()


