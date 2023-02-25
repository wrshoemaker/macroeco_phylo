import config
import diversity_utils
import dbd_utils
import matplotlib.pyplot as plt
import plot_utils
import numpy
import scipy.stats as stats


environment = 'human gut metagenome'

color_environment_taxon = plot_utils.get_custom_cmap_taxon(environment)

dbd_richness_dict = dbd_utils.load_richness_slm_dbd_dict()
dbd_diversity_dict = dbd_utils.load_diversity_slm_dbd_dict()


#dbd_richness_dict = dbd_richness_dict['human_gut_n']

fig = plt.figure(figsize = (8, 4)) #
fig.subplots_adjust(bottom= 0.1,  wspace=0.15)

ax_observed = plt.subplot2grid((1, 2), (0,0))
ax_predicted = plt.subplot2grid((1, 2), (0,1))


focal_coarse_all = list(dbd_diversity_dict[environment]['taxon']['genus']['focal_coarse'].keys())

for f in focal_coarse_all:

    if (f in dbd_richness_dict[environment]['taxon']['genus']['focal_coarse']) and (f in dbd_diversity_dict[environment]['taxon']['genus']['focal_coarse']):

        
        richness_fine = dbd_richness_dict[environment]['taxon']['genus']['focal_coarse'][f]['measure_fine']
        richness_coarse = dbd_richness_dict[environment]['taxon']['genus']['focal_coarse'][f]['measure_coarse']
        richness_fine_predicted = dbd_richness_dict[environment]['taxon']['genus']['focal_coarse'][f]['measure_fine_predicted']
        richness_coarse_predicted = dbd_richness_dict[environment]['taxon']['genus']['focal_coarse'][f]['measure_coarse_predicted']

        diversity_fine = dbd_diversity_dict[environment]['taxon']['genus']['focal_coarse'][f]['measure_fine']
        diversity_coarse = dbd_diversity_dict[environment]['taxon']['genus']['focal_coarse'][f]['measure_coarse']
        diversity_fine_predicted = dbd_diversity_dict[environment]['taxon']['genus']['focal_coarse'][f]['measure_fine_predicted']
        diversity_coarse_predicted = dbd_diversity_dict[environment]['taxon']['genus']['focal_coarse'][f]['measure_coarse_predicted']

        color_taxon = color_environment_taxon[5]
        
        #richness_coarse_mean = numpy.mean(richness_coarse)
        #richness_coarse_predicted_mean = numpy.mean(richness_coarse_predicted)
        #diversity_coarse_mean = numpy.mean(diversity_coarse)
        #diversity_coarse_predicted_mean = numpy.mean(diversity_coarse_predicted)


        ax_observed.scatter(richness_fine, diversity_fine, alpha=0.9, s=8, color=color_taxon)
        ax_predicted.scatter(richness_fine_predicted, diversity_fine_predicted, alpha=0.5, s=9,  color=color_taxon)




ax_observed.set_xlabel("Observed richness, focal, genus", fontsize = 12)
ax_observed.set_ylabel("Observed diversity, focal, genus", fontsize = 12)

ax_predicted.set_xlabel("Predicted richness, focal, genus", fontsize = 12)
ax_predicted.set_ylabel("Predicted diversity, focal, genus", fontsize = 12)


fig.subplots_adjust(hspace=0.2,wspace=0.25)
fig_name = "%srichness_vs_diversity.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
