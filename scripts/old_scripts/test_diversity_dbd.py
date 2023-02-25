import dbd_utils
import config
import pickle
import numpy
import matplotlib.pyplot as plt
import diversity_utils
import scipy.stats as stats

diversity_dbd_dict_path = config.data_directory + "diversity_dbd_dict.pickle"




dbd_utils.make_richness_diversity_prediction_phylo_dict('soil metagenome')

# running
#dbd_utils.make_diversity_slm_dbd_dict()
#dbd_utils.make_diversity_slm_dbd_simulation_taxon_dict()

# to do
#dbd_utils.make_diversity_slm_dbd_simulation_dict()

#dbd_utils.make_diversity_slm_dbd_simulation_phylo_dict()



with open(diversity_dbd_dict_path, 'rb') as handle:
    dbd_slope_dict = pickle.load(handle)


fig, ax = plt.subplots(figsize=(4,4))


environment_all = list(dbd_slope_dict.keys())

slope_all = []
slope_slm_all = []

for e in environment_all:

    for key, value in dbd_slope_dict[e].items():

        for rank_idx, rank in enumerate(diversity_utils.taxa_ranks):

            if rank in value:

                if 'slope_slm_all' in value[rank]:

                    slope_all.extend(value[rank]['slope_all'])
                    slope_slm_all.extend(value[rank]['slope_slm_all'])

                    mean_slm_slople_all = value[rank]['slope_slm_all']

                    mean_mean_slm_slople_all = numpy.mean(mean_slm_slople_all)
                    mean_slope_all = numpy.mean( value[rank]['slope_all'])

                    rho_rank = numpy.corrcoef(value[rank]['slope_all'], mean_slm_slople_all)[0,1]

                    #error_rank = numpy.mean(numpy.absolute( numpy.asarray(value[rank]['slope_all']) - numpy.asarray(mean_slm_slople_all) )/ numpy.absolute(numpy.asarray(value[rank]['slope_all']) ))
                    #print(error_rank)

                    #print(rank, len(mean_slm_slople_all), rho_rank**2)

                    ax.scatter(mean_slope_all, mean_mean_slm_slople_all, c='k', alpha=0.5)



#slope, intercept, r_valuer_value, p_value, std_err = stats.linregress(slope_all, slope_slm_all)
rho = numpy.corrcoef(slope_all, slope_slm_all)[0,1]

print('r^2 = ', rho**2)



#ax.set_ylabel("Mean relative errror", fontsize = 10)

ax.plot([-0.2,0.7], [-0.2,0.7], lw=2,ls='--',c='k',zorder=1, label='1:1')
ax.set_xlim([-0.2,0.7])
ax.set_ylim([-0.2,0.7])
ax.set_xlabel("Observed diversity DBD slope", fontsize = 12)
ax.set_ylabel("Predicted diversity DBD slope", fontsize = 12)





fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%sdiversity_dbd_slope_slm_test.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()


