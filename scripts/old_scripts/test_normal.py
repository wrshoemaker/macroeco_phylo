import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import config
import diversity_utils

# Generate some data for this demonstration.
data = norm.rvs(10.0, 2.5, size=10000)

# Fit a normal distribution to the data:
mu, std = norm.fit(data)

data = data[data>0]
log_data = np.log(data)
rescaled_log_data = (log_data - np.mean(log_data))/np.std(log_data)
# Plot the histogram.

fig, ax = plt.subplots(figsize=(4,4))

hist_to_plot, bins_mean_to_plot = diversity_utils.get_hist_and_bins(rescaled_log_data)

ax.plot(bins_mean_to_plot, hist_to_plot, color='k', alpha=1, lw=1)

# Plot the PDF.
#xmin, xmax = plt.xlim()
#x = np.linspace(xmin, xmax, 100)
#p = norm.pdf(x, mu, std)
#plt.plot(x, p, 'k', linewidth=2)
#title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
#plt.title(title)

ax.set_yscale('log', base=10)

fig.subplots_adjust(hspace=0.3,wspace=0.35)
fig_name = "%stest_normal.png" % config.analysis_directory
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
