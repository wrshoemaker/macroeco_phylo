import numpy
import diversity_utils
from matplotlib import cm
import matplotlib as mpl
from matplotlib import colors
from matplotlib.lines import Line2D



sub_plot_labels = ['a','b','c', 'd','e','f', 'g','h','i', 'j','k','l', 'm','n','o']

# taxonomic hierarchy colors
cmap_offset = int(0.2*16)
# +cmap_offset
rgb_blue_taxon = cm.Blues(numpy.linspace(0,1,len(diversity_utils.taxa_ranks)+3))
rgb_blue_taxon = mpl.colors.ListedColormap(rgb_blue_taxon[cmap_offset:,:-1])


rgb_blue_phylo = cm.Blues(numpy.linspace(0,1,100+10))
rgb_blue_phylo = mpl.colors.ListedColormap(rgb_blue_phylo[cmap_offset:,:-1])








environment_color_map = {'marine metagenome': colors.hex2color(colors.cnames['seagreen']),
                    'human gut metagenome': colors.hex2color(colors.cnames['aquamarine']),
                    'human oral metagenome': colors.hex2color(colors.cnames['dodgerblue']),
                    'freshwater sediment metagenome': colors.hex2color(colors.cnames['cadetblue']),
                    'microbial mat metagenome': colors.hex2color(colors.cnames['darkorange']),
                    'human skin metagenome': colors.hex2color(colors.cnames['sandybrown']),
                    'freshwater metagenome': colors.hex2color(colors.cnames['orchid']),
                    'soil metagenome': colors.hex2color(colors.cnames['darkturquoise']),
                    'marine sediment metagenome': colors.hex2color(colors.cnames['darkred'])}









#label='abc', color='red', linewidth=1.5, marker='o',
#              markerfacecolor='yellow', markeredgewidth=1.5, markersize=16




def make_blue_cmap(n):

    rgb_blue_phylo = cm.Blues(numpy.linspace(0,1,n+3))
    rgb_blue_phylo = mpl.colors.ListedColormap(rgb_blue_phylo[cmap_offset:,:-1])

    return rgb_blue_phylo



def CustomCmap(from_rgb,to_rgb):

    # from color r,g,b
    r1,g1,b1 = from_rgb

    # to color r,g,b
    r2,g2,b2 = to_rgb

    cdict = {'red': ((0, r1, r1),
                   (1, r2, r2)),
           'green': ((0, g1, g1),
                    (1, g2, g2)),
           'blue': ((0, b1, b1),
                   (1, b2, b2))}

    cmap = colors.LinearSegmentedColormap('custom_cmap', cdict)
    return cmap



def get_custom_cmap_taxon(environment):

    color = environment_color_map[environment]

    custom_cmap = CustomCmap([1.00, 1.00, 1.00], color)

    rgb_ = custom_cmap(numpy.linspace(0, 1, len(diversity_utils.taxa_ranks)+2))

    return rgb_



def get_custom_cmap_phylo(environment, n):

    color = environment_color_map[environment]

    custom_cmap = CustomCmap([1.00, 1.00, 1.00], color)

    rgb_ = custom_cmap(numpy.linspace(0, 1, n+1))

    return rgb_






# https://github.com/weecology/macroecotools/blob/master/macroecotools/macroecotools.py
# code to cluster points
def count_pts_within_radius(x, y, radius, logscale=0):
    """Count the number of points within a fixed radius in 2D space"""
    #TODO: see if we can improve performance using KDTree.query_ball_point
    #http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query_ball_point.html
    #instead of doing the subset based on the circle
    unique_points = set([(x[i], y[i]) for i in range(len(x))])
    count_data = []
    logx, logy, logr = numpy.log10(x), numpy.log10(y), numpy.log10(radius)
    for a, b in unique_points:
        if logscale == 1:
            loga, logb = numpy.log10(a), numpy.log10(b)
            num_neighbors = len(x[((logx - loga) ** 2 +
                                   (logy - logb) ** 2) <= logr ** 2])
        else:
            num_neighbors = len(x[((x - a) ** 2 + (y - b) ** 2) <= radius ** 2])
        count_data.append((a, b, num_neighbors))
    return count_data



def plot_color_by_pt_dens(x, y, radius, loglog=0):
    """Plot bivariate relationships with large n using color for point density

    Inputs:
    x & y -- variables to be plotted
    radius -- the linear distance within which to count points as neighbors
    loglog -- a flag to indicate the use of a loglog plot (loglog = 1)

    The color of each point in the plot is determined by the logarithm (base 10)
    of the number of points that occur with a given radius of the focal point,
    with hotter colors indicating more points. The number of neighboring points
    is determined in linear space regardless of whether a loglog plot is
    presented.
    """
    plot_data = count_pts_within_radius(x, y, radius, loglog)
    sorted_plot_data = numpy.array(sorted(plot_data, key=lambda point: point[2]))

    return sorted_plot_data




def get_scatter_density_arrays_for_linearlog(x, y, color_radius=2):


    #idx_to_keep = (x>0) & (y > 0)
    #x = x[idx_to_keep]
    #y = y[idx_to_keep]

    x = numpy.log10(x)

    x_and_y = numpy.concatenate((x,y),axis=0)
    #min_ = min(x_and_y)
    #max_ = max(x_and_y)

    sorted_plot_data = plot_color_by_pt_dens(x, y, radius=color_radius, loglog=0)
    x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]

    return x, y, z



def get_scatter_density_arrays_for_loglog(x, y, color_radius=2):


    idx_to_keep = (x>0) & (y > 0)
    x = x[idx_to_keep]
    y = y[idx_to_keep]

    x_and_y = numpy.concatenate((x,y),axis=0)
    min_ = min(x_and_y)
    max_ = max(x_and_y)

    sorted_plot_data = plot_color_by_pt_dens(x, y, radius=color_radius, loglog=1)
    x,y,z = sorted_plot_data[:, 0], sorted_plot_data[:, 1], sorted_plot_data[:, 2]

    return x, y, z



def get_bin_mean_x_y(x, y, bins=20, min_n_bin=5):

    x_log10 = numpy.log10(x)
    y_log10 = numpy.log10(y)

    hist_all, bin_edges_all = numpy.histogram(x_log10, density=True, bins=bins)
    #bins_x = [0.5 * (bin_edges_all[i] + bin_edges_all[i+1]) for i in range(0, len(bin_edges_all)-1 )]
    bins_x_to_keep = []
    bins_y = []
    for i in range(0, len(bin_edges_all)-1 ):
        y_log10_i = y_log10[(x_log10>=bin_edges_all[i]) & (x_log10<bin_edges_all[i+1])]

        if len(y_log10_i) >= min_n_bin:
            bins_x_to_keep.append(bin_edges_all[i])
            bins_y.append(numpy.median(y_log10_i))


    bins_x_to_keep = numpy.asarray(bins_x_to_keep)
    bins_y = numpy.asarray(bins_y)

    bins_x_to_keep_no_nan = bins_x_to_keep[(~numpy.isnan(bins_x_to_keep)) & (~numpy.isnan(bins_y))]
    bins_y_no_nan = bins_y[(~numpy.isnan(bins_x_to_keep)) & (~numpy.isnan(bins_y))]

    bins_x_to_keep_no_nan = 10**bins_x_to_keep_no_nan
    bins_y_no_nan = 10**bins_y_no_nan

    return bins_x_to_keep_no_nan, bins_y_no_nan




legend_elements = [Line2D([0], [0], marker='o', color='w', label='Marine', markerfacecolor=get_custom_cmap_phylo('marine metagenome', 20)[15], markersize=7, markeredgewidth=1.4, markeredgecolor='k'),
                   Line2D([0], [0], marker='o', color='w', label='Marine sediment', markerfacecolor=get_custom_cmap_phylo('marine sediment metagenome', 20)[15], markersize=7, markeredgewidth=1.4, markeredgecolor='k'),
                   Line2D([0], [0], marker='o', color='w', label='Human gut', markerfacecolor=get_custom_cmap_phylo('human gut metagenome', 20)[15], markersize=7, markeredgewidth=1.4, markeredgecolor='k'),
                   Line2D([0], [0], marker='o', color='w', label='Human oral', markerfacecolor=get_custom_cmap_phylo('human oral metagenome', 20)[15], markersize=7, markeredgewidth=1.4, markeredgecolor='k'),
                   Line2D([0], [0], marker='o', color='w', label='Human skin', markerfacecolor=get_custom_cmap_phylo('human skin metagenome', 20)[15], markersize=7, markeredgewidth=1.4, markeredgecolor='k'),
                   Line2D([0], [0], marker='o', color='w', label='Freshwater sediment', markerfacecolor=get_custom_cmap_phylo('freshwater sediment metagenome', 20)[15], markersize=7, markeredgewidth=1.4, markeredgecolor='k'),
                   Line2D([0], [0], marker='o', color='w', label='Microbial mat', markerfacecolor=get_custom_cmap_phylo('microbial mat metagenome', 20)[15], markersize=7, markeredgewidth=1.4, markeredgecolor='k'),
                   Line2D([0], [0], marker='o', color='w', label='Freshwater', markerfacecolor=get_custom_cmap_phylo('freshwater metagenome', 20)[15], markersize=7, markeredgewidth=1.4, markeredgecolor='k'),
                   Line2D([0], [0], marker='o', color='w', label='Soil', markerfacecolor=get_custom_cmap_phylo('soil metagenome', 20)[15], markersize=7, markeredgewidth=1.4, markeredgecolor='k')]




