import  matplotlib.pyplot as plt
import config

evenness = []
symmetry = []


file = open('%stest_tree_sad_shape.csv' % (config.data_directory), 'r')
file.readline() # header
for line in file:
    line_split = line.strip().split(',')
    evenness.append(float(line_split[1]))
    symmetry.append(float(line_split[2]))

file.close()




fig = plt.figure(figsize = (4, 4)) #
fig.subplots_adjust(bottom= 0.15)

plt.scatter(evenness, symmetry)


plt.xlabel("Evenness", fontsize = 20)
plt.ylabel("Symmetry", fontsize = 20)

plt.xscale('log', basex=10)
plt.yscale('log', basey=10)


fig.subplots_adjust() #hspace=0.3, wspace=0.5
fig_name = '%stest_figure.png' % (config.analysis_directory)
fig.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
