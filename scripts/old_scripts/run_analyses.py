import os
#from Bio import Phylo
import random
import copy
import sys
import numpy
import random

import tree_utils

#from ete3 import Tree

import ete3

import diversity_utils
import tree_utils
import config

#import treeswift

random.seed(123456789)

sys.stderr.write("Subsetting samples...\n")


samples = diversity_utils.subset_observations()

sys.stderr.write("Getting SADs...\n")

SADs, taxonomy_names, samples_keep = diversity_utils.get_SADs(samples)


SAD_idx_list = random.sample(range(len(SADs)), 30)


SADs_sample = [SADs[i] for i in SAD_idx_list]
taxonomy_names_sample = [taxonomy_names[i] for i in SAD_idx_list]
samples_keep = [samples_keep[i] for i in SAD_idx_list]


sys.stderr.write("Loading tree...\n")

tree = ete3.Tree('%sclosed_ref_silva/97_otus.tre' % config.data_directory, quoted_node_names=True, format=1)

sys.stderr.write("Calculating diversity and tree features...\n")

record_strs = [",".join(['SampleID', 'evenness', 'symmetry'])]

for i in range(len(SADs_sample)):

    SADs_sample_i = SADs_sample[i]
    taxonomy_names_sample_i = taxonomy_names_sample[i]
    samples_keep_i = samples_keep[i]

    tree_subset = tree_utils.subset_tree(taxonomy_names_sample_i, tree)

    symmetry = tree_utils.tree_symmetry(tree_subset).weighted_symmetry()

    evenness = diversity_utils.calculate_pielou_evenness(SADs_sample_i)

    record_str = ",".join([samples_keep_i, str(evenness), str(symmetry)])

    record_strs.append(record_str)



sys.stderr.write("Writing output file...\n")
filename = config.data_directory + 'test_tree_sad_shape.csv'
file = open(filename,"w")
record_str = "\n".join(record_strs)
file.write(record_str)
file.close()
sys.stderr.write("Done!\n")


# ete test

#t = Tree(directory+'test.tre', quoted_node_names=True, format=1)
#labels_sample = random.sample(tree_labels, 10)





#tree_labels = [t_i.name for t_i in t.get_leaves()]
#labels_sample = random.sample(tree_labels, 10)
