from __future__ import division
import string
import config
import random, sys
import pandas
import numpy
from ete3 import Tree
import ete3
from itertools import combinations


random.seed(123456789)



sys.setrecursionlimit(15000)

#tree = get_emp_tree()
#print(max(calculate_phylogenetic_distance_all(tree)))
# returns = 3.410779999999999


def subset_tree(species_list, tree):

    # deep copy. slowest copy option, can change later if this is too slow
    tree_copy = tree.copy()

    #tree_labels = [t_i.name for t_i in tree_copy.get_leaves()]

    tree_copy.prune(species_list, preserve_branch_length=True)


    return tree_copy


def calculate_edge_length_abundance_distribution(tree):

    def search_by_size(tree, size):
        "Finds nodes with a given number of leaves"
        matches = []
        for n in tree.traverse():
           if len(n) == size:
              matches.append(n)
        return matches

    number_tips = len([tip for tip in tree])

    ead_range = list( range(1, number_tips))
    ead = []

    for i in ead_range:

        edges = search_by_size(tree, size=i)

        edge_lengths = [edge.dist for edge in edges]

        edge_sum = sum(edge_lengths)

        ead.append(edge_sum)

    ead_range = numpy.asarray(ead_range)
    ead = numpy.asarray(ead)

    return ead_range, ead



def subsample_phylogenetic_diversity(tree, species_abundance_distribution, sample_size):

    # assumes that tree has already been subset
    # phylogenetic diversity, eqn. 1, O'Dwyer et al., 2015

    N = sum(species_abundance_distribution)

    ead_range, ead = calculate_edge_length_abundance_distribution(tree)

    phylogenetic_diversity = sum([ ead[i] * (1 - (((1 - sample_size) / N)**k)) for i in range(ead_range) ])

    return phylogenetic_diversity




def cache_distances(tree):

    ''' precalculate distances of all nodes to the root'''

    node2rootdist = {tree:0}
    for node in tree.iter_descendants('preorder'):
        node2rootdist[node] = node.dist + node2rootdist[node.up]
    return node2rootdist



def collapse(tree_copy, node2tips, root_distance,  min_dist):

    # collapses all nodes within a given distance for the tree
    # cache the tip content of each node to reduce the number of times the tree_copy is traversed

    for node in tree_copy.get_descendants('preorder'):
        if not node.is_leaf():
            avg_distance_to_tips = numpy.mean([root_distance[tip]-root_distance[node]
                                         for tip in node2tips[node]])

            if avg_distance_to_tips < min_dist:
                # do whatever, ete support node annotation, deletion, labeling, etc.
                # rename
                # ' COLLAPSED'
                node.name += 'COLLAPSED-%g-%s' %(avg_distance_to_tips,
                                                 ','.join([tip.name for tip in node2tips[node]]))
                # label
                node.add_features(collapsed=True)
                # set drawing attribute so they look collapsed when displayed with tree_copy.show()
                #node.img_style['draw_descendants'] = False

    return tree_copy



def get_emp_tree():
    tree = ete3.Tree('%semp/otu_info/silva_123/97_otus.tre' % config.data_directory, quoted_node_names=True, format=1)
    return tree



def calculate_phylogenetic_distance_all(tree):

    tree_copy = tree.copy()

    root_distance = cache_distances(tree_copy)

    distance_all = list(root_distance.values())

    return distance_all



def get_dalbello_tree():
    tree = ete3.Tree('%sdalbello/seqtab-nochim_muscle_fasttree.tre' % config.data_directory, quoted_node_names=True, format=1)
    return tree




def make_distance_collapsed_dict(tree, distances, min_n_clades=5):

    distance_collapsed_dict = {}

    n_tips = len(tree.get_leaf_names())

    for distance in distances:

        tree_copy = tree.copy()
        node2tips = tree_copy.get_cached_content()
        root_distance = cache_distances(tree_copy)
        tree_coarse_grained = collapse(tree_copy, node2tips, root_distance, distance)

        # collapsed nodes are labeled, so you locate them and prune them
        for n in tree_coarse_grained.search_nodes(collapsed=True):
            for descendents in n.get_children():
                descendents.detach()

        n_total = 0
        collapsed_node_names_all = []
        for t in tree_coarse_grained.get_leaf_names():

            if 'COLLAPSED-' in t:
                collapsed_node_names = t.split('-')[-1].split(',')

            else:
                collapsed_node_names = [t]

            collapsed_node_names_all.append(collapsed_node_names)
            n_total+= len(collapsed_node_names)



        if (len(collapsed_node_names_all) < min_n_clades) or (len(collapsed_node_names_all) == n_tips):
            continue

        print(distance, len(collapsed_node_names_all))

        distance_collapsed_dict[distance] = collapsed_node_names_all



    return distance_collapsed_dict






class tree_symmetry():

    '''
    Shape stastic for tree symmetry proposed in Siegel and Sugihara 1983
    Used in Sugihara et al., 2003
    '''

    def __init__(self, tree):

        self.tree = tree

    ## Create random string (for naming nodes in min/max trees)
    def rand_str(self, n):
        return(''.join([random.choice(string.ascii_lowercase) for i in range(n)]))


    ## Maximum symmetry tree w/ num_leaves number of leaves
    def max_sym_tree(self, num_leaves):
        t1 = Tree()

        N_log2 = int(numpy.log2(num_leaves))
        mod = num_leaves % (2**N_log2)
        for i in range(N_log2):
            for node in t1.get_leaves():
                node.populate(2,names_library=[self.rand_str(10),self.rand_str(10)])

        chosen = numpy.linspace(0,(2**N_log2)-1,mod)
        nodes = [t1.get_leaves()[int(i)] for i in chosen]
        for node in nodes:
            node.populate(2,names_library=[self.rand_str(10),self.rand_str(10)])

        return t1


    ## Minimum symmetry tree w/ num_leaves number of leaves
    def min_sym_tree(self, num_leaves):

        t = Tree()
        num_leaves = num_leaves
        N = 2*num_leaves
        node1 = t

        for i in range(N):
            node1 = node1.add_child(name=self.rand_str(10))
            node2 = node1.add_child(name=self.rand_str(10))

        return t


    ## unnormalized symmetry
    def raw_symmetry(self, tree):

        all_nodes = [node.name for node in tree.get_leaves()]
        df_tree = pandas.DataFrame(columns=all_nodes,index=all_nodes)

        node_pairs = list(combinations(all_nodes, r = 2))

        distance_dict = {}

        E = []

        for node_pair in node_pairs:
            node1 = node_pair[0]
            node2 = node_pair[1]
            A = tree & node1
            C = tree & node2
            P = tree.get_common_ancestor(A,C)
            # so we only consider the distance of one node to the common ancestor?
            d_P_A = A.get_distance(P, topology_only=True)
            d_P_C = C.get_distance(P, topology_only=True)

            #distance_dict[node_pair] = (d_P_A, d_P_C)

            E.append(d_P_A + d_P_A*d_P_C - d_P_A**2)
            E.append(d_P_C + d_P_A*d_P_C - d_P_C**2)

        #for node1 in all_nodes:
        #    print(node1)
        #    for node2 in all_nodes:
        #        A = tree & node1
        #        C = tree & node2
        #        P = tree.get_common_ancestor(A,C)
        #        d = A.get_distance(P, topology_only=True)
        #        df_tree[node1][node2] = d

        #for node1 in df_tree.index:
        #    for node2 in df_tree.columns[df_tree.columns != node1]:
        #        N1 = df_tree[node1][node2]
        #        N2 = df_tree[node2][node1]
        #        E.append(N1 + N1*N2 - N1**2)

        return numpy.mean(E)

    ## normalized symmetry
    def weighted_symmetry(self):

        n = len(self.tree.get_leaves())

        min_tree = self.min_sym_tree(n)
        max_tree = self.max_sym_tree(n)

        min_symmetry = self.raw_symmetry(min_tree)
        max_symmetry = self.raw_symmetry(max_tree)

        tree_symmetry = self.raw_symmetry(self.tree)

        return (tree_symmetry-min_symmetry)/(max_symmetry - min_symmetry)
