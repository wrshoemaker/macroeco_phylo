
import sys
import diversity_utils
import dbd_utils
import numpy

environments_to_keep = diversity_utils.environments_to_keep
rarefied = False

for environment in environments_to_keep:

    #coarse_grained_tree_dict = coarse_grained_tree_dict_all[environment]
    #coarse_grained_tree_dict = dbd_utils.load_coarse_grained_tree_no_subsampling_dict(environment=environment, rarefied=rarefied)

    #sys.stderr.write("Subsetting samples for %s...\n" % environment)
    pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_phylo_dict(environment, rarefied = False)
    samples = pres_abs_dict['samples']
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
    
    sys.stderr.write("Getting site-by-species matrix...\n")
    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])


    n_reads = numpy.sum(s_by_s, axis=0)
    n_otus = len(taxa)
    #print('Phylogenetic', environment)
    #print(n_otus)
    #print(min(n_reads.tolist()))


    # print for taxonomic
    pres_abs_dict = dbd_utils.load_sad_annotated_no_subsampling_taxon_dict(environment, rarefied = False)
    samples = pres_abs_dict['samples']
    taxa = numpy.asarray(list(pres_abs_dict['taxa'].keys()))
    
    sys.stderr.write("Getting site-by-species matrix...\n")
    s_by_s = numpy.asarray([pres_abs_dict['taxa'][t]['abundance'] for t in taxa])


    n_reads = numpy.sum(s_by_s, axis=0)
    n_otus = len(taxa)

    #print('Taxonomic', environment)
    #print(n_otus)
    print(min(n_reads.tolist()))