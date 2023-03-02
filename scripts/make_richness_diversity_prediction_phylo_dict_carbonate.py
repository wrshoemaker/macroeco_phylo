import sys
import dbd_utils



environment = str(sys.argv[1])
environment = environment.replace('_', ' ')


dbd_utils.make_richness_diversity_prediction_phylo_dict(environment, iter_=100)