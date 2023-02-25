import sys
import dbd_utils



environment = str(sys.argv[1])
environment = environment.replace('_', ' ')


dbd_utils.run_svd_and_make_diversity_slm_dbd_simulation_phylo_dict(environment, iter_=500)