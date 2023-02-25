import sys
import dbd_utils



environment = str(sys.argv[1])
environment = environment.replace('_', ' ')


dbd_utils.make_diversity_slm_phlyo_integral_dict(environment)