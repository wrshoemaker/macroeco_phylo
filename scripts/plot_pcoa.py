import utils
import skbio



for c in utils.carbon_sources:

    s_by_s, asv, community, inocula, od600 = utils.get_s_by_s(c, 'M9', remove_inocula_pres_once=True)
