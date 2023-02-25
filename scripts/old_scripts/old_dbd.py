

presence_absence_path = "%spresence_absence_annotated.pickle" % (config.data_directory)


def make_presence_absence_annotated_dict():

    sys.stderr.write("Subsetting samples...\n")

    samples = diversity_utils.subset_observations()

    sys.stderr.write("Getting SADs...\n")

    SADs, taxonomy_names, samples_keep = diversity_utils.get_SADs(samples)

    sys.stderr.write("Making presence-absence dictionary...\n")
    # D_5
    taxonomy_names_flat = numpy.concatenate(taxonomy_names).ravel()
    taxonomy_names_unique = list(set(taxonomy_names_flat.tolist()))

    afd_dict = {}
    taxonomy_dict = {}
    for line in open("%semp/otu_info/silva_123/taxonomy/consensus_taxonomy_7_levels.txt" % config.data_directory, 'r'):
        line_split = line.strip().split('\t')
        taxon = line_split[0]
        genus = line_split[1].split(';')[5].split('D_5__')[-1]

        if genus in diversity_utils.genera_to_ignore:
            continue

        if taxon in taxonomy_names_unique:
            #taxonomy_dict[taxon] = genus
            #afd_dict[taxon] = []
            afd_dict[taxon] = {}
            afd_dict[taxon]['genus'] = genus
            afd_dict[taxon]['presence_absence'] = []

        #all.append(taxon, genus)

    #taxa_to_keep = list(taxonomy_dict.keys())
    taxa_to_keep = list(afd_dict.keys())
    for i in range(len(SADs)):

        SAD_i = SADs[i]
        taxonomy_names_i = taxonomy_names[i]

        for t in taxa_to_keep:
            if t in taxonomy_names_i:
                afd_dict[t]['presence_absence'].append(1)
            else:
                afd_dict[t]['presence_absence'].append(0)

    sys.stderr.write("Saving presence-absence dictionary...\n")
    with open(presence_absence_path, 'wb') as handle:
        pickle.dump(afd_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




def load_presence_absence_annotated_dict():

    with open(presence_absence_path, 'rb') as handle:
        presence_absence_dict = pickle.load(handle)
    return presence_absence_dict






def make_presence_absence_slope_dict():

    sys.stderr.write("Loading presence-absence dict...\n")

    pres_abs_dict = load_presence_absence_annotated_dict()

    sys.stderr.write("Estimating dbd slope...\n")

    taxa = numpy.asarray(list(pres_abs_dict.keys()))
    all_genera = list(set([pres_abs_dict[t]['genus'] for t in taxa]))

    dbd_slope_dict = {}

    for focal_genus in all_genera:

        presence_absence_all = []
        n_asvs_focal_genus = 0
        n_asvs_non_focal_genus = []

        genus_to_taxa_dict = {}

        for t in taxa:
            genus_t = pres_abs_dict[t]['genus']

            if genus_t == focal_genus:
                n_asvs_focal_genus += 1

            else:
                n_asvs_non_focal_genus.append(genus_t)
                if genus_t not in genus_to_taxa_dict:
                    genus_to_taxa_dict[genus_t] = []
                genus_to_taxa_dict[genus_t].append(t)

            presence_absence_all.append(pres_abs_dict[t]['presence_absence'])

        # dont look at genera with fer species.
        if n_asvs_focal_genus < 10:
            continue

        unique_non_focal_genera, counts_non_focal_genera = numpy.unique(n_asvs_non_focal_genus, return_counts=True)
        coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_non_focal_genera))[:-1]
        presence_absence_array = numpy.asarray(presence_absence_all)

        # each element is the number of ASVS belonging to the non-focal genus in a given site
        richness_focal = numpy.zeros(presence_absence_array.shape[1])
        richness_non_focal = numpy.zeros(presence_absence_array.shape[1])

        for t in taxa:
            if pres_abs_dict[t]['genus'] == focal_genus:
                richness_focal += numpy.asarray(pres_abs_dict[t]['presence_absence'])


        for g in unique_non_focal_genera:

            g_taxa = genus_to_taxa_dict[g]
            g_taxa_idx = numpy.asarray([numpy.where(taxa == g_taxa_i)[0][0] for g_taxa_i in g_taxa])
            # get total number of ASVs belonging to a given genus
            g_sum = numpy.sum(presence_absence_array[g_taxa_idx,:], axis=0)
            # make presence-absence for a single genues
            g_sum[g_sum>0] = 1
            #richness_non_focal += numpy.sum(presence_absence_array[g_taxa_idx,:], axis=0)
            richness_non_focal += g_sum

        slope, intercept, r_value, p_value, std_err = stats.linregress(richness_non_focal, richness_focal)
        slope_null_all = []
        for i in range(iter):
            numpy.random.shuffle(presence_absence_array)
            # each element is the number of ASVS belonging to the non-focal genus in a given site
            presence_absence_focal_null = presence_absence_array[:n_asvs_focal_genus,:]
            presence_absence_non_focal_null = presence_absence_array[n_asvs_focal_genus:,:]
            richness_focal_null = numpy.sum(presence_absence_focal_null, axis=0)

            presence_absence_non_focal_coarse_null = numpy.add.reduceat(presence_absence_non_focal_null, coarse_grain_by_genus_idx, axis=0)
            # make presence-absence again
            presence_absence_non_focal_coarse_null[presence_absence_non_focal_coarse_null>0] = 1
            richness_non_focal_null = numpy.sum(presence_absence_non_focal_coarse_null, axis=0)

            slope_null, intercept_null, r_value_null, p_value_null, std_err_null = stats.linregress(richness_non_focal_null, richness_focal_null)
            slope_null_all.append(slope_null)

        slope_null_all = numpy.asarray(slope_null_all)
        slope_null_all = numpy.sort(slope_null_all)

        slope_null_median = numpy.median(slope_null_all)
        if slope >= slope_null_median:
            p_value = sum(slope_null_all > slope)/iter
        else:
            p_value = sum(slope_null_all < slope)/iter

        #p_value = sum(slope_null_all > slope)/iter
        lower_ci = slope_null_all[int(iter*0.025)]
        upper_ci = slope_null_all[int(iter*0.975)]

        dbd_slope_dict[focal_genus] = {}
        dbd_slope_dict[focal_genus]['slope'] = slope
        dbd_slope_dict[focal_genus]['null_dist'] = slope_null_all.tolist()
        dbd_slope_dict[focal_genus]['p_value'] = p_value
        dbd_slope_dict[focal_genus]['lower_ci'] = lower_ci
        dbd_slope_dict[focal_genus]['upper_ci'] = upper_ci
        dbd_slope_dict[focal_genus]['n_asvs_focal_genus'] = n_asvs_focal_genus
        dbd_slope_dict[focal_genus]['standard_score'] =  (slope - numpy.mean(slope_null_all)) / numpy.std(slope_null_all)


        print(focal_genus, slope, lower_ci, upper_ci, p_value)


    sys.stderr.write("Saving slope dictionary...\n")
    with open(slope_path, 'wb') as handle:
        pickle.dump(dbd_slope_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
