from __future__ import division
import config
import biom
import numpy
import os
os.environ["SYMPY_USE_CACHE"]="no"
import sympy
#import math
import simulation_utils
import scipy.integrate as integrate
import scipy.special as special
#import dbd_utils


#SampleID       BarcodeSequence LinkerPrimerSequence    Description     host_subject_id study_id        title   principal_investigator  doi     ebi_accession   target_gene     target_subfragment      pcr_primers     illumina_technology     extraction_center       run_center      run_date        read_length_bp  sequences_split_libraries       observations_closed_ref_greengenes      observations_closed_ref_silva   observations_open_ref_greengenes        observations_deblur_90bp
#observations_deblur_100bp       observations_deblur_150bp       emp_release1    qc_filtered     subset_10k      subset_5k       subset_2k       sample_taxid    sample_scientific_name  host_taxid      host_common_name_provided       host_common_name        host_scientific_name    host_superkingdom       host_kingdom    host_phylum     host_class      host_order      host_family     host_genus      host_species    collection_timestamp    country latitude_deg    longitude_deg   depth_m
#altitude_m      elevation_m     env_biome       env_feature     env_material    envo_biome_0    envo_biome_1    envo_biome_2    envo_biome_3    envo_biome_4    envo_biome_5    empo_0  empo_1  empo_2  empo_3  adiv_observed_otus      adiv_chao1      adiv_shannon    adiv_faith_pd   temperature_deg_c       ph      salinity_psu    oxygen_mg_per_l phosphate_umol_per_l    ammonium_umol_per_l     nitrate_umol_per_l      sulfate_umol_per_l


#taxa_ranks = ['phylum', 'class', 'order', 'family', 'genus']
#taxa_ranks_label = ['Phylum', 'Class', 'Order', 'Family', 'Genus']

taxa_ranks = ['genus', 'family', 'order',  'class',  'phylum']
taxa_ranks_label = ['Genus', 'Family', 'Order', 'Class', 'Phylum']


taxa_ranks_with_asv = ['OTU', 'genus', 'family', 'order',  'class',  'phylum']
taxa_ranks_label_with_asv = ['OTU', 'Genus', 'Family', 'Order', 'Class', 'Phylum']


taxa_ranks_label_with_asv = ['OTU', 'Genus', 'Family', 'Order', 'Class', 'Phylum']

taxa_ranks_label_capital_dict = {'OTU': 'OTU', 'genus':'Genus', 'family':'Family', 'order':'Order', 'class':'Class', 'phylum':'Phylum'}



def get_rarefied_label(rarefied):

    if rarefied == True:
        rarefied_label = '_rarefied'
    else:
        rarefied_label = ''

    return rarefied_label


def get_environment_label(environment):

    environment_label = environment.replace(' ', '_')

    return environment_label



def format_environment_label(environment):

    environment_label = ' '.join(environment.split(' ')[:-1]).capitalize()

    return environment_label


def get_label(environment, rarefied):

    label =  get_environment_label(environment) + get_rarefied_label(rarefied)

    return label


#environments_to_keep = ['root metagenome', 'air metagenome', 'freshwater metagenome', 'coral metagenome',\
#                        'fish metagenome', 'rhizosphere metagenome', 'bovine gut metagenome', 'marine sediment metagenome',\
#                        'marine metagenome', 'soil metagenome', 'human oral metagenome', 'human skin metagenome',\
#                        'primate metagenome', 'microbial mat metagenome', 'freshwater sediment metagenome', 'human gut metagenome']


environments_to_keep = ['freshwater sediment metagenome', 'freshwater metagenome', 'marine sediment metagenome',\
                        'marine metagenome', 'human oral metagenome', 'human skin metagenome',\
                        'soil metagenome', 'human gut metagenome', 'microbial mat metagenome']



genera_to_ignore = ['Ambiguous_taxa', 'uncultured', 'uncultured marine archaeon', \
                    'bacterium enrichment culture clone E27', 'uncultured Stenotrophomonas sp.', \
                    "bacterium 'New Zealand C'", 'uncultured bacterium', 'uncultured archaeon', \
                    'Ambiguous_taxa', 'uncultured organism', 'uncultured Verrucomicrobia bacterium', \
                    'uncultured Acidobacteria bacterium', 'uncultured rumen bacterium', 'uncultured actinobacterium', \
                    'uncultured candidate division WS6 bacterium', 'uncultured Pseudonocardiaceae bacterium', \
                    'uncultured Chloroflexi bacterium', 'uncultured alpha proteobacterium', 'Rubus hybrid cultivar', \
                    'Arachis hypogaea (peanut)', 'I-8', 'archaeon GW2011_AR16', 'uncultured rumen methanogen', \
                    'uncultured Crater Lake bacterium CL500-15', 'uncultured soil bacterium', 'uncultured euryarchaeote', \
                    'uncultured planctomycete', 'uncultured bacterium #0319-7F19', 'uncultured bacterium adhufec279',\
                    'Phyllocladus trichomanoides (celery pine)', 'uncultured Bacilli bacterium', 'uncultured Cytophagales bacterium',\
                    'uncultured soil bacterium PBS-II-5', 'gamma proteobacterium endosymbiont of Astomonema sp.', 'Musa acuminata subsp. malaccensis',\
                    'Hordeum vulgare subsp. vulgare (domesticated barley)', 'Fragaria vesca subsp. vesca', 'Bremia lactucae (lettuce downy mildew)',\
                    'uncultured bacterium #0319-6C24', 'uncultured crenarchaeote pBRKC125', 'SV1-3', 'endosymbiont of Pogonophora sp. JT-1',\
                    'bacterium enrichment culture clone JCA1', 'uncultured bacterium UASB_TL56', 'unclassified Pseudomonadales (miscellaneous)',\
                    'uncultured archaeon WCHA1-57', 'Iris sp. Qiu 95091', 'Phacus sp. 5 JIK-2013', 'uncultured Acidimicrobiaceae bacterium',\
                    'uncultured soil bacterium PBS-II-1', 'uncultured Sphingobacteriales bacterium', 'uncultured gamma proteobacterium HF4000_19M20',\
                    'Silene noctiflora (night-flowering catchfly)', 'uncultured sludge bacterium A21b', 'bacterium enrichment culture clone auto8_4W',\
                    'uncultured delta proteobacterium', 'uncultured gamma proteobacterium CHAB-XI-27', 'uncultured Gemmatimonadetes bacterium',\
                    'uncultured Microgenomates bacterium', 'uncultured bacterium TA08', 'uncultured delta proteobacterium', 'uncultured compost bacterium',\
                    'uncultured Bacteroidetes bacterium', 'uncultured Endomicrobia bacterium', 'uncultured deep-sea bacterium', 'uncultured Spartobacteria bacterium',\
                    'uncultured gamma proteobacterium', 'uncultured Sinobacteraceae bacterium', 'uncultured Rhodocyclaceae bacterium', 'uncultured crenarchaeote',\
                    'uncultured marine bacterium', 'uncultured forest soil bacterium', 'uncultured Parcubacteria bacterium', 'uncultured methanogenic archaeon',\
                    'uncultured Nitrospirae bacterium', 'uncultured archaeon APA1-0cm', 'uncultured Candidatus Saccharibacteria bacterium', 'uncultured bacterium gp10',\
                    'uncultured Myxococcales bacterium', 'uncultured microorganism', 'uncultured Latescibacteria bacterium', 'uncultured prokaryote', 'uncultured thaumarchaeote',\
                    'uncultured verrucomicrobium DEV05', 'uncultured Sphingomonadaceae bacterium', 'uncultured sludge bacterium', 'uncultured Gram-negative bacterium', \
                    'uncultured Acetobacteraceae bacterium', 'uncultured Acidimicrobidae bacterium', 'uncultured diatom', 'uncultured cyanobacterium', 'uncultured Kofleriaceae bacterium',\
                    'uncultured epsilon proteobacterium', 'uncultured SAR11 cluster alpha proteobacterium', 'uncultured anaerobic bacterium', 'uncultured Burkholderiaceae bacterium',\
                    'uncultured proteobacterium', 'uncultured Chlorobi bacterium', 'uncultured Rhizobiales bacterium', 'uncultured eukaryote', 'uncultured haloarchaeon', 'uncultured phototrophic eukaryote',\
                    'uncultured Thermoprotei archaeon', 'uncultured Acidobacteriales bacterium', 'uncultured Firmicutes bacterium', 'uncultured sediment bacterium', 'uncultured soil bacterium PBS-25',\
                    'uncultured Clostridia bacterium', 'uncultured Green Bay ferromanganous micronodule bacterium MND8', 'uncultured Firmicutes bacterium', 'uncultured low G+C Gram-positive bacterium',\
                    'uncultured marine group II euryarchaeote', 'bacterium Ellin6530', 'Ricinus communis (castor bean)', 'Termite Treponema cluster', 'Triticum aestivum (bread wheat)', \
                    'Ipomoea batatas (sweet potato)', 'Nicotiana tabacum/Hyoscyamus niger cybrid', 'Pseudotsuga menziesii (Douglas-fir)', 'castanea mollissima (chinese chestnut)',\
                    'myristica fragrans (nutmeg)', 'medicago truncatula (barrel medic)', 'vitis vinifera (wine grape)', 'punica granatum (pomegranate)', 'proteobacterium ls-2',\
                    'mwh-ta3', 'rs-h88 termite group', 'eubostrichus topiarius associated bacterium 1', 'proteobacterium ls-2', 'rhizanthella gardneri',  'candidatus methanomethylophilus',\
                    'possible genus 05', 'proteobacterium fa350', 'parcubacteria bacterium raac4_od1_1', 'pir4 lineage', 'trachelomonas hispida', 'candidatus evansia',\
                    'tetraselmis cordiformis', 'bacterium whc7-12', 'alpha proteobacterium 047', 'om43 clade', 'prevotellaceae yab2003 group', 'candidatus phytoplasma', 'Firmicutes bacterium CAG:822',\
                    'X', 'x']


# '(' these are the host spececies names
# domain level taxonomy checked with
#https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi


phylum_partial_match_to_ignore = []
class_partial_match_to_ignore = ['south african gold mine gp 1(sagmcg-1)', 'oscillatoria sp. ntgm13', 'unknown class', 'uncultured']
order_partial_match_to_ignore = ['order iii', 'uncultured', 'km16', 'botryococcus braunii', 'chakia ciliosa 27', 'subgroup 3', 'bdallophytum americanum', 'subsectionv', 'unknown order', 'subsectioniv',\
                                'nitzschia sp. iriss06', 'hydrogenispora ethanolica', 'stigeoclonium helveticum', 'bathycoccus prasinos', 'nitzschia sp. iriil01', 'weingartia kargliana', 'order iv',\
                                'cytophaga sp. prpr22', 'order ii', 'virgulinella fragilis', 'cytophaga sp. dex80-43', 'aquamonas haywardensis', '4-15', 'nemalionopsis tortuosa', 'subsectioniii', 'subgroup 4',\
                                'aaa34a10', 'incertae sedis', 'subsectionii', 'dsm 5130', 'b38', 'loriellopsis cavernicola lf-b5']


family_partial_match_to_ignore = ['uncultured', 'incertae sedis', 'family', 'cluster', ' sp. ', 'mn 122.2a', 'ma-28-i98c', '64k2', 'lwsr-14', '1174-901-12', 'orca-3n101', 'gom arc i', 'ferrovum', 'mitochondria',\
                                    'c2u', 'surface 1', 'blfdi19', 'akau3744', 'mob164', 'p-102', 'wchb1-69', 'tbz33', 'kcm-b-112']



genera_partial_match_to_ignore = ['uncultured', 'unidentified', 'peperomia', 'enrichment culture clone', 'symbiont', '(', 'chlamydomonas', 'trachelomonas', \
                                'discoplastis', 'phacus', 'didymeles', 'pyramimonas', 'maytenus', 'filamentous cyanobacterium odo1mo59', 'marine metagenome',\
                                'plantago', 'seep-srb1', 'archaeon gw2011_ar20', 'blvii28 wastewater-sludge group', 'dictyochloropsis', 'acholeplasmatales bacterium canine oral taxon 375',
                                'euglenaria', 'actinobacterium yim 75507', 'termite', 'kartchner caverns bacterium pf-h', 'clostridium sensu stricto 3', 'stellarima',\
                                'pir2 lineage', 'glarea lozoyensis', 'u29-b03', 'hotseep-1', 'bacterium str. 77003', 'alpha proteobacterium scgc aaa288-e13', 'euglena',\
                                'saltmarsh clone lcp-67', 'family xiii ucg-002', 'phormidiaceae', 'alpha proteobacterium', 'fs140-16b-02 marine group', 'om75 clade', 'ruminococcaceae',\
                                'lepocinclis', 'bacterium str. 77003', 'possible genus sk003-sk004', 'planctonema', 'arachis', 'acidobacteria bacterium wx27', 'family xiii ucg-001',\
                                'hydnora', 'metagenome', 'p-1088-a5 gut group', 'gelidium', 'w4', 'hypericum', 'planctomycete', 'hoa5-07d05 gut group', 'watanabea', 'gal15',\
                                'sclerotinia', 'galdieria', 'idiospermum', 'heliobacteriaceae bacterium', 'agricultural soil bacterium sc-i-84', 'xanthomonadaceae bacterium wwh73',\
                                'sva0081 sediment group', 'clostridiales bacterium canine oral taxon 260', 'nicotiana', 'clostridium sensu stricto 2', 'lachnospiraceae ucg-003',\
                                'jakoba', 'rice cluster', 'lycopodium', 'sciadopitys', 'aquilaria', 'reticulomyxa', 'xylochloris', 'plasmodium', 'tsukubamonas', 'gracilariopsis',\
                                'oophila', 'ps-b30', 'mwh-ta3', 'aenigmarchaeota archaeon jgi 0000106-f11', 'sm1a02', 'nephroselmis', 'cb31g03', 'methanobacteriaceae archaeon 15az',\
                                'oryza', 'epixenosomes of euplotidium arenarium', 'bacterium 09.96.18', 'gammaproteobacteria bacterium q1', 'clostridium sensu stricto 10', ' clone ',\
                                'houttuynia', 'epipogium', 'sensu stricto', 'tetraselmis', 'planctomycete', 'bacterium gla1', 'vaccinium', 'candidate division tm7 bacterium ly2',\
                                'prasinophyceae', 'eubostrichus topiarius associated bacterium 1', 'lachnospiraceae ucg-010', 'hymenophyton', 'eimeria', 'phalacroma', 'polycephalomyces',\
                                'karenia', 'gossypium', 'ettlia', 'methanogenic archaeon ch1270', 'aristolochia', 'flavobacteriaceae bacterium mola 32', 'p131-4', 'soil bacterium wwh121',\
                                'solanum', 'monomastix', 'peumus', 'nannochloropsis', 'parcubacteria', 'em3', 'sp3-e08', 'prasinoderma', 'bacterium wh3-1', 'cryptoglena', 'mucus bacterium 80',\
                                'elodea', 'bacterium em-19', 'archaeon gw2011_ar5', 'arthroderma', 'corynaea', 'lathyrus', 'cordyceps', 'lachnospiraceae ucg-009', 'isaria', 'verrucomicrobia bacterium lx181',\
                                'oedogonium', 'cylindrospermosis', 'cylindrotheca', 'seculamonas', 'bacterium whc4-8', 'neosiphonia', 'picochlorum', 'escherichia-shigella', 'bacillaceae bacterium efn-4',\
                                'flavobacteriaceae bacterium dokdo 020', 'bacteroidetes bacterium', 'babesia', 'pterocladiella', 'erysipelotrichaceae ucg-006', 'rhizobiales bacterium ume16',\
                                'lachnospiraceae nk3a20 group', 'amborella', 'withania', 'lachnospiraceae ucg-004', 'archaeon gw2011_ar17', 'galeola', 'actinobacterium msi70', 'ml602j-51',\
                                '12up', 'wildemania', 'colletotrichum', 'hylocomium', 'prevotellaceae ucg-004', 'wastewater-sludge group', 'chlorarachnion', 'euptilota', 'azolla', 'lotus ',\
                                'bd1-7 clade', 'archaeon gw2011_ar15', 'austrobaileya', 'passiflora', 'prevotellaceae ucg-001', 'costus', 'coccomyxa', 'incertae', 'gks98 freshwater group',\
                                'scaevola', 'aemula', 'salvia', 'eutreptia', 'karlodinium', 'desmochloris', 'lachnospiraceae nd3007 group', 'chaetoceros', 'chroomonas', 'prototheca', 'orchidantha', \
                                'codium', 'ptilocladia', 'microgenomates bacterium scgc aaa011-a19', 'deltaproteobacteria bacterium canine oral taxon 266', 'aureoumbra', 'betula', 'bdellocephala',\
                                'rhodobacteraceae bacterium b62ydz-zz', 'mesotaenium', 'gracilibacteria bacterium canine oral taxon 364', 'probable genus ', 'dga-11 gut group', 'actinobacterium gp-6',\
                                'guillardia', 'buxbaumia', 'corethron', 'defluviitaleaceae ', 'actinobacterium ', 'palmaria', 'microgenomates bacterium ', 'bacterium ellin6529', 'cl500-29 marine group',\
                                'diplosphaera', 'chlamydiales bacterium ', 'crib 32', 'attheya', 'gleichenia', 'strombomonas', 'bacterium whc7-12', 'hypnum', 'isoetes', 'glacial ice ', 'najas',\
                                ' marine group', 'rhipidosiphon', 'hgci clade', 'heveochlorella', 'acetobacteraceae bacterium ', 'wx59', 'noccaea', 'elaeis', 'erysipelotrichaceae', 'pleea',\
                                'microthamnion', 'gracilaria', 'coriobacteriaceae ', 'eutreptiella', 'athetis', 'sporothrix', 'silene', 'marine ', 'braarudosphaera', 'thermoactinomycetaceae bacterium',\
                                'gracilariophila', 'thecamonas', 'lachnospiraceae ucg-007', 'porphyra', 'possible genus ', 'acidobacteria bacterium ', 'desulfovibrionales bacterium ',\
                                'ds001', 'proteobacterium ', 'aeginetia', 'entransia', 'capsicum', 'rhodobacteraceae bacterium ', 'bryopsis', 'cyanobacterium ', 'porites', 'extubocellulus',\
                                'bacillales bacterium ', 'mi4', 'sporolithon', 'piper betle', 'proteobacteria bacterium ', 'dehalococcoidia bacterium ', ' bacterium ', ' group', 'cn77', 'p1rh1',\
                                'anme-3', 'jl-etnp-f27', 'nitella', 'msbl7', 'ls-2', 'cryptomonas', 'ptilocladiopsis', 'ephemera', 'thottea', 'chorispora', 'pseudendoclonium', 'c1-b045', ' clade ',\
                                'scherffelia', 'nuphar', 'actinidia', 'boechera', 'papuacedrus', 'pb79', ' str. ', ' clade', 'colobanthus', 'seep-srb2', 'resultomonas', 'flabellia', 'bacterium ellin6504',\
                                'scgc ab-539-j10', 'panax', 'lobochlamys', 'porphyridium', 'lolium', 'cryptocarya', 'utricularia', 'sc103', 'odontella', 'glaucocystis', 'bacterium id4391', 'pythium',\
                                'aegilops', 'colacium', 'euryarchaeota archaeon scgc aaa286-e23', 'bacterium ellin6543', 'pt46', 'centroceras', 'abies', 'cyanidium', 'halobacteriaceae archaeon ',\
                                'bacterium ellin517', 'bacterium ellin6519', 'bacterium ', 'crenarchaeote ', 'gleocapsa', 'pabia', 'rhipocephalus', 'monoraphidium', 'rhipilia', 'pyropia', 'lachnospiraceae ucg-002',\
                                'crouania', 'ecdeiocolea', 'talaromyces', 'tetraphis', 'prasiolopsis', 'centrolepis', 'metzgeria', 'cucurbita', 'bangiopsis', 'rebecca', 'camelina', 'spirogyra',\
                                'takakia', 'polytoma', 'gloeochaete', 'oltmannsiellopsis', 'akyg587', 'lachnospiraceae ucg-005', 'trebouxiophyceae', 'curcuma', 'pilostyles', 'asterionella', 'cymbomonas', 'hypseocharis',\
                                'cynomorium', 'succinivibrionaceae ', 'conopholis', 'trichocolea', 'echinogammarus', 'mesostigma', 'hyalella', 'spermatozopsis', 'pir4 lineage', 'cassytha', 'spi55',\
                                'thorea', 'anme-2b', 'coriobacteriaceae ', 'ucg-003', 'archaeon gw2011_ar13', 'z195mb87', 'lachnospiraceae ', 'cl500-3', 'syngonanthus', 'wy67', 'verticillium', 'hypnea',\
                                'frullania', 'zygnema', 'mamiella', 'brasenia', 'pir1 lineage', 'chlorella', 'tdx16', 'heteroptilon', 'humulus', 'eucalyptus', 'triglochin', 'liriodendron', 'closteriopsis',\
                                'juniperus', 'flintiella', 'malawimonas', 'grateloupia', 'parachlorella', 'phaeoceros', 'fucus', 'prd01a011b', 'pedinomonas', 'chloroidium', 'polytrichum', 'rhodochaete',\
                                'helosis', 'reclinomonas', 'cheilotheca', 'haemonchus', 'potamogeton', 'dixoniella', 'nanoarchaeota archaeon scgc aaa011-d5', 'genlisea', 'trimenia', 'ahnfeltia', 'archaeon gw2011_ar18',\
                                'georgeantha', 'nonionella', 'laurus', 'chaetosphaeridium', 'sphagnum', 'carneum', 'pseudochloris', 'phaeodactylum', 'exanthemachrysis', 'phyllocladus', 'plocamium', 'asahi brw2',\
                                'marsupiomonas', 'zymoseptoria', 'chrysodidymus', 'beauveria', 'koliella', 'nelumbo', 'fritillaria', 'tomentella', 'asclepias', 'symphyogyna', 'glaucosphaera', 'carpothamnion',\
                                'mitrastemon', 'zea mays', 'coleochaete', 'trifolium', 'bonnemaisonia', 'pterocladia', 'lucida', 'turnera', 'erythrotrichia', 'orontium', 'megaceros', 'dicloster', 'monomorphina',\
                                'juncus', 'persicaria', 'hedyosmum', 'compsopogon', 'gulsonia', 'gerbera', 'sphagnicola', 'andalucia', 'pterosperma', 'dendrobium', 'aspergillus', 'floydiella', 'archaeon gw2011_ar11',\
                                'puccinia', 'ophiostoma', 'dianthus', 'aglaothamnion', 'chlorella', 'closterium', 'caulerpa ', 'berberis', 'meliosma', 'leptographium', 'cytinus', 'cymbella', 'morus ', 'pir3 lineage',\
                                'polypodium', 'vertebrata', 'rhodella', 'prevotellaceae ucg-003', 'jasminum', 'planctomycete', 'carludovica', 'seep-srb4', 'avrainvillea', 'fusarium', 'ziziphus', 'chlorokybus',\
                                'ceramium', 'gymnodinium', 'authm297', 'anemia ', 'polysphondylium', 'xylosandrus', 'xanthorhiza', 'selaginella', 'helicosporidium', 'cyanoptyche', 'pachysandra', 'lecanicillium',\
                                'monsonia', 'haloarchaeon', 'paracoccidioides', 'histiona', 'chlorosarcina', 'anaerolineaceae', 'neocystis', 'grevillea', 'choreocolax', 'penicillium', 'lagarostrobos', 'rhodymenia',\
                                'marvania', 'crustomastix', 'annulohypoxylon', 'acremonium', 'botrychium', 'parathyasira', 'kappaphycus', 'acrosiphonia', 'siparuna', 'anemopsis', 'anurida', 'monotropastrum', 'cunninghamia',\
                                'catanema', 'leptosira', 'veillonellaceae', 'ptilophora', 'gelidiella', 'dolichomastix', 'pluvialis', 'leyanella', 'delphineis', 'zantedeschia', 'tigriopus', 'lactoris', 'sargentodoxa',\
                                'ostreobium', 'portulaca', 'cercis', 'hafnia', 'bacillaria', 'illicium', 'roya ', 'pandanus', 'mixed culture isolate', 'triticum', 'spongospora', 'asterionellopsis', 'psammodictyon',\
                                'baetis', 'fusochloris', 'glycine', 'microcachrys ', 'rhizanthella ', 'azadirachta ', 'podophyllum ', 'cuscuta ', 'schizomeris', 'cyperus ', 'paralemanea ', 'dunaliella ', 'geranium ',\
                                'argulus ', 'cyanophora ', 'hydrodictyon']




def all_observations():

    sample_IDs = []

    mapping_file_path = '%semp/mapping_files/emp_qiime_mapping_qc_filtered.tsv' % config.data_directory
    #mapping_file_path = '%semp/mapping_files/emp_qiime_mapping_subset_2k.tsv' % config.data_directory
    #sample_ids_all = []
    for line in open(mapping_file_path, 'r'):
        line = line.strip().split('\t')

        #sample_ids_all.append(line[31])

        #if 'human gut metagenome' in line[31]:
        sample_IDs.append(line[0])

    return sample_IDs




def subset_observations(environment='human gut metagenome'):

    sample_IDs = []

    mapping_file_path = '%semp/mapping_files/emp_qiime_mapping_qc_filtered.tsv' % config.data_directory
    #mapping_file_path = '%semp/mapping_files/emp_qiime_mapping_subset_2k.tsv' % config.data_directory
    sample_ids_all = []
    for line in open(mapping_file_path, 'r'):
        line = line.strip().split('\t')

        sample_ids_all.append(line[31])

        if environment == 'all':
            sample_IDs.append(line[0])

        else:
            if environment in line[31]:
                sample_IDs.append(line[0])

    return sample_IDs





def get_sample_to_environment_dict():

    dict_ = {}
    mapping_file_path = '%semp/mapping_files/emp_qiime_mapping_qc_filtered.tsv' % config.data_directory
    sample_ids_all = []
    for line in open(mapping_file_path, 'r'):
        line = line.strip().split('\t')

        sample_id = line[0]
        environment = line[31]

        dict_[sample_id] = environment


    return dict_




def get_SADs(samples, S_min = 10, S_max=500, rarefied=True):
    # returns list of SADs and list of taxon names for a list of EMP observations

    #biom_path = '%semp/otu_tables/closed_ref_silva/emp_cr_silva_16S_123.qc_filtered.biom' % config.data_directory
    #biom_path = '%semp/otu_tables/closed_ref_silva/emp_cr_silva_16S_123.subset_2k.rare_10000.biom' % config.data_directory
    if rarefied == True:
        biom_path = '%semp/otu_tables/closed_ref_silva/emp_cr_silva_16S_123.qc_filtered.rare_10000.biom' % config.data_directory

    else:
        biom_path = '%semp/otu_tables/closed_ref_silva/emp_cr_silva_16S_123.qc_filtered.biom' % config.data_directory


    biom_table = biom.load_table(biom_path)

    taxonomy = biom_table.metadata_to_dataframe('observation')

    SADs = []
    taxonomy_names = []

    # remove few samples that aren't in biom file
    samples_in_biom = [sample for sample in samples if sample in  biom_table.ids() ]

    biom_table_subset = biom_table.filter(samples_in_biom, axis='sample', inplace=False)
    #dense=True
    biom_table_subset_df = biom_table_subset.to_dataframe(dense=True)
    biom_table_subset_values = biom_table_subset_df.values
    taxa_labels = biom_table_subset_df.index.values
    biom_samples = biom_table_subset_df.columns.values

    samples_keep = []

    for SAD_idx, SAD in enumerate(biom_table_subset_values.transpose()):

        SAD_nonzero_idx = SAD > 0
        SAD_nonzero_taxa = taxa_labels[SAD_nonzero_idx]
        SAD_nonzero = SAD[SAD_nonzero_idx]

        if len(SAD_nonzero) < S_min:
            continue

        #if len(SAD_nonzero) > S_max:
        #    continue

        taxonomy_names.append(SAD_nonzero_taxa)
        SADs.append(SAD_nonzero)

        samples_keep.append(biom_samples[SAD_idx])


    return SADs, taxonomy_names, samples_keep



def get_s_by_s(samples, S_min = 10, rarefied=False):

    if rarefied == True:
        biom_path = '%semp/otu_tables/closed_ref_silva/emp_cr_silva_16S_123.qc_filtered.rare_10000.biom' % config.data_directory

    else:
        biom_path = '%semp/otu_tables/closed_ref_silva/emp_cr_silva_16S_123.qc_filtered.biom' % config.data_directory


    #biom_path = '%semp/otu_tables/closed_ref_silva/emp_cr_silva_16S_123.subset_2k.rare_10000.biom' % config.data_directory

    biom_table = biom.load_table(biom_path)

    taxonomy = biom_table.metadata_to_dataframe('observation')

    SADs = []
    taxonomy_names = []

    # remove few samples that aren't in biom file
    samples_in_biom = [sample for sample in samples if sample in  biom_table.ids() ]

    biom_table_subset = biom_table.filter(samples_in_biom, axis='sample', inplace=False)
    #dense=True
    biom_table_subset_df = biom_table_subset.to_dataframe(dense=True)
    # remove absent species
    biom_table_subset_df = biom_table_subset_df[(biom_table_subset_df.T != 0).any()]

    biom_table_subset_values = biom_table_subset_df.values
    taxa_labels = biom_table_subset_df.index.values
    biom_samples = biom_table_subset_df.columns.values


    return biom_table_subset_values, taxa_labels, biom_samples





def get_s_by_s_deblur_subsample(samples):

    biom_path = '%semp/otu_tables/deblur/emp_deblur_90bp.subset_2k.rare_5000.biom' % config.data_directory

    biom_table = biom.load_table(biom_path)

    taxonomy = biom_table.metadata_to_dataframe('observation')

    SADs = []
    taxonomy_names = []

    # remove few samples that aren't in biom file
    samples_in_biom = [sample for sample in samples if sample in  biom_table.ids() ]

    biom_table_subset = biom_table.filter(samples_in_biom, axis='sample', inplace=False)
    biom_table_subset_df = biom_table_subset.to_dataframe(dense=True)
    # remove absent species
    biom_table_subset_df = biom_table_subset_df[(biom_table_subset_df.T != 0).any()]

    biom_table_subset_values = biom_table_subset_df.values
    taxa_labels = biom_table_subset_df.index.values
    biom_samples = biom_table_subset_df.columns.values


    return biom_table_subset_values, taxa_labels, biom_samples




def get_s_by_s_deblur(samples):

    biom_path = '%semp/otu_tables/deblur/emp_deblur_90bp.qc_filtered.biom' % config.data_directory

    biom_table = biom.load_table(biom_path)

    taxonomy = biom_table.metadata_to_dataframe('observation')

    SADs = []
    taxonomy_names = []

    # remove few samples that aren't in biom file
    samples_in_biom = [sample for sample in samples if sample in  biom_table.ids() ]

    biom_table_subset = biom_table.filter(samples_in_biom, axis='sample', inplace=False)
    biom_table_subset_df = biom_table_subset.to_dataframe(dense=True)
    # remove absent species
    biom_table_subset_df = biom_table_subset_df[(biom_table_subset_df.T != 0).any()]

    biom_table_subset_values = biom_table_subset_df.values
    taxa_labels = biom_table_subset_df.index.values
    biom_samples = biom_table_subset_df.columns.values


    return biom_table_subset_values, taxa_labels, biom_samples








def calculate_shannon_diversity(sad):

    relative_sad = sad / sum(sad)
    relative_sad = relative_sad[relative_sad>0]
    shannon_diversity = -1*sum(relative_sad*numpy.log(relative_sad) )

    #shannon_diversity = -1*sum(relative_sad_i * numpy.log(relative_sad_i) for relative_sad_i in relative_sad)

    return shannon_diversity


def calculate_pielou_evenness(sad):

    return calculate_shannon_diversity(sad) / numpy.log(len(sad))



def calculate_sparsity(sad):

    return sum(sad==0)/len(sad)


def calculate_relative_richness(sad):

    return sum(sad>0)/len(sad)


def calculate_richness(sad):

    return sum(sad>0)


def predict_occupancy(s_by_s, species, totreads=numpy.asarray([])):

    # get squared inverse cv
    # assume that entries are read counts.
    rel_s_by_s_np = (s_by_s/s_by_s.sum(axis=0))

    beta_all = []
    mean_all = []

    for s in rel_s_by_s_np:

        var = numpy.var(s)
        mean = numpy.mean(s)

        beta = (mean**2)/var

        mean_all.append(mean)
        beta_all.append(beta)

    beta_all = numpy.asarray(beta_all)
    mean_all = numpy.asarray(mean_all)


    s_by_s_presence_absence = numpy.where(s_by_s > 0, 1, 0)

    occupancies = s_by_s_presence_absence.sum(axis=1) / s_by_s_presence_absence.shape[1]

    # calcualte total reads if no argument is passed
    # sloppy quick fix
    if len(totreads) == 0:
        totreads = s_by_s.sum(axis=0)

    # calculate mean and variance excluding zeros
    # tf = mean relative abundances
    tf = []
    for afd in s_by_s:
        afd_no_zeros = afd[afd>0]
        tf.append(numpy.mean(afd_no_zeros/ totreads[afd>0]))

    tf = numpy.asarray(tf)
    # go through and calculate the variance for each species

    tvpf_list = []
    for afd in s_by_s:
        afd_no_zeros = afd[afd>0]

        N_reads = s_by_s.sum(axis=0)[numpy.nonzero(afd)[0]]
        tvpf_list.append(numpy.mean(  (afd_no_zeros**2 - afd_no_zeros) / (totreads[afd>0]**2) ))

    tvpf = numpy.asarray(tvpf_list)

    f = occupancies*tf
    vf= occupancies*tvpf

    # there's this command in Jacopo's code %>% mutate(vf = vf - f^2 )%>%
    # It's applied after f and vf are calculated, so I think I can use it
    # This should be equivalent to the mean and variance including zero
    vf = vf - (f**2)

    beta = (f**2)/vf
    theta = f/beta

    predicted_occupancies = []
    # each species has it's own beta and theta, which is used to calculate predicted occupancy
    for beta_i, theta_i in zip(beta,theta):
        predicted_occupancies.append(1 - numpy.mean( ((1+theta_i*totreads)**(-1*beta_i ))   ))

    predicted_occupancies = numpy.asarray(predicted_occupancies)

    species = numpy.asarray(species)
    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
    mad = numpy.mean(rel_s_by_s, axis=1)

    return occupancies, predicted_occupancies, mad, beta, species




def predict_mean_richness(s_by_s, species, totreads=numpy.asarray([])):

    # get squared inverse cv
    # assume that entries are read counts.
    rel_s_by_s_np = (s_by_s/s_by_s.sum(axis=0))

    beta_all = []
    mean_all = []

    for s in rel_s_by_s_np:

        var = numpy.var(s)
        mean = numpy.mean(s)

        beta = (mean**2)/var

        mean_all.append(mean)
        beta_all.append(beta)

    beta_all = numpy.asarray(beta_all)
    mean_all = numpy.asarray(mean_all)


    s_by_s_presence_absence = numpy.where(s_by_s > 0, 1, 0)
    occupancies = s_by_s_presence_absence.sum(axis=1) / s_by_s_presence_absence.shape[1]

    # calcualte total reads if no argument is passed
    # sloppy quick fix
    if len(totreads) == 0:
        totreads = s_by_s.sum(axis=0)

    # calculate mean and variance excluding zeros
    # tf = mean relative abundances
    tf = []
    for afd in s_by_s:
        afd_no_zeros = afd[afd>0]
        tf.append(numpy.mean(afd_no_zeros/ totreads[afd>0]))

    tf = numpy.asarray(tf)
    # go through and calculate the variance for each species

    tvpf_list = []
    for afd in s_by_s:
        afd_no_zeros = afd[afd>0]

        N_reads = s_by_s.sum(axis=0)[numpy.nonzero(afd)[0]]
        tvpf_list.append(numpy.mean(  (afd_no_zeros**2 - afd_no_zeros) / (totreads[afd>0]**2) ))

    tvpf = numpy.asarray(tvpf_list)

    f = occupancies*tf
    vf= occupancies*tvpf

    # there's this command in Jacopo's code %>% mutate(vf = vf - f^2 )%>%
    # It's applied after f and vf are calculated, so I think I can use it
    # This should be equivalent to the mean and variance including zero
    vf = vf - (f**2)

    beta = (f**2)/vf
    theta = f/beta

    richness_observed = s_by_s_presence_absence.sum(axis=0)
    richness_predicted = numpy.asarray([sum(1-((1+theta*totreads_i)**(-1*beta))) for  totreads_i in totreads])

    to_keep = ((~numpy.isnan(richness_observed)) & (~numpy.isnan(richness_predicted)))

    #richness_observed = richness_observed[to_keep]
    #richness_predicted = richness_predicted[to_keep]

    mean_richness_observed = numpy.mean(richness_observed[to_keep])
    mean_richness_predicted = numpy.mean(richness_predicted[to_keep])


    # variance
    #richness_predicted_second_moment = numpy.asarray([sum( (1-((1+theta*totreads_i)**(-1*beta)) ) **2 ) for totreads_i in totreads])
    
    #richness_predicted_variance = numpy.asarray([sum( (1+theta*totreads_i)**(-1*beta) - ((1+theta*totreads_i)**(-1*beta) ) **2 ) for totreads_i in totreads])


    #to_keep = ((~numpy.isnan(richness_observed)) & (~numpy.isnan(richness_predicted_variance)) )

    #richness_predicted_second_moment = richness_predicted_second_moment[to_keep]
    #richness_predicted_variance = richness_predicted_variance[to_keep]
    #richness_observed = richness_observed[to_keep]

    #variance_richness_observed = numpy.var(richness_observed)
    #variance_richness_predicted = numpy.mean(richness_predicted_variance)

    return mean_richness_observed, mean_richness_predicted





def predict_var_richness(s_by_s, species, totreads=numpy.asarray([])):

    rel_s_by_s_np = (s_by_s/s_by_s.sum(axis=0))

    beta_all = []
    mean_all = []

    for s in rel_s_by_s_np:

        var = numpy.var(s)
        mean = numpy.mean(s)

        beta = (mean**2)/var

        mean_all.append(mean)
        beta_all.append(beta)

    beta_all = numpy.asarray(beta_all)
    mean_all = numpy.asarray(mean_all)


    s_by_s_presence_absence = numpy.where(s_by_s > 0, 1, 0)
    occupancies = s_by_s_presence_absence.sum(axis=1) / s_by_s_presence_absence.shape[1]

    # calcualte total reads if no argument is passed
    # sloppy quick fix
    if len(totreads) == 0:
        totreads = s_by_s.sum(axis=0)

    # calculate mean and variance excluding zeros
    # tf = mean relative abundances
    tf = []
    for afd in s_by_s:
        afd_no_zeros = afd[afd>0]
        tf.append(numpy.mean(afd_no_zeros/ totreads[afd>0]))

    tf = numpy.asarray(tf)
    # go through and calculate the variance for each species

    tvpf_list = []
    for afd in s_by_s:
        afd_no_zeros = afd[afd>0]

        N_reads = s_by_s.sum(axis=0)[numpy.nonzero(afd)[0]]
        tvpf_list.append(numpy.mean(  (afd_no_zeros**2 - afd_no_zeros) / (totreads[afd>0]**2) ))

    tvpf = numpy.asarray(tvpf_list)

    f = occupancies*tf
    vf= occupancies*tvpf

    # there's this command in Jacopo's code %>% mutate(vf = vf - f^2 )%>%
    # It's applied after f and vf are calculated, so I think I can use it
    # This should be equivalent to the mean and variance including zero
    vf = vf - (f**2)

    beta = (f**2)/vf
    theta = f/beta

    richness_predicted_variance_term_one = numpy.asarray([sum( (1+theta*totreads_i)**(-1*beta) - ((1+theta*totreads_i)**(-1*beta) ) **2 ) for totreads_i in totreads])

    #prob_absence_over_samples = numpy.mean([])
    #1-((1+theta*totreads_i)**(-1*beta)) for 
    prob_absence_over_samples = []
    for i in range(len(beta)):
        prob_absence_over_samples.append(numpy.mean((1+theta[i]*totreads)**(-1*beta[i])))
    
    prob_absence_over_samples = numpy.asarray(prob_absence_over_samples)

    term_ij_all = 0
    for i in range(len(beta)):

        for j in range(i):

            prob_absence_i = (1+theta[i]*totreads)**(-1*beta[i])
            prob_absence_j = (1+theta[j]*totreads)**(-1*beta[j])

            term_ij = numpy.mean((1-prob_absence_i)*(1-prob_absence_j)) - (1-prob_absence_over_samples[i])*(1-prob_absence_over_samples[j])
            term_ij_all += term_ij
    

    richness_predicted_variance = numpy.mean(richness_predicted_variance_term_one)
    richness_predicted_variance = richness_predicted_variance + 2*term_ij_all

    richness_observed_variance = numpy.var(numpy.sum(s_by_s_presence_absence, axis=0))

    cov_matrix = numpy.cov(s_by_s_presence_absence)
    cov_sum = 2*sum(cov_matrix[numpy.triu_indices(cov_matrix.shape[0], k = 1)])
    richness_predicted_variance_plus_covariance = richness_predicted_variance + cov_sum
    

    return richness_observed_variance, richness_predicted_variance, richness_predicted_variance_plus_covariance





def predict_mean_and_var_richness_and_diversity_using_gamma_rv_with_corr(s_by_s, iter_=100):

    rel_s_by_s_np = (s_by_s/s_by_s.sum(axis=0))

    s_by_s_np_pres_abs = (s_by_s>0)
    richness_observed = s_by_s_np_pres_abs.sum(axis=0)
    mean_richness_observed = numpy.mean(richness_observed)
    var_richness_observed = numpy.var(richness_observed)

    diversity_observed = numpy.apply_along_axis(calculate_shannon_diversity, 0, rel_s_by_s_np)
    mean_diversity_observed = numpy.mean(diversity_observed)
    var_diversity_observed = numpy.var(diversity_observed)

    n_reads = s_by_s.sum(axis=0)

    mean_rel_s_by_s = numpy.mean(rel_s_by_s_np, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s_np, axis=1)

    corr_rel_s_by_s = numpy.corrcoef(rel_s_by_s_np)

    mean_richness_all = []
    mean_richness_corr_all = []

    var_richness_all = []
    var_richness_corr_all = []

    mean_diversity_all = []
    mean_diversity_corr_all = []

    var_diversity_all = []
    var_diversity_corr_all = []

    #import timeit

    for i in range(iter_):

        #start = timeit.default_timer()        
        
        simulation_utils.genrate_community_from_mean_and_var_rv(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads))

        # SVD for identity matrix
        
        rel_abundances_gamma, read_counts_gamma = simulation_utils.genrate_community_from_mean_and_var_svd(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads), svd)
        # genrate_community_from_mean_and_var_svd(mean, var, N, n_sites, svd)

        # All the program statements
        #stop = timeit.default_timer()
        #execution_time = stop - start
        #print(execution_time, "seconds")


        #start = timeit.default_timer()     

        #rel_abundances_gamma_corr, read_counts_gamma_corr = simulation_utils.genrate_community_from_mean_and_var(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads), corr_matrix=corr_rel_s_by_s)

        #stop = timeit.default_timer()
        #execution_time = stop - start
        #print(execution_time, "seconds")



        read_counts_gamma_pres_abs = (read_counts_gamma>0)
        read_counts_gamma_corr_pres_abs = (read_counts_gamma_corr>0)

        richness = read_counts_gamma_pres_abs.sum(axis=0)
        richness_corr = read_counts_gamma_corr_pres_abs.sum(axis=0)

        mean_richness_all.append(numpy.mean(richness))
        mean_richness_corr_all.append(numpy.mean(richness_corr))

        var_richness_all.append(numpy.var(richness))
        var_richness_corr_all.append(numpy.var(richness_corr))

        # diversity
        diversity_null = numpy.apply_along_axis(calculate_shannon_diversity, 0, read_counts_gamma)
        diversity_null_corr = numpy.apply_along_axis(calculate_shannon_diversity, 0, read_counts_gamma_corr)
        
        mean_diversity_all.append(numpy.mean(diversity_null))
        mean_diversity_corr_all.append(numpy.mean(diversity_null_corr))

        var_diversity_all.append(numpy.var(diversity_null))
        var_diversity_corr_all.append(numpy.var(diversity_null_corr))



    return mean_richness_observed, var_richness_observed, \
            numpy.mean(mean_richness_all), numpy.mean(var_richness_all), \
            numpy.mean(mean_richness_corr_all), numpy.mean(var_richness_corr_all), \
            mean_diversity_observed, var_diversity_observed, \
            numpy.mean(mean_diversity_all), numpy.mean(var_diversity_all), \
            numpy.mean(mean_diversity_corr_all), numpy.mean(var_diversity_corr_all)  






def predict_richness_and_diversity_using_gamma_rv(s_by_s, iter_=100):

    rel_s_by_s_np = (s_by_s/s_by_s.sum(axis=0))

    s_by_s_np_pres_abs = (s_by_s>0)
    richness_observed = s_by_s_np_pres_abs.sum(axis=0)
    mean_richness_observed = numpy.mean(richness_observed)
    var_richness_observed = numpy.var(richness_observed)

    diversity_observed = numpy.apply_along_axis(calculate_shannon_diversity, 0, rel_s_by_s_np)
    mean_diversity_observed = numpy.mean(diversity_observed)
    var_diversity_observed = numpy.var(diversity_observed)

    n_reads = s_by_s.sum(axis=0)

    mean_rel_s_by_s = numpy.mean(rel_s_by_s_np, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s_np, axis=1)

    mean_richness_all = []
    var_richness_all = []
    mean_diversity_all = []
    var_diversity_all = []

    # make SVD of identity matrix

    u = numpy.identity(len(mean_rel_s_by_s))
    s = numpy.asarray([1]*len(mean_rel_s_by_s))
    v = numpy.identity(len(mean_rel_s_by_s))
    svd = (u, s, v)

    for i in range(iter_):
        
        rel_abundances_gamma, read_counts_gamma = simulation_utils.genrate_community_from_mean_and_var_svd(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads), svd)

        read_counts_gamma_pres_abs = (read_counts_gamma>0)

        richness = read_counts_gamma_pres_abs.sum(axis=0)

        mean_richness_all.append(numpy.mean(richness))
        var_richness_all.append(numpy.var(richness))

        # diversity
        diversity_null = numpy.apply_along_axis(calculate_shannon_diversity, 0, read_counts_gamma)
       
        mean_diversity_all.append(numpy.mean(diversity_null))
        var_diversity_all.append(numpy.var(diversity_null))



    return mean_richness_all, var_richness_all, mean_diversity_all, var_diversity_all








def predict_mean_and_var_richness_and_diversity_using_gamma_rv(s_by_s, iter_=100):

    rel_s_by_s_np = (s_by_s/s_by_s.sum(axis=0))

    s_by_s_np_pres_abs = (s_by_s>0)
    richness_observed = s_by_s_np_pres_abs.sum(axis=0)
    mean_richness_observed = numpy.mean(richness_observed)
    var_richness_observed = numpy.var(richness_observed)

    diversity_observed = numpy.apply_along_axis(calculate_shannon_diversity, 0, rel_s_by_s_np)
    mean_diversity_observed = numpy.mean(diversity_observed)
    var_diversity_observed = numpy.var(diversity_observed)

    n_reads = s_by_s.sum(axis=0)

    mean_rel_s_by_s = numpy.mean(rel_s_by_s_np, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s_np, axis=1)

    mean_richness_all = []
    var_richness_all = []
    mean_diversity_all = []
    var_diversity_all = []
    var_diversity_poisson_all = []

    # make SVD of identity matrix

    u = numpy.identity(len(mean_rel_s_by_s))
    s = numpy.asarray([1]*len(mean_rel_s_by_s))
    v = numpy.identity(len(mean_rel_s_by_s))
    svd = (u, s, v)

    second_moment_diversity_all = []

    for i in range(iter_):

        rel_abundances_gamma, read_counts_gamma = simulation_utils.genrate_community_from_mean_and_var_rvs(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads))
        
        #rel_abundances_gamma, read_counts_gamma = simulation_utils.genrate_community_from_mean_and_var_svd(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads), svd)

        read_counts_gamma_pres_abs = (read_counts_gamma>0)

        richness = read_counts_gamma_pres_abs.sum(axis=0)

        mean_richness_all.append(numpy.mean(richness))
        var_richness_all.append(numpy.var(richness))

        # diversity
        diversity_null = numpy.apply_along_axis(calculate_shannon_diversity, 0, read_counts_gamma)
       
        mean_diversity_all.append(numpy.mean(diversity_null))
        var_diversity_all.append(numpy.var(diversity_null))

        # repeat for poisson sampling
        #diversity_poisson_null = numpy.apply_along_axis(calculate_shannon_diversity, 0, read_counts_gamma_poisson)
        #var_diversity_poisson_all.append(numpy.var(diversity_poisson_null))





    return mean_richness_observed, var_richness_observed, \
            numpy.mean(mean_richness_all), numpy.mean(var_richness_all), \
            mean_diversity_observed, var_diversity_observed, \
            numpy.mean(mean_diversity_all), numpy.mean(var_diversity_all)  







def prob_n_reads(n, N, mean_, beta_):

    # exp( gammaln(beta+n) - gammaln(n+1) - gammaln(beta) )
    # gamma of factorial results in numerical overflow, do logamma trick instead
    # gamma(beta+n) and gamma(n+1) are large, but their ratio is not, so gammaln(beta+n) - gammaln(n+1) is ok and can be exponentiated

    return numpy.exp( special.gammaln(beta_+n) - special.gammaln(n+1) - special.gammaln(beta_) )   * (((mean_*N)/(beta_ + mean_*N))**n) * ((beta_/(beta_ + mean_*N))**beta_)


def integrand_first_moment(n, N, mean_, beta_):
    return (n/N)*numpy.log(n/N) * prob_n_reads(n, N, mean_, beta_)


def integrand_second_moment(n, N, mean_, beta_):
    return (((n/N)*numpy.log(n/N))**2) * prob_n_reads(n, N, mean_, beta_)



def sum_summation_first_moment(N, mean_, beta_):

    sum_ = 0

    for n in range(1, N+1):
        sum_ += (n/N)*numpy.log(n/N)*prob_n_reads(n, N, mean_, beta_)

    return sum_





def sum_summation_second_moment(N, mean_, beta_):

    sum_ = 0

    for n in range(1, N+1):
        sum_ += (((n/N)*numpy.log(n/N))**2)*prob_n_reads(n, N, mean_, beta_)

    return sum_




def summation_first_moment(n_range, N, mean_, beta_):

    summation = sum((n_range/N)*numpy.log(n_range/N)*numpy.exp( special.gammaln(beta_+n_range) - special.gammaln(n_range+1) - special.gammaln(beta_) )   * (((mean_*N)/(beta_ + mean_*N))**n_range) * ((beta_/(beta_ + mean_*N))**beta_))

    return summation


def summation_second_moment(n_range, N, mean_, beta_):

    summation = sum((((n_range/N)*numpy.log(n_range/N))**2)*numpy.exp( special.gammaln(beta_+n_range) - special.gammaln(n_range+1) - special.gammaln(beta_) )   * (((mean_*N)/(beta_ + mean_*N))**n_range) * ((beta_/(beta_ + mean_*N))**beta_))

    return summation


#def double_integral_second_moment(n, x, N, mean_, beta_):
#    return (((n/N)*numpy.log(n/N))) * special.comb(N, n) * (x**n) * ((1-x)**(N-n)) * (x**(beta_-1)) * numpy.exp(-1*(beta_/mean_)*x) 
    


#N = 10000
#x = 0.0005
#beta=2
#x = summation_first_moment(N, x, beta)

#integrand_second_moment_result = integrate.quad(integrand_first_moment, 0, N, args=(N, x, beta),  epsabs=1e-25)

#print(x, integrand_second_moment_result[0])



def predict_second_moment(s_by_s, iter_=100):

    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))

    n_reads = s_by_s.sum(axis=0)

    mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
    beta_rel_s_by_s = (mean_rel_s_by_s**2)/var_rel_s_by_s

    diversity_second_moment_all = []
    diversity_first_moment_all = []
    # dict with integral for each species for each sample
    product_pairs = []
    for m in range(len(n_reads)):

        N_m = n_reads[m]

        diversity_second_moment_m = []
        diversity_second_moment_second_term_m = []
        diversity_first_moment_m = []
        product_pairs_m = []
        
        for i in range(len(mean_rel_s_by_s)):

            mean_i = mean_rel_s_by_s[i]
            beta_i = beta_rel_s_by_s[i]

            integrand_second_moment_result = integrate.quad(integrand_second_moment, 0, N_m, args=(N_m, mean_i, beta_i),  epsabs=1e-25)
            diversity_second_moment_m.append(numpy.absolute(integrand_second_moment_result[0]))

            integrand_first_moment_result = integrate.quad(integrand_first_moment, 0, N_m, args=(N_m, mean_i, beta_i),  epsabs=1e-25)
            diversity_second_moment_second_term_m.append(integrand_first_moment_result[0])
            diversity_first_moment_m.append(integrand_first_moment_result[0])

        diversity_second_moment_all.append(diversity_second_moment_m)
        diversity_first_moment_all.append(diversity_first_moment_m)

        for i in range(len(mean_rel_s_by_s)):
            for j in range(i):

                product_pairs_m.append(diversity_second_moment_second_term_m[i]*diversity_second_moment_second_term_m[j])

        product_pairs.append(product_pairs_m)


    product_pairs = numpy.asarray(product_pairs)
    expected_product_pairs = numpy.mean(product_pairs, axis=0)

    diversity_second_moment_all = numpy.asarray(diversity_second_moment_all)
    diversity_first_moment_all = numpy.asarray(diversity_first_moment_all)
    expected_diversity_second_moment = numpy.mean(diversity_second_moment_all, axis=0)

    u = numpy.identity(len(mean_rel_s_by_s))
    s = numpy.asarray([1]*len(mean_rel_s_by_s))
    v = numpy.identity(len(mean_rel_s_by_s))
    svd = (u, s, v)

    expected_diversity_second_moment_rvs_all = []
    product_pairs_null = []
    var_diversity_rvs_all = []
    for i in range(iter_):
        
        rel_abundances_gamma, read_counts_gamma = simulation_utils.genrate_community_from_mean_and_var_svd(mean_rel_s_by_s, var_rel_s_by_s, n_reads, len(n_reads), svd)
        rel_read_counts_gamma = read_counts_gamma/numpy.sum(read_counts_gamma, axis=0)

        measure_rel_read_counts_gamma = (rel_read_counts_gamma*numpy.log(rel_read_counts_gamma))
        measure_rel_read_counts_gamma[numpy.isnan(measure_rel_read_counts_gamma)] = 0

        #rel_read_counts_gamma_first_moment = measure_rel_read_counts_gamma*1
        #rel_read_counts_gamma_first_moment[numpy.isnan(rel_read_counts_gamma_first_moment)] = 0

        # diversity
        var_diversity_rvs = numpy.var(numpy.apply_along_axis(calculate_shannon_diversity, 0, rel_read_counts_gamma))
        var_diversity_rvs_all.append(var_diversity_rvs)

        expected_diversity_first_moment_rvs = numpy.mean(measure_rel_read_counts_gamma, axis=1)
        
        rel_read_counts_gamma_second_moment = measure_rel_read_counts_gamma**2
        rel_read_counts_gamma_second_moment[numpy.isnan(rel_read_counts_gamma_second_moment)] = 0

        expected_diversity_second_moment_rvs = numpy.mean(rel_read_counts_gamma_second_moment, axis=1)
        expected_diversity_second_moment_rvs_all.append(expected_diversity_second_moment_rvs)

        product_pairs_null_i = []
        # product pairs
        for i in range(len(expected_diversity_first_moment_rvs)):
            for j in range(i):

                product_pairs_null_i.append(expected_diversity_first_moment_rvs[i] * expected_diversity_first_moment_rvs[j])

        product_pairs_null.append(product_pairs_null_i)

    
    product_pairs_null = numpy.asarray(product_pairs_null)
    expected_product_pairs_null_null = numpy.mean(product_pairs_null, axis=0)

    expected_diversity_second_moment_rvs_all = numpy.asarray(expected_diversity_second_moment_rvs_all)
    expected_diversity_second_moment_rvs = numpy.mean(expected_diversity_second_moment_rvs_all, axis=0)

    
    print(len(numpy.sum(expected_diversity_second_moment_rvs_all, axis=1)))
    print('first term')
    print(sum(expected_diversity_second_moment), numpy.mean(numpy.sum(expected_diversity_second_moment_rvs_all, axis=1)))

    # print sum
    print('second term')
    print(sum(expected_product_pairs), numpy.mean(numpy.sum(product_pairs_null, axis=1)))


    #print('second moment')
    # variance rvs
    #variance_rvs = numpy.sum(expected_diversity_second_moment_rvs_all, axis=1) + numpy.sum(product_pairs_null, axis=1)
    #print(sum(expected_diversity_second_moment) + sum(expected_product_pairs), numpy.mean(variance_rvs), numpy.mean(var_diversity_rvs_all))

    print('variance')

    #sum(expected_diversity_second_moment) + sum(expected_product_pairs)

    var_predicted = numpy.mean(numpy.sum(diversity_second_moment_all, axis=1) + numpy.sum(product_pairs, axis=1) - (numpy.sum(diversity_first_moment_all, axis=1)**2))

    #(numpy.sum(diversity_first_moment_all, axis=1)**2)

    print(var_predicted, numpy.mean(var_diversity_rvs_all))

    
    
    return expected_diversity_second_moment, expected_diversity_second_moment_rvs, expected_product_pairs, expected_product_pairs_null_null







def predict_mean_and_var_diversity_analytic(s_by_s, species):

    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
    diversity_observed = numpy.apply_along_axis(calculate_shannon_diversity, 0, rel_s_by_s)
    mean_diversity_observed = numpy.mean(diversity_observed)
    var_diversity_observed = numpy.var(diversity_observed)

    x_log_x = rel_s_by_s*numpy.log(rel_s_by_s)
    x_log_x[numpy.isnan(x_log_x)] = 0
    cov_x_log_x = numpy.cov(x_log_x)
    cov_sum = 2*sum(cov_x_log_x[numpy.triu_indices(cov_x_log_x.shape[0], k = 1)])
    
    n_reads = s_by_s.sum(axis=0)

    mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
    beta_rel_s_by_s = (mean_rel_s_by_s**2)/var_rel_s_by_s

    #diversity_first_moment_all = []
    #diversity_second_moment_all = []

    mean_all = []
    var_all = []

    # dict with integral for each species for each sample
    for m in range(len(n_reads)):

        N_m = int(n_reads[m])

        diversity_first_moment_m = 0
        diversity_second_moment_m = 0

        diversity_first_moment_sum_m = 0
        #diversity_second_moment_sum_m = 0

        integrand_first_moment_all = []
        integrand_second_moment_all = []
        for i in range(len(mean_rel_s_by_s)):
            
            mean_i = mean_rel_s_by_s[i]
            beta_i = beta_rel_s_by_s[i]

            integral_first_moment_result = integrate.quad(integrand_first_moment, 0, N_m, args=(N_m, mean_i, beta_i), epsabs=1e-20)
            #diversity_first_moment_m += integral_first_moment_result[0]
            integrand_first_moment_all.append(integral_first_moment_result[0])

            integrand_second_moment_result = integrate.quad(integrand_second_moment, 0, N_m, args=(N_m, mean_i, beta_i),  epsabs=1e-25)
            integrand_second_moment_all.append(integrand_second_moment_result[0])
            #diversity_second_moment_m += numpy.absolute(integrand_second_moment_result[0])

            # summation
            #summation_first_moment_result = sum_summation_first_moment(N_m, mean_i, beta_i)
            #diversity_first_moment_sum_m += summation_first_moment_result
            #integrand_first_moment_all.append(summation_first_moment_result)

            #summation_second_moment_result = sum_summation_second_moment(N_m, mean_i, beta_i)
            #diversity_second_moment_m += numpy.absolute(summation_second_moment_result)
        
        integrand_first_moment_all = numpy.absolute(integrand_first_moment_all)
        integrand_second_moment_all = numpy.absolute(integrand_second_moment_all)

        diversity_second_moment_second_term_m = 0
        for i in range(len(mean_rel_s_by_s)):
            for j in range(i):
                diversity_second_moment_second_term_m += integrand_first_moment_all[i]*integrand_first_moment_all[j]

        mean_m = sum(integrand_first_moment_all)
        var_m = sum(integrand_second_moment_all) + 2*diversity_second_moment_second_term_m - (mean_m**2)

        mean_all.append(mean_m)
        var_all.append(var_m)

        #diversity_second_moment_m = diversity_second_moment_m+(2*diversity_second_moment_second_term_m)
        #diversity_first_moment_all.append(diversity_first_moment_m)
        #diversity_second_moment_all.append(diversity_second_moment_m+(2*diversity_second_moment_second_term_m))

        #diversity_first_moment_sum_all.append(diversity_first_moment_sum_m)
        #diversity_second_moment_sum_all.append(diversity_second_moment_sum_m+(2*diversity_second_moment_second_term_sum_m))



    mean_diversity_predicted = numpy.mean(mean_all) 
    var_diversity_predicted = numpy.mean(var_all)
    #var_diversity_predicted = numpy.mean(diversity_second_moment_all - (mean_diversity_predicted**2))

    #mean_diversity_predicted_sum = numpy.mean(numpy.absolute(diversity_first_moment_sum_all)) 
    #var_diversity_predicted_sum = numpy.mean(diversity_second_moment_sum_all) - (mean_diversity_predicted_sum**2)

    var_diversity_predicted_plus_covariance = var_diversity_predicted+cov_sum

    return mean_diversity_observed, mean_diversity_predicted, var_diversity_observed, var_diversity_predicted, var_diversity_predicted_plus_covariance




def predict_mean_and_var_diversity_analytic_with_integral(s_by_s, species, first_moment_integral, second_moment_integral):

    rel_s_by_s = (s_by_s/s_by_s.sum(axis=0))
    diversity_observed = numpy.apply_along_axis(calculate_shannon_diversity, 0, rel_s_by_s)
    mean_diversity_observed = numpy.mean(diversity_observed)
    var_diversity_observed = numpy.var(diversity_observed)

    x_log_x = rel_s_by_s*numpy.log(rel_s_by_s)
    x_log_x[numpy.isnan(x_log_x)] = 0
    cov_x_log_x = numpy.cov(x_log_x)
    cov_sum = 2*sum(cov_x_log_x[numpy.triu_indices(cov_x_log_x.shape[0], k = 1)])
    
    n_reads = s_by_s.sum(axis=0)

    mean_rel_s_by_s = numpy.mean(rel_s_by_s, axis=1)
    var_rel_s_by_s = numpy.var(rel_s_by_s, axis=1)
    beta_rel_s_by_s = (mean_rel_s_by_s**2)/var_rel_s_by_s

    diversity_first_moment_all = []
    diversity_second_moment_all = []
    # dict with integral for each species for each sample
    for m in range(len(n_reads)):

        integrand_first_moment_all = first_moment_integral[:,m]
        integrand_second_moment_all = second_moment_integral[:,m]

        diversity_first_moment_m = numpy.absolute(numpy.sum(integrand_first_moment_all))
        diversity_second_moment_m = numpy.absolute(numpy.sum(integrand_second_moment_all))

        diversity_second_moment_second_term_m = 0
        for i in range(len(mean_rel_s_by_s)):
            for j in range(i):
                diversity_second_moment_second_term_m += numpy.absolute(integrand_first_moment_all[i]*integrand_first_moment_all[j])

        diversity_second_moment_m = diversity_second_moment_m+(2*diversity_second_moment_second_term_m)

        diversity_first_moment_all.append(diversity_first_moment_m)
        diversity_second_moment_all.append(diversity_second_moment_m)


    mean_diversity_predicted = numpy.mean(numpy.absolute(diversity_first_moment_all)) 
    var_diversity_predicted = numpy.mean(diversity_second_moment_all) - (mean_diversity_predicted**2)

    print(numpy.mean(diversity_second_moment_all))

    var_diversity_predicted_plus_covariance = var_diversity_predicted+cov_sum

    return mean_diversity_observed, mean_diversity_predicted, var_diversity_observed, var_diversity_predicted, var_diversity_predicted_plus_covariance











def calculate_error(observed, predicted):

    # remove error of zero since we look at the log

    #idx_to_keep = (observed>0) & (predicted > 0)
    #observed = observed[idx_to_keep]
    #predicted = predicted[idx_to_keep]

    error = numpy.absolute(observed - predicted)/observed
    error = error[(~(numpy.isnan(error))) &  (~(numpy.isneginf(error)))  &  (~(numpy.isinf(error))) & (error>0)]

    return error







def get_hist_and_bins(flat_array, bins=20):

    # make sure its an array
    flat_array = numpy.asarray(flat_array)

    flat_array = flat_array[~numpy.isnan(flat_array)]

    # null is too large, so we are binning it for the plot in this script
    hist_, bin_edges_ = numpy.histogram(flat_array, density=True, bins=bins)
    bins_mean_ = numpy.asarray([0.5 * (bin_edges_[i] + bin_edges_[i+1]) for i in range(0, len(bin_edges_)-1 )])
    hist_to_plot = hist_[hist_>0]
    bins_mean_to_plot = bins_mean_[hist_>0]

    return hist_to_plot, bins_mean_to_plot



#def lognormal_probability(cumK, mu, s):



def Klogn(cumK, c, mu0=-19, s0=5):
    # This function estimates the parameters (mu, s) of the lognormal distribution of K
    m1 = numpy.mean(numpy.log(cumK[cumK>c]))
    m2 = numpy.mean(numpy.log(cumK[cumK>c])**2)
    xmu = sympy.symbols('xmu')
    xs = sympy.symbols('xs')
    eq1 = -m1+xmu + sympy.sqrt(2/sympy.pi)*xs*sympy.exp(-((sympy.log(c)-xmu)**2)/2/(xs**2))/(sympy.erfc((sympy.log(c)-xmu)/sympy.sqrt(2)/xs))
    eq2 = -m2+xs**2+m1*xmu+sympy.log(c)*m1-xmu*sympy.log(c)

    sol = sympy.nsolve([eq1,eq2],[xmu,xs],[mu0,s0])

    return float(sol[0]), float(sol[1])


def get_best_Klogn(cumK):

    #min_c = 0.0000001
    # pad max_c because high c values fail to converge
    #max_c = (1/min(s_by_s.sum(axis=0)))*100

    count = 0
    #signal.signal(signal.SIGALRM, timeout_handler)
    #while (count_less < iter) or (count_greater < iter):
    while (count < iter):

        c_i = numpy.random.uniform(numpy.log10(min_c), numpy.log10(max_c))
        c_i = 10**c_i
        sigma_i = numpy.random.uniform(0.2, 1.99)

        try:
            n = RunWithTimeout(Klogn, (cumK, c_i))
            mu, s = n.run(10)

        except:
            continue

        S_tot = int(2*S_obs / erf((numpy.log(c_i)-mu) / numpy.sqrt(2)*s ))

        # calculate likelihood


    # return parameter estimates for c that maximes the likelihood



def sort_s_by_s_by_rank_labels(sad_annotated_dict, rank):

    taxa = list(sad_annotated_dict['taxa'].keys())
    s_by_s = numpy.asarray([sad_annotated_dict['taxa'][t]['abundance'] for t in taxa])
    all_coarse = list([sad_annotated_dict['taxa'][t][rank] for t in taxa])

    # sort taxa by coarse labels
    zip_sorted_coarse_taxa = sorted(zip(all_coarse, taxa), key=lambda x: x[0])
    zip_sorted_coarse = numpy.asarray([x[0] for x in zip_sorted_coarse_taxa])
    zip_sorted_taxa = numpy.asarray([x[1] for x in zip_sorted_coarse_taxa])

    idx_sort = numpy.asarray([taxa.index(x) for x in zip_sorted_taxa])
    s_by_s_sorted = s_by_s[idx_sort,:]

    # now get counts for each item
    unique_coarse, counts_coarse = numpy.unique(zip_sorted_coarse, return_counts=True)
    coarse_grain_by_genus_idx = numpy.append([0], numpy.cumsum(counts_coarse))[:-1]

    return s_by_s_sorted, idx_sort, coarse_grain_by_genus_idx




