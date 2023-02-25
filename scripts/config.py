import os

#directory = os.path.expanduser("~/GitHub/macroeco_phylo/")

#data_directory = '%sdata/' % directory


if os.geteuid() == 501:
    data_directory = os.path.expanduser("~/GitHub/macroeco_phylo/data/")
    analysis_directory = os.path.expanduser("~/GitHub/macroeco_phylo/analyses/")

elif os.geteuid() == 1378857:
    data_directory = os.path.expanduser("/N/u/wrshoema/Carbonate/macroeco_phylo/data/")
    analysis_directory = os.path.expanduser("/N/u/wrshoema/Carbonate/macroeco_phylo/analyses/")

else:
    data_directory = os.path.expanduser("/u/project/ngarud/wrshoema/macroeco_phylo/data/")
    analysis_directory = os.path.expanduser("/u/project/ngarud/wrshoema/macroeco_phylo/analyses/")
