# panel_config.R
# Auth: james.opzoomer@kcl.ac.uk
# Desc: Config file to set up the FCS panel parameters for clustering
# Date: 29/01/2018

# To do:
# Full autoscript test
# Remove the dependancy on exact columns

# projectID
pid = '01'

# This is the ddirectory where the data is stored
gswd = '~/Desktop/Arnold/Data/2018_01_Louvain_clustering/second_validation/'

# This is the filepath to your FCS csv file
fcs_file = "scale_F480_TAM.csv"

## This is the panel paramter setup

# These are a list of proteins stained for i.e CD34 . . 
# These must MATCH THE ORDER OF THE FLUOROPHORES LISTED BELOW
markers = c("FSC_A", "SSC_A", "CD206", "CD11c", "CD11b", "MHCII", "F480")
n_markers = length(markers)

# These are the list of Fluorphores used
# This must MATCH THE ORDER OF THE MARKERS LISTED ABOVE
# They must MATCH THE COLUMN ORDER OF THE FLUOROPHRES IN THE FCS FILE/XLS
stains = c("FSC_A", "SSC_A", "APC", "APC-Cy7", "Amcyan", "PE", "Pacific Blue")

dfPanel = data.frame('markers' = markers, 'stains' = stains) 

# These are the markers that you want to run into the clustering algorithm
# They can be in any order but must match the name of the markers
tSNE_markers = c()

# The number of cytometry events to take for tSNE if the file is too big
subsample = 10000
