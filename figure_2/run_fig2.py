import lomap

"""
Running this script with HiMap installed will reproduce the results in Fig. 2.
Must vary epsilon in the interactive run to output D and E.
"""
# *****************************************************************************
# This example file was written by Dr. Mary Pitman. 2022
# *****************************************************************************

#-------------------------------------------------------#
# Define input files, read data.
#-------------------------------------------------------#
# Input files for weight scores and ligand names.
sim_scores_in = 'similarities.csv'
IDs_in = 'mol_names.txt'

# Read files
n_arr, ID_list = lomap.read_data(sim_scores_in, IDs = IDs_in)

#-------------------------------------------------------#
# Clustering.
#-------------------------------------------------------#
# Perform clustering.
#   sub_arr, sub_ID:   the n_arr and ID_list subdivided by clusters
#   selected_clusters: user selected clusters during interaction.
sub_arr, sub_ID, selected_clusters = lomap.cluster_interactive(n_arr, ID_list)
