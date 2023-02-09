import lomap
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
"""
This script performs analysis of cluster similarities.
"""
# *****************************************************************************
# This example file was written by Dr. Mary Pitman. 2023
# *****************************************************************************

#-------------------------------------------------------#
# Define input files, read data.
#-------------------------------------------------------#
# Input files for weight scores and ligand names.
sim_scores_in = '../test/optimize/sim_scores.csv'
IDs_in = '../test/optimize/mol_names.txt'

# Read files, clean any potential NaN scores.
#   Added optional parameter:
#             delimiter: default is ','
n_arr, ID_list = lomap.read_data(sim_scores_in, IDs = IDs_in)


# Get indices of labels for pandas df
def output_nums(ID_list, cluster):
    ID_numbers = []
    for j in cluster:
        num = ID_list.index(j)
        ID_numbers.append(num)
    return ID_numbers

'''
#-------------------------------------------------------#
# Clustering. Uncomment to output "cluster_IDs.json"
#-------------------------------------------------------#
# Perform clustering.
   sub_arr, sub_ID:   the n_arr and ID_list subdivided by clusters
   selected_clusters: user selected clusters during interaction.
sub_arr, sub_ID, selected_clusters = lomap.cluster_interactive(n_arr, ID_list)
'''

file = open("cluster_IDs.json", 'r')
json_data = json.load(file)

df = pd.DataFrame(data = n_arr,
                  index = range(len(ID_list)),
                  columns = ID_list)

# This LOMAP fix is typically applied prior to clustering
#np.fill_diagonal(df.values, 1.0)

# get the column names for clusters.
cluster0 = json_data['0']
cluster1 = json_data['1']
cluster2 = json_data['2']

# get the index numbers for clusters.
cluster0_nums = output_nums(ID_list, cluster0)
cluster1_nums = output_nums(ID_list, cluster1)
cluster2_nums = output_nums(ID_list, cluster2)

# Subset cluster 0
df_col = df[cluster0]
df_0 = df_col.iloc[cluster0_nums]

# Subset cluster 1
df_col = df[cluster1]
df_1 = df_col.iloc[cluster1_nums]

# Subset cluster 2
df_col = df[cluster2]
df_2 = df_col.iloc[cluster2_nums]

# Combine subsets
cluster01 = cluster0 + cluster1
cluster01_nums = cluster0_nums +  cluster1_nums

cluster12 = cluster1 + cluster2
cluster12_nums = cluster1_nums +  cluster2_nums

cluster20 = cluster2 + cluster0
cluster20_nums = cluster2_nums +  cluster0_nums

# Subset cluster 0, 1
df_col = df[cluster01]
df_01 = df_col.iloc[cluster01_nums]
# Subset cluster 1, 2
df_col = df[cluster12]
df_12 = df_col.iloc[cluster12_nums]
# Subset cluster 2, 0
df_col = df[cluster20]
df_20 = df_col.iloc[cluster20_nums]

# Print the shapes. These are the dfs to work with in joyplots
print(df_0.shape)
print(df_1.shape)
print(df_2.shape)
print(df_01.shape)
print(df_12.shape)
print(df_20.shape)


# (sum(larger block A+B) - (sum(A)+sum(B)))/(sum(A)+sum(B))

def get_sum(df):
   sum = (df.to_numpy().sum())#/(len(df.index)*len(df.index))
   print(sum)
   return sum

# Sums of clusters
sum0 = get_sum(df_0)
sum1 = get_sum(df_1)
sum2 = get_sum(df_2)

# off diagonal Sums of groupings
sum01 = get_sum(df_01) - sum0 - sum1
sum12 = get_sum(df_12) - sum1 - sum2
sum20 = get_sum(df_20) - sum2 - sum0

print(sum01)
print(sum12)
print(sum20)

# NUmber of entries
n0 = len(df_0.index)*len(df_0.index)
n1 = len(df_1.index)*len(df_1.index)
n2 = len(df_2.index)*len(df_2.index)

# off diagonal
n01 = len(df_01.index)*len(df_01.index) - (n0+n1)
n12 = len(df_12.index)*len(df_12.index) - (n1+n2)
n20 = len(df_20.index)*len(df_20.index) - (n0+n2)
'''
#Average entries
av0 = sum0 / n0
av1 = sum1 / n1
av2 = sum2 / n2

av01 = sum01 / n01
av12 = sum12 / n12
av20 = sum20 / n20

print(f'The average value for cluster 0, 1, 2 are {av0}, {av1}, and {av2}')
print(f'The average value for off diagonal values between cluster 0 to 1, 1 to 2, 2 to 0 are {av01}, {av12}, and {av20}')

d = {'0': [av0, av01, av20], '1': [av01, av1, av12], '2': [av20, av12, av2]}
df_heatmap = pd.DataFrame(data=d)

plt.imshow(df_heatmap, cmap ="inferno")

plt.colorbar()
plt.show()
'''


import os
import subprocess
import glob

from joypy import joyplot
from matplotlib import cm

#--------------------------------------------#
# Hard coded variables.
#--------------------------------------------#

mol_names = ['0', '1', '2', '0, 1', '1, 2', '2, 0']

print(df_0.shape)
print(df_1.shape)
print(df_2.shape)
print(df_01.shape)
print(df_12.shape)
print(df_20.shape)



#--------------------------------------------#
#  Driver code for joyplots.
#--------------------------------------------#
#I don't have the proper off diagonals
df_col = df[cluster0]
df_01 = df_col.iloc[cluster1_nums]

df_col = df[cluster1]
df_12 = df_col.iloc[cluster2_nums]

df_col = df[cluster2]
df_20 = df_col.iloc[cluster0_nums]


df0 = pd.DataFrame(df_0.to_numpy().flatten())
df1 = pd.DataFrame(df_1.to_numpy().flatten())
df2 = pd.DataFrame(df_2.to_numpy().flatten())
df01 = pd.DataFrame(df_01.to_numpy().flatten())
df12= pd.DataFrame(df_12.to_numpy().flatten())
df20 = pd.DataFrame(df_20.to_numpy().flatten())


#-----------------------------
# Uncomment for maximum scores
#-----------------------------

df0_max = df_0.max(axis=1)
print(df0_max)
df1_max = df_1.max(axis=1)
df2_max = df_2.max(axis=1)
df01_max = df_01.max(axis=1)
df12_max = df_12.max(axis=1)
df20_max = df_20.max(axis=1)

df0 = pd.DataFrame(df0_max)
df1 = pd.DataFrame(df1_max)
df2 = pd.DataFrame(df2_max)
df01 = pd.DataFrame(df01_max)
df12= pd.DataFrame(df12_max)
df20 = pd.DataFrame(df20_max)

#-----------------------------
# Max heatmap
#-----------------------------

#The maximums of the whole groups
av0 = df0_max.max()
print(av0)
av1 = df1_max.max()
av2 = df2_max.max()

av01 = df01_max.max()
av12 = df12_max.max()
av20 = df20_max.max()

print(f'The max value for cluster 0, 1, 2 are {av0}, {av1}, and {av2}')
print(f'The max value for off diagonal values between cluster 0 to 1, 1 to 2, 2 to 0 are {av01}, {av12}, and {av20}')

d = {'0': [av0, av01, av20], '1': [av01, av1, av12], '2': [av20, av12, av2]}
df_heatmap = pd.DataFrame(data=d)

plt.imshow(df_heatmap, cmap ="inferno")

plt.colorbar()
plt.show()


'''
#-----------------------------

df = pd.concat([df0, df1, df2, df01, df12, df20], ignore_index=True, axis=1, names = mol_names)

df.columns = mol_names
print(df.tail)
plt.figure()
joyplot(df, figsize=(10,12), alpha=0.6, x_range = [0, 1], overlap = 0.8, kind="normalized_counts", bins=25, colormap=cm.inferno, linewidth=1)
#plt.xlim(0.8, 1)


plt.tight_layout()
#plt.title('Maximum Similarity Score Per Ligand', fontsize=18)
plt.title('Similarity Scores', fontsize=18)
plt.show()
#plt.savefig('confidence_values_joyplot.png')
'''
