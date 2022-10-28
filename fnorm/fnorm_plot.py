# Code author: Mary Pitman 2022; mpitman@uci.edu

import numpy as np
import warnings
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import pandas as pd
from scipy.ndimage.filters import gaussian_filter1d

###################################
# File import, parse.
###################################
plt.rcParams['font.size'] = '16'

# Load File
fobj = np.loadtxt('fnorm_data_cv_n_10.txt', comments='#')
fobj2 = np.loadtxt('fnorm_data_cv_n_20.txt', comments='#')
fobj3 = np.loadtxt('fnorm_data_cv_n_40.txt', comments='#')

#Parse files
#X vals
k_10 = fobj[:,0]
k_20 = fobj2[:,0]
k_40 = fobj3[:,0]
# Y vals
d_10 = fobj[:,1]
d_20 = fobj2[:,1]
d_40 = fobj3[:,1]


# Group data by the k index
# Find indices where k changes values
n_disect = np.where(k_10[:-1] != k_10[1:])[0]
n_bounds = np.append(n_disect+1, len(k_10)+1)
n_groups = len(n_bounds)

n_disect2 = np.where(k_20[:-1] != k_20[1:])[0]
n_bounds2 = np.append(n_disect2+1, len(k_20)+1)
n_groups2 = len(n_bounds2)

n_disect3 = np.where(k_40[:-1] != k_40[1:])[0]
n_bounds3 = np.append(n_disect3+1, len(k_40)+1)
n_groups3 = len(n_bounds3)


# Generate arrays of averages and stdevs
# format: [n_edges averages stdev]
avg_arr = np.zeros(0)
stdev_arr = np.zeros(0)
k_list = np.zeros(0)
init = 0
for i in n_bounds:
    a = d_10[init:i]
    k = k_10[init]
    entry = np.average(a)
    stdev = np.std(a)
    init = init + 20
    avg_arr = np.append(avg_arr, entry)
    stdev_arr = np.append(stdev_arr, stdev)
    k_list = np.append(k_list, k)
     
avg_arr2 = np.zeros(0)
stdev_arr2 = np.zeros(0)
k_list2 = np.zeros(0)
init = 0
for i in n_bounds2:
    a = d_20[init:i]
    k = k_20[init]
    entry = np.average(a)
    stdev = np.std(a)
    init = init + 20
    avg_arr2 = np.append(avg_arr2, entry)
    stdev_arr2 = np.append(stdev_arr2, stdev)
    k_list2 = np.append(k_list2, k)
    
avg_arr3 = np.zeros(0)
stdev_arr3 = np.zeros(0)
k_list3 = np.zeros(0)
init = 0
for i in n_bounds3:
    a = d_40[init:i]
    k = k_40[init]
    entry = np.average(a)
    stdev = np.std(a)
    init = init + 20
    avg_arr3 = np.append(avg_arr3, entry)
    stdev_arr3 = np.append(stdev_arr3, stdev)
    k_list3 = np.append(k_list3, k)
    

###################################
# Functions.
###################################
cmap = cm.get_cmap('inferno')
max_colors = n_groups + 24
color_number = 0

def restart_colors():
    global color_number
    color_number = 0
    #np.random.seed(1)

def next_color():
    global color_number
    color_number += 1
    #color = tuple(np.random.uniform(0.0, 0.5, 3))
    color = cmap( ((5 * color_number) % max_colors) / max_colors )
    return color

###################################
# Plotting of curves
###################################
# Calculate averages and shaded error regions.
avg_arr += np.random.normal(0, stdev_arr / 10, avg_arr.size)
upper = gaussian_filter1d(avg_arr + stdev_arr, sigma=1.2)
lower = gaussian_filter1d(avg_arr - stdev_arr, sigma=1.2)

avg_arr2 += np.random.normal(0, stdev_arr2 / 5, avg_arr2.size)
upper2 = gaussian_filter1d(avg_arr2 + stdev_arr2, sigma=1.2)
lower2 = gaussian_filter1d(avg_arr2 - stdev_arr2, sigma=1.2)

avg_arr3 += np.random.normal(0, stdev_arr3 / 5, avg_arr3.size)
upper3 = gaussian_filter1d(avg_arr3 + stdev_arr3, sigma=1.2)
lower3 = gaussian_filter1d(avg_arr3 - stdev_arr3, sigma=1.2)

# Three subplots.
fig, ax = plt.subplots(ncols=3, sharey=False, figsize=(11, 4))

# Plot curves and errors.
ax[0].plot(k_list, avg_arr, marker='.', linestyle = 'dashed', color='#cf7202')
ax[0].fill_between(k_list, upper, lower, color='#fca612', alpha=0.5)

ax[1].plot(k_list2, avg_arr2, marker='.', linestyle = 'dashed', color='#cf7202')
ax[1].fill_between(k_list2, upper2, lower2, color='#fca612', alpha=0.5)

ax[2].plot(k_list3, avg_arr3, marker='.', linestyle = 'dashed', color='#cf7202')
ax[2].fill_between(k_list3, upper3, lower3, color='#fca612', alpha=0.5)

ax[0].plot(k_list[4], avg_arr[4], marker='o', color='black')
ax[1].plot(k_list2[5], avg_arr2[5], marker='o', color='black')
ax[2].plot(k_list3[6], avg_arr3[6], marker='o', color='black', linestyle='none', label='k = $n\cdot\ln(n)$')
ax[2].legend(loc="upper right")

# Format plots.
ax[0].set_title('n = 10')
ax[1].set_title('n = 20')
ax[2].set_title('n = 40')

# Labels
ax[0].set_xlabel("edges, k")
ax[1].set_xlabel("edges, k")
ax[2].set_xlabel("edges, k")
ax[0].set_ylabel("Distance to $D_{full}$")

# Limits.
ax[0].set_ylim(-5, (max(upper)+10))
ax[1].set_ylim(-5, (max(upper2)+10))
ax[2].set_ylim(-15, (max(upper3)+20))

# Display, save.
#plt.savefig('cycles_surface.png', dpi=500)
plt.tight_layout()
plt.show()



