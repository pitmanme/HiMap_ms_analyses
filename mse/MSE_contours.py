# Code author: Mary Pitman 2022; mpitman@uci.edu

import numpy as np
import warnings
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
import pandas as pd

###################################
# File import, parse.
###################################
# Load File
fobj = np.loadtxt('MSE_sim_data.txt', comments='#')

dat = ["n", "k", "A", "D", "rand"]
columns = np.arange(len(dat))

# Parse columns of text data file
for i,j in zip(dat, columns):
    globals()[i] = fobj[:,j]

# Group data by the n index
# Find row indices where n changes values
n_disect = np.where(n[:-1] != n[1:])[0]
n_bounds = np.append(n_disect+1, len(n)+1)
n_groups = len(n_bounds)

# Parse rows of text data file
# Parse edge number data
k1 = k[0:n_bounds[0]]
k2 = k[n_bounds[0]:n_bounds[1]]
k3 = k[n_bounds[1]:n_bounds[2]]
k4 = k[n_bounds[2]:n_bounds[3]]

# Parse node number data.
n1 = n[0:n_bounds[0]]
n2 = n[n_bounds[0]:n_bounds[1]]
n3 = n[n_bounds[1]:n_bounds[2]]
n4 = n[n_bounds[2]:n_bounds[3]]

# A-optimal graphs subplot
A1 = A[0:n_bounds[0]]
A2 = A[n_bounds[0]:n_bounds[1]]
A3 = A[n_bounds[1]:n_bounds[2]]
A4 = A[n_bounds[2]:n_bounds[3]]

# D-optimal graphs subplot
D1 = D[0:n_bounds[0]]
D2 = D[n_bounds[0]:n_bounds[1]]
D3 = D[n_bounds[1]:n_bounds[2]]
D4 = D[n_bounds[2]:n_bounds[3]]

# Random graphs subplot
rand1 = rand[0:n_bounds[0]]
rand2 = rand[n_bounds[0]:n_bounds[1]]
rand3 = rand[n_bounds[1]:n_bounds[2]]
rand4 = rand[n_bounds[2]:n_bounds[3]]


# Dividing points in the data file
A_splits = [A1, A2, A3, A4]
D_splits = [D1, D2, D3, D4]
rand_splits = [rand1, rand2, rand3, rand4]
k_splits = [k1, k2, k3, k4]
n_splits = [n4, n3, n2, n1]

###################################
# Functions.
###################################
cmap = cm.get_cmap('inferno')
max_colors = n_groups + 24
color_number = 0

def restart_colors():
    global color_number
    color_number = 0

def next_color():
    global color_number
    color_number += 1
    color = cmap( ((5 * color_number) % max_colors) / max_colors )
    return color

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx

def curve_plot(ax, Y, label):
    # Reset colors between plots.
    restart_colors()
    ind = 0
    d = []
    # Plot each isonode.
    for i, j in zip(reversed(k_splits), reversed(Y)):
        color = next_color()
        # Scatter plot of data
        ax.scatter(i, j, color=color, marker='.', alpha=0.8)
        
        # Interpolate to smooth plot and make curve
        X_Y_Spline = make_interp_spline(i, j)
        X_ = np.linspace(i.min(), i.max(), 2000)
        Y_ = X_Y_Spline(X_)
        
        # Calculate value of n of isonode
        n_value = n_splits[ind]
        n_i = n_value[0]
        # Define edge counts of interet
        nlnn = n_i*np.log(n_i)
        nx2 = n_i*2.0
        nx1 = n_i
        nx1p5 = n_i*1.5
        nx3 = n_i*3

        # Find value in spline closest to edge count.
        xln, idx = find_nearest(X_, nlnn)
        x_2, idx2 = find_nearest(X_, nx2)
        x_3, idx3 = find_nearest(X_, nx1)
        x_4, idx4 = find_nearest(X_, nx1p5)
        x_5, idx5 = find_nearest(X_, nx3)
        
        # Capture the index in the spline of nearest point
        yln = Y_[idx]
        y_2 = Y_[idx2]
        y_3 = Y_[idx3]
        y_4 = Y_[idx4]
        y_5 = Y_[idx5]
        y_6 = min(Y_)
        ax.scatter(xln, yln, color='k', marker='o', alpha=0.8)
        ax.scatter(x_2, y_2, color='b', marker='o', alpha=0.8) #2x
        
        # Calculate relative MSEs for stacked bar plots
        d.append(
            {
                'nlnn': yln - y_6,
                '3n': y_5 - y_6,
                '2n': y_2 - y_6,
                '1.5n': y_4 - y_6,
                '1n': y_3 - y_6
            }
        )
        df = pd.DataFrame(d)
        ind = ind + 1
        
        # Plot smoothed plot
        ax.plot(X_, Y_, label=label, color = color)
    return df
    

###################################
# Plotting of curves
###################################
# plot raw data 1
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=False, figsize=(14, 5))

# Plot.
df_A = curve_plot(ax1, A_splits, "A")
df_D = curve_plot(ax2, D_splits, "D")
df_rand = curve_plot(ax3, rand_splits, "Random")

# Format plots.
ax1.set_title('A-optimal')
ax2.set_title('D-optimal')
ax3.set_title('Random')
ax3.legend( loc="upper right")

ax1.set_xlabel("number of edges, k")
ax2.set_xlabel("number of edges, k")
ax3.set_xlabel("number of edges, k")
ax1.set_ylabel("<MSE \u0394G> (kcal/mol)$^2$")

#ax1.set_ylim(0.5, 1.0)
#ax2.set_ylim(0.5, 1.0)
#ax3.set_ylim(0.5, 1.0)

# Display, save.
plt.savefig('cycles_surface.png', dpi=500)

######################################################################
# Plotting of stacked bar,
# this section requires removing commenting to plot each stacked bar plot
######################################################################

import random
# Get the pandas dataframe output
#df = pd.DataFrame(np.abs(np.random.randn(10,10)),columns=['A','B','C','D','E','F','G','H','I','J'], index=range(10))
colors2 = cmap(np.linspace(0, 0.8, 5))

#plt.figure(1)
#plt.subplot(1,1,1)

### Plot the stacked bar plots for A-optimal:

#df_A.plot( kind='bar', stacked=True, color=colors2, logy=False, legend=False) #no need to specify for first axis
#plt.xlabel("number of ligands, n")
#plt.ylabel("relative <MSE \u0394G> (kcal/mol)$^2$")
#plt.title("A-optimal")

### Plot the stacked bar plots for D-optimal:

#plt.subplot(1,2,1)
#df_D.plot( kind='bar', stacked=True, color=colors2, logy=False, legend=False)
#plt.xlabel("number of ligands, n")
#plt.title("D-optimal")

### Plot the stacked bar plots for Random:

#plt.subplot(1,3,1)
df_rand.plot( kind='bar', stacked=True, color=colors2, logy=False)
plt.xlabel("number of ligands, n")
plt.title("Pseudo-random")
plt.savefig('stacked_random.pdf', dpi=500)

plt.show()



