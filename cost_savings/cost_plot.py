# Code author: Mary Pitman 2022; mpitman@uci.edu

# Import libraries
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

from scipy import interpolate
from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.colors as colors
import matplotlib.ticker as ticker

plt.rcParams.update({'font.size': 10})

#---------------------------------------
# Functions
#---------------------------------------
def min_cost(x, y):
    n_nodes = (x/y)
    n_edges = n_nodes*np.log(n_nodes)
    cost = n_edges*24.0
    return cost
    
def max_cost(x, y):
    # All clusters are size 2 except one cluster
    # Number cluster with only 2 ligands
    small_cs = (y-1.0)*2.0
    # Get number of smaller clusters
    num_small = y - 1.0
    # Number of edges from small clusters
    edges_small = num_small * 2.0 * np.log(2.0)
    # Remaining nodes
    remaining = x - small_cs
    edges_large = remaining * np.log(remaining)
    # Total number of edges
    n_edges = edges_small + edges_large
    cost = n_edges*24.0
    return cost
    
#---------------------------------------
# Script
#---------------------------------------
fig = plt.figure(figsize =(6, 10))
ax = plt.subplot(211, projection='3d')
ax2 = plt.subplot(212, projection='3d')

# Generate data points for mesh
x = np.arange(2.0, 150.0, 0.5)
y = x/2.0

X, Y = np.meshgrid(x, y)

# Calculate z values.
# Min.
zs = np.array(min_cost(X, Y))
Z = zs.reshape(X.shape)
# Max.
zs2 = np.array(max_cost(X, Y))
Z2 = zs2.reshape(X.shape)

# Creating color map
my_cmap = plt.get_cmap('inferno')

# Change color bar to max Z of all Z
# Plot.
surf = ax.plot_surface(X, Y, Z, cmap = my_cmap, norm=colors.LogNorm(vmin=1, vmax=Z.max()), edgecolor ='none', linewidth=0, antialiased=False, rcount=200, ccount=200)

surf2 = ax2.plot_surface(X, Y, Z2, cmap = my_cmap, norm=colors.LogNorm(vmin=1, vmax=Z.max()), edgecolor ='none', linewidth=0, antialiased=False, rcount=200, ccount=200)

# Colorbar controls
axs =[ax, ax2]
cbar = fig.colorbar(surf, ax = axs,  location='bottom', shrink= 0.7)
cbar.set_label('Cost (USD)')

# Plotting aesthetic controls.

# Set axis labels
ax.set_xlabel('N ligands', labelpad=9)
ax.set_ylabel('N clusters', labelpad=9)
ax.set_zlabel('Cost (USD)', labelpad=9)

ax2.set_xlabel('N ligands', labelpad=9)
ax2.set_ylabel('N clusters', labelpad=9)
ax2.set_zlabel('Cost (USD)', labelpad=9)

# Subplot title
ax.set_title('Minimum')
ax2.set_title('Maximum')

# Control bounds
lig_max = np.max(x)

# Axis wall colors.
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.04))
ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.09))

ax2.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax2.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.04))
ax2.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.09))

# Control tick spacking
xtick_spacing = 50
ax.xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))
ax2.xaxis.set_major_locator(ticker.MultipleLocator(xtick_spacing))

ytick_spacing = 20
ax.yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))
ax2.yaxis.set_major_locator(ticker.MultipleLocator(ytick_spacing))

ztick_spacing = 4000
ax.zaxis.set_major_locator(ticker.MultipleLocator(ztick_spacing))
ax2.zaxis.set_major_locator(ticker.MultipleLocator(ztick_spacing))

# Control view.
ax.view_init(25, 160)
ax2.view_init(25, 160)


plt.savefig('cost_savings.png', dpi=500)
plt.show()



