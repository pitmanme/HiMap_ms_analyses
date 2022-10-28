# Code author: Mary Pitman 2022; mpitman@uci.edu

# imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib import tri, cm

from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.path import Path
import matplotlib.patches as patches
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


###################################
# File import
###################################
fobj = np.loadtxt('cycle_sim_data.txt', comments='#')

#Parse file
n = fobj[:,0]   # Y valz
k = fobj[:,1]   # X vals
A = fobj[:,2] *10**-4  # Z's
D = fobj[:,3] *10**-4
rand = fobj[:,4]*10**-4

ln_k = np.log(k)
ln_n = np.log(n)
ln_A = np.log(A)
ln_D = np.log(D)
ln_rand = np.log(rand)


# Do triangular interpolation 1
triang = tri.Triangulation(k, n)
refiner = tri.UniformTriRefiner(triang)
new, new_z = refiner.refine_field(rand, subdiv=4)

norm = plt.Normalize(vmax=abs(D).max(), vmin=abs(D).min()-0.6)

# Do triangular interpolation 2
new2, new_z2 = refiner.refine_field(D, subdiv=4)
kwargs = dict(triangles=new.triangles, norm=norm, cmap='inferno', linewidth=0.2, alpha=1)


# Initiate subplots
fig = plt.figure(figsize=(12, 5))
ax = plt.subplot(121, projection='3d')
ax2 = plt.subplot(122, projection='3d')

# Plot datapoints
ax.scatter(k, n, rand, color='k', marker='o', alpha=1.0)
ax2.scatter(k, n, D, color='k', marker='o', alpha=1.0)

# Get axis limit variables
xlim = ax2.get_xlim()
ylim = ax2.get_ylim()
zlim = ax2.get_zlim()


###################################
# Do 3D plane plotting with projection
###################################
# Generate surfaces.
surf = ax.plot_trisurf(new.x, new.y, new_z, **kwargs)
surf2 = ax2.plot_trisurf(new2.x, new2.y, new_z2, **kwargs)
# Set plot limits.
ax.set_zlim(np.min(new_z2), np.max(new_z2))
ax2.set_zlim(np.min(new_z2), np.max(new_z2))

# Refine axis wall colors. Makes less dark than default.
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.04))
ax.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.09))

ax2.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax2.w_yaxis.set_pane_color((0.0, 0.0, 0.0, 0.04))
ax2.w_zaxis.set_pane_color((0.0, 0.0, 0.0, 0.09))

# Set labels.
ax.set_xlabel('k, edges')
ax.set_ylabel('n, nodes')
ax.set_zlabel('$10^4$ cycles')

ax2.set_xlabel('k, edges')
ax2.set_ylabel('n, nodes')

# Titles.
ax.set_title('Random')
ax2.set_title('D-optimal')

# Color bar controls.
#cbar = plt.colorbar(surf2, shrink=0.6, aspect=15)
#cbar.set_label('$10^4$ cycles', labelpad = 20.0, rotation=270)

# Control view.
ax.view_init(10, 230)
ax2.view_init(10, 230)

plt.savefig('cycles_surface_rand.png', dpi=500)
plt.show()
