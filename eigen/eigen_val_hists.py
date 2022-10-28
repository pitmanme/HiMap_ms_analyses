# Code author: Mary Pitman 2022; mpitman@uci.edu

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import matplotlib.cm as cm
import re
#from StringIO import StringIO
from numpy import array
import itertools
import matplotlib.mlab as mlab
import seaborn as sns
import pandas as pd
#-------------------------------------------------------------------------------------#
#Global font controls
font = {'family' : 'arial'}

#global sns control
sns.set(style="white", palette="muted", color_codes=True)

#File to load
fobj1 = np.loadtxt('jdata_n_20.txt', comments='#')

# Parse input data file.
# Histogram 1
y1a=fobj1[:, 0]
y2a=fobj1[:, 1]
y3a=fobj1[:, 2]
y1 = y1a + y2a + y3a

# Histogram 2
y1b=fobj1[:, 3]
y2b=fobj1[:, 4]
y3b=fobj1[:, 5]
y2 = y1b + y2b + y3b

# Histogram 3
y1c=fobj1[:, 6]
y2c=fobj1[:, 7]
y3c=fobj1[:, 8]
y3 = y1c + y2c + y3c

# Histogram 4
y1d=fobj1[:, 9]
y2d=fobj1[:, 10]
y3d=fobj1[:, 11]
y4 = y1d + y2d + y3d

# Calculate the average determinants (product of eigenvalues)
av1 = (np.prod(y1a) + np.prod(y1b) + np.prod(y1c))/3.0
av2 = (np.prod(y1b) + np.prod(y2b) + np.prod(y3b))/3.0
av3 = (np.prod(y1c) + np.prod(y2c) + np.prod(y3c))/3.0
av4 = (np.prod(y1d) + np.prod(y2d) + np.prod(y3d))/3.0

# Set to no subplots
f, ax1 = plt.subplots(ncols=1, sharey=True)

# Plot.
sns.kdeplot(y1, shade=True, clip = (0, 20), ax=ax1, color='#fca612', label='k = 20', bw_adjust = 0.5);
sns.kdeplot(y2, shade=True, clip = (0, 20), ax=ax1, color='#cf7202', label='k = 40', bw_adjust = 0.5);
sns.kdeplot(y3, shade=True, clip = (0, 20), ax=ax1, color='#824801', label='k = 60', bw_adjust = 0.5);
sns.kdeplot(y4, shade=True, clip = (0, 20), ax=ax1, color='#331c01', label='k = 80', bw_adjust = 0.5);

#Plotting controls
axes = [ax1]
for ax in axes:
    ax.set_xlabel('Eigenvalues of $\mathscr{F}$', fontsize=20)
    ax.set_ylabel(' ', fontsize=20)
    ax.set_xlim([0,20])
    ax.set_ylim([0,0.35])
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)

ax1.set_ylabel('Probability Density', fontsize=20)

# Include legend
ax1.legend(prop={'size': 20})

# Fine-tune figure;
f.tight_layout()
f.subplots_adjust(wspace=0.20)
f.set_size_inches(6.5, 6.5)
f.savefig('eigen_hists.png', dpi=400)
plt.show()
