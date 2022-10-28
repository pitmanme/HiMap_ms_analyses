# Code author: Mary Pitman 2022; mpitman@uci.edu

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import matplotlib.cm as cm
import re

from numpy import array
import itertools
import matplotlib.mlab as mlab
import seaborn as sns
import pandas as pd
#-------------------------------------------------------------------------------------#
# Global font controls
font = {'family' : 'arial'}
# Global sns controls
sns.set(style="white", palette="muted", color_codes=True)

y1 = []
y2 = []
y3 = []
y1b = []
y2b = []
y3b = []
y1c = []
y2c = []
y3c = []

#File to load
fobj1 = np.loadtxt('leverages_1.txt', comments='#')
fobj2 = np.loadtxt('leverages_2.txt', comments='#')
fobj3 = np.loadtxt('leverages_3.txt', comments='#')
fobj4 = np.loadtxt('leverages_4.txt', comments='#')
fobj5 = np.loadtxt('leverages_5.txt', comments='#')

#selection of data ranges from file
# Plot 1
y1=fobj1[:,0]
y2=fobj1[:,1]
y3=fobj1[:,2]
av1 = np.average(y1)

# Plot 2
y1b=fobj2[:,0]
y2b=fobj2[:,1]
y3b=fobj2[:,2]
av2 = np.average(y1b)

#Plot 3
y1c=fobj3[:,0]
y2c=fobj3[:,1]
y3c=fobj3[:,2]
av3 = np.average(y1c)

#Plot 4
y1d=fobj4[:,0]
y2d=fobj4[:,1]
y3d=fobj4[:,2]
av4 = np.average(y1d)

#Plot 5
y1e=fobj5[:,0]
y2e=fobj5[:,1]
y3e=fobj5[:,2]
av5 = np.average(y1e)

# Three subplots, #4878CF color #8EC4DE 
f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(ncols=5, sharey=True)


sns.kdeplot(y1, shade=True, clip = (0, 1),  ax=ax1, color='#fc9112');
sns.kdeplot(y2, shade=True, clip = (0, 1), ax=ax1, color='#3e007f');
sns.kdeplot(y3, shade=True, clip = (0, 1), ax=ax1, color='#db3a3f');

sns.kdeplot(y1b, shade=True, clip = (0, 1),  ax=ax2, color='#fc9112');
sns.kdeplot(y2b, shade=True, clip = (0, 1), ax=ax2, color='#3e007f');
sns.kdeplot(y3b, shade=True, clip = (0, 1), ax=ax2, color='#db3a3f');

sns.kdeplot(y1c, shade=True, clip = (0, 1), ax=ax3, color='#fc9112');
sns.kdeplot(y2c, shade=True, clip = (0, 1), ax=ax3, color='#3e007f');
sns.kdeplot(y3c, shade=True, clip = (0, 1), ax=ax3, color='#db3a3f');

sns.kdeplot(y1d, shade=True, clip = (0, 1), ax=ax4, color='#fc9112');
sns.kdeplot(y2d, shade=True, clip = (0, 1), ax=ax4, color='#3e007f');
sns.kdeplot(y3d, shade=True, clip = (0, 1), ax=ax4, color='#db3a3f');

sns.kdeplot(y1e, shade=True, clip = (0, 1), ax=ax5, color='#fc9112');
sns.kdeplot(y2e, shade=True, clip = (0, 1), ax=ax5, color='#3e007f');
sns.kdeplot(y3e, shade=True, clip = (0, 1), ax=ax5, color='#db3a3f');

#plot verticle average lines
ax1.axvline(x=av1, color='k', linestyle='--', alpha = 0.8)
ax2.axvline(x=av2, color='k', linestyle='--', alpha = 0.8)
ax3.axvline(x=av3, color='k', linestyle='--', alpha = 0.8)
ax4.axvline(x=av4, color='k', linestyle='--', alpha = 0.8)
ax5.axvline(x=av5, color='k', linestyle='--', alpha = 0.8)

#Label the x and y axis
ax1.set_xlabel('Leverage', fontsize=16)
ax2.set_xlabel('Leverage', fontsize=16)
ax3.set_xlabel('Leverage', fontsize=16)
ax4.set_xlabel('Leverage', fontsize=16)
ax5.set_xlabel('Leverage', fontsize=16)

ax1.set_ylabel('Density', fontsize=16)
ax2.set_ylabel(' ', fontsize=16)
ax3.set_ylabel(' ', fontsize=16)
ax4.set_ylabel(' ', fontsize=16)
ax5.set_ylabel(' ', fontsize=16)


#Set titles for each subplot
ax1.set_title("cdk2", fontsize=16, color='#5A5A5A')
ax2.set_title("tyk2", fontsize=16, color='#5A5A5A')
ax3.set_title("tnks2", fontsize=16, color='#5A5A5A')
ax4.set_title("ptp1b", fontsize=16, color='#5A5A5A')
ax5.set_title("pfkfb3", fontsize=16, color='#5A5A5A')

#Define x-axis ranges
ax1.set_xlim([0,1])
ax2.set_xlim([0,1])
ax3.set_xlim([0,1])
ax4.set_xlim([0,1])
ax5.set_xlim([0,1])

ax1.set_ylim([0,5.1])
ax2.set_ylim([0,5.1])
ax3.set_ylim([0,5.1])
ax4.set_ylim([0,5.1])
ax5.set_ylim([0,5.1])

# Axis tick controls
for ax in [ax1, ax2, ax3, ax4, ax5]:
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)

#Set legend labels on plot
ax5.text(0.05, 4.8, 'A-opt', ha='left', va='center',  fontsize=16, color='#fc9112')
ax5.text(0.05, 4.4, 'D-opt', ha='left', va='center', fontsize=16, color='#3e007f')
ax5.text(0.05, 4.0, 'LOMAP', ha='left', va='center',  fontsize=16, color='#db3a3f')

# Fine-tune figure;
f.tight_layout()
f.subplots_adjust(wspace=0.40)

#f.savefig('test2png.png', dpi=300)
plt.show()
