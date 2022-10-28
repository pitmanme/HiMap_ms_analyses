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
#Global font cenpaols

font = {'family' : 'arial'}
#matplotlib.rc('font',**font)

#global sns cenpaol
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

fobj1 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_1.txt', comments='#')
fobj2 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_2.txt', comments='#')
fobj3 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_3.txt', comments='#')
fobj4 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_4.txt', comments='#')
fobj5 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_5.txt', comments='#')

fobj6 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_6.txt', comments='#')
fobj7 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_7.txt', comments='#')
fobj8 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_8.txt', comments='#')
fobj9 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_9.txt', comments='#')
fobj10 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_10.txt', comments='#')

fobj11 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_11.txt', comments='#')
fobj12 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_12.txt', comments='#')
fobj13 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_13.txt', comments='#')
fobj14 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_14.txt', comments='#')
fobj15 = np.loadtxt('/Users/mpitman/work/gits/packages/Lomap-main/examples/paper_analysis/leverages/levs_A_15.txt', comments='#')


#selection of data ranges from file
# Plot 1
y1=fobj1[:]
y2=fobj2[:]
y3=fobj3[:]
av1 = np.average(y1)
av12 = np.average(y2)
av13 = np.average(y3)

# Plot 2
y1b=fobj4[:]
y2b=fobj5[:]
y3b=fobj6[:]
av2 = np.average(y1b)
av22 = np.average(y2b)
av23 = np.average(y3b)

#Plot 3
y1c=fobj7[:]
y2c=fobj8[:]
y3c=fobj9[:]
av3 = np.average(y1c)
av32 = np.average(y2c)
av33 = np.average(y3c)

#Plot 4
y1d=fobj10[:]
y2d=fobj11[:]
y3d=fobj12[:]
av4 = np.average(y1d)
av42 = np.average(y2d)
av43 = np.average(y3d)

#Plot 5
y1e=fobj13[:]
y2e=fobj14[:]
y3e=fobj15[:]
av5 = np.average(y1e)
av52 = np.average(y2e)
av53 = np.average(y3e)

# Three subplots, #4878CF color #8EC4DE 
f, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(ncols=5, sharey=True)
#fig, axes = plt.subplots(ncols=5, sharex=False, sharey=True)
#sns.despine()

sns.kdeplot(y1, shade=True, clip = (0, 1),  ax=ax1, color='#fc9112');
sns.kdeplot(y2, shade=True, clip = (0, 1), ax=ax1, color='#b16202');
sns.kdeplot(y3, shade=True, clip = (0, 1), ax=ax1, color='#4c2a01');

sns.kdeplot(y1b, shade=True, clip = (0, 1),  ax=ax2, color='#fc9112');
sns.kdeplot(y2b, shade=True, clip = (0, 1), ax=ax2, color='#b16202');
sns.kdeplot(y3b, shade=True, clip = (0, 1), ax=ax2, color='#4c2a01');

sns.kdeplot(y1c, shade=True, clip = (0, 1), ax=ax3, color='#fc9112');
sns.kdeplot(y2c, shade=True, clip = (0, 1), ax=ax3, color='#b16202');
sns.kdeplot(y3c, shade=True, clip = (0, 1), ax=ax3, color='#4c2a01');

sns.kdeplot(y1d, shade=True, clip = (0, 1), ax=ax4, color='#fc9112');
sns.kdeplot(y2d, shade=True, clip = (0, 1), ax=ax4, color='#b16202');
sns.kdeplot(y3d, shade=True, clip = (0, 1), ax=ax4, color='#4c2a01');

sns.kdeplot(y1e, shade=True, clip = (0, 1), ax=ax5, color='#fc9112', label='control');
sns.kdeplot(y2e, shade=True, clip = (0, 1), ax=ax5, color='#b16202', label='2n');
sns.kdeplot(y3e, shade=True, clip = (0, 1), ax=ax5, color='#4c2a01', label='n$\cdot$ln(n)');

#plot verticle average lines
ax1.axvline(x=av1, color='#fc9112', linestyle='--', alpha = 0.6)
ax1.axvline(x=av12, color='#b16202', linestyle='--', alpha = 0.6)
ax1.axvline(x=av13, color='#4c2a01', linestyle='--', alpha = 0.6)

ax2.axvline(x=av2, color='#fc9112', linestyle='--', alpha = 0.6)
ax2.axvline(x=av22, color='#b16202', linestyle='--', alpha = 0.6)
ax2.axvline(x=av23, color='#4c2a01', linestyle='--', alpha = 0.6)

ax3.axvline(x=av3, color='#fc9112', linestyle='--', alpha = 0.6)
ax3.axvline(x=av32, color='#b16202', linestyle='--', alpha = 0.6)
ax3.axvline(x=av33, color='#4c2a01', linestyle='--', alpha = 0.6)

ax4.axvline(x=av4, color='#fc9112', linestyle='--', alpha = 0.6)
ax4.axvline(x=av42, color='#b16202', linestyle='--', alpha = 0.6)
ax4.axvline(x=av43, color='#4c2a01', linestyle='--', alpha = 0.6)

ax5.axvline(x=av5, color='#fc9112', linestyle='--', alpha = 0.6)
ax5.axvline(x=av52, color='#b16202', linestyle='--', alpha = 0.6)
ax5.axvline(x=av53, color='#4c2a01', linestyle='--', alpha = 0.6)

#Plotting controls
axes = [ax1, ax2, ax3, ax4, ax5]
for ax in axes:
    ax.set_xlabel('Leverage', fontsize=16)
    ax.set_ylabel(' ', fontsize=16)
    ax.set_xlim([0,1])
    ax.set_ylim([0,9.5])
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)

ax1.set_ylabel('Density', fontsize=16)

#Set titles for each subplot
ax1.set_title("cdk2", fontsize=16, color='#5A5A5A')
ax2.set_title("tyk2", fontsize=16, color='#5A5A5A')
ax3.set_title("tnks2", fontsize=16, color='#5A5A5A')
ax4.set_title("ptp1b", fontsize=16, color='#5A5A5A')
ax5.set_title("pfkfb3", fontsize=16, color='#5A5A5A')

ax5.legend()


# Fine-tune figure;
f.tight_layout()
f.subplots_adjust(wspace=0.40)
f.set_size_inches(12, 4)
f.savefig('testA.png', dpi=300)
#plt.show()
