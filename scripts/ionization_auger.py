#!/usr/bin/env python2

import numpy as np
from matplotlib import pyplot as plt

import matplotlib.patches as mpatches

import plot
import plot_params
from fleches import *
from photon import *

r_core  = 1.0
r_e     = 0.25

el_N      = np.array([2, 8, 18, 18, 8])
el_r      = r_core * np.arange(2.0, len(el_N)+2)

photon_length = 4.0*r_core

# Inner shell ionization: 4th n = 4 shell electron absorbs a photon
auger_1_n = 4 # Principal number
auger_1_i = 4 # Electron id
auger_1_xy = np.array([None]*2)

# Outer shell decay: 3rd n = 5 shell electron
auger_2_n = 5 # Principal number
auger_2_i = 3 # Electron id
auger_2_xy = np.array([None]*2)

# Outer shell ionization: 4th n = 5 shell electron
auger_3_n = 5 # Principal number
auger_3_i = 4 # Electron id
auger_3_xy = np.array([None]*2)

#fig = plot.figure()
#ax1 = fig.add_subplot(1,2,1, aspect = 'equal', xticks=[], yticks=[], frameon=False)
#ax2 = fig.add_subplot(1,2,2, aspect = 'equal', xticks=[], yticks=[], frameon=False)
fig1 = plot.figure()
fig2 = plot.figure()
ax1 = fig1.add_subplot(1,1,1, aspect = 'equal', xticks=[], yticks=[], frameon=False)
ax2 = fig2.add_subplot(1,1,1, aspect = 'equal', xticks=[], yticks=[], frameon=False)


# *******************************************************************************************
# Generic

# Atom core
for ax in [ax1, ax2]:
    ax.add_patch(mpatches.Circle([0,0], r_core, color = 'm', ec="none", alpha=0.6))
    ax.text(0, 0, "Xe", ha="center", va="center", size=20)

    # Loop over principal numbers (orbitals)
    for si in xrange(len(el_N)):
        # Circle around core
        ax.add_patch(mpatches.Circle([0,0], el_r[si], color='none', ec="black", alpha=0.4))
        # Orbital principal number
        theta = 5.3*np.pi/32.0
        ax.text(el_r[si]*np.cos(theta), el_r[si]*np.sin(theta), str(si+1),
                ha="center", va="center",
                size=20, backgroundcolor = 'w')
        # Every electrons
        dtheta = 2.0*np.pi / float(el_N[si])
        for ni in xrange(el_N[si]):
            x = el_r[si] * np.cos(dtheta*float(ni))
            y = el_r[si] * np.sin(dtheta*float(ni))
            # Store Auger electrons position
            if ((si == auger_1_n-1) and (ni == auger_1_i-1)):
                auger_1_xy[0] = x
                auger_1_xy[1] = y
            elif ((si == auger_2_n-1) and (ni == auger_2_i-1)):
                auger_2_xy[0] = x
                auger_2_xy[1] = y
            elif ((si == auger_3_n-1) and (ni == auger_3_i-1)):
                auger_3_xy[0] = x
                auger_3_xy[1] = y
            else:
                ax.add_patch(mpatches.Circle([x,y], r_e, color='blue', ec="none", alpha=0.6))


# *******************************************************************************************
# Auger 1st step

# Outer shells still intact
ax1.add_patch(mpatches.Circle(auger_2_xy, r_e, color='blue', ec="none", alpha=0.6))
ax1.add_patch(mpatches.Circle(auger_3_xy, r_e, color='blue', ec="none", alpha=0.6))

# New hole in inner shell filling
ax1.add_patch(mpatches.Circle(auger_1_xy, r_e, color='none', ec="blue", lw = 2, ls = 'dotted'))

# Photon hitting inner shell
photon(ax1, [auger_1_xy[0]-photon_length, auger_1_xy[1]], auger_1_xy, 1.0)

# Inner shell electron leaving
auger_1_xy_new = np.array([None]*2)
angle = np.arctan2(auger_1_xy[1], auger_1_xy[0]) * 0.8
auger_1_xy_new[0] = 1.3 * el_r.max() * np.cos(angle)
auger_1_xy_new[1] = 1.3 * el_r.max() * np.sin(angle)
ax1.add_patch(mpatches.Circle(auger_1_xy_new, r_e, color='blue', ec="none", alpha=1.0))

# Arrow showing it
ar = arrow('AuE', auger_1_xy, auger_1_xy_new)    # Vertical (Ke)
ar.Plot(ax1, color = 'g')


# *******************************************************************************************
# Auger 2nd step

# New outer shell holes
ax2.add_patch(mpatches.Circle(auger_2_xy, r_e, color='none', ec="blue", lw = 2, ls = 'dotted'))
ax2.add_patch(mpatches.Circle(auger_3_xy, r_e, color='none', ec="blue", lw = 2, ls = 'dotted'))

# Inner shell filling
ax2.add_patch(mpatches.Circle(auger_1_xy, r_e, color='blue', ec="none", alpha = 1.0))

# Arrow for inner shell filling
arrow_end   = (auger_1_xy[0], auger_1_xy[1])
arrow_start = (auger_2_xy[0], auger_2_xy[1])
ax2.annotate('', xycoords='data',
             xy = arrow_end, xytext = arrow_start, textcoords='data',
             size = 20, arrowprops=dict(arrowstyle="simple",
                                        fc="g", ec="none",
                                        connectionstyle="arc3,rad=-0.3"))

# Outer shell electron leaving
auger_3_xy_new = np.array([None]*2)
angle = np.arctan2(auger_3_xy[1], auger_3_xy[0])
auger_3_xy_new[0] = 1.4 * el_r.max() * np.cos(angle)
auger_3_xy_new[1] = 1.4 * el_r.max() * np.sin(angle)
ax2.add_patch(mpatches.Circle(auger_3_xy_new, r_e, color='blue', ec="none", alpha=1.0))

# Arrow showing it
ar = arrow('AuE', auger_3_xy, auger_3_xy_new)    # Vertical (Ke)
ar.Plot(ax2, color = 'g')

# *******************************************************************************************
for ax in [ax1, ax2]:
    ax.set_xlim((-1.1*el_r.max(), 1.1*el_r.max()))
    ax.set_ylim((-1.1*el_r.max(), 1.1*el_r.max()))

for ext in ['pdf', 'svg']:
    plot.savefig(['auger_step_1.' + ext, 'auger_step_2.' +  ext])
plot.show()
