#!/usr/bin/env python2

import numpy as np
from matplotlib import pyplot as plt

import plot
import plot_params

linewidth = 3

s = 1.0     # Vertices size
d = 0.25    # Distance from two planes
diag = np.sqrt(2.0*d**2)


# Vertices:
#       6-------7
#     / |     / |
#   2-------3   |
#   |   |   |   |
#   |   4---|---5   -
#   | /     | /     | d
#   0-------1       -
#
#           |---|
#             d

# [x,y]
vertices = np.asarray([
                [0.0,   0.0],   # 0
                [s,     0.0],   # 1
                [0.0,   s],     # 2
                [s,     s],     # 3
                [d,    d],      # 4
                [s+d,  d],      # 5
                [d,    s+d],    # 6
                [s+d,  s+d]     # 7
            ])

edges_full = np.asarray([
                [0, 1],
                [1, 3],
                [3, 2],
                [3, 7],
                [2, 0],
                [1, 5],
                [5, 7],
                [7, 6],
                [6, 2]
            ])
edges_shaded = np.asarray([
                [0, 4],
                [4, 5],
                [4, 6]
            ])

def plot_cell(ax):
    for ci, edge in enumerate(edges_full):
        c0 = edge[0]
        c1 = edge[1]
        ax.plot([vertices[c0][0], vertices[c1][0]], [vertices[c0][1], vertices[c1][1]], '-ok', linewidth = linewidth, markersize = 12)
    for ci, edge in enumerate(edges_shaded):
        c0 = edge[0]
        c1 = edge[1]
        ax.plot([vertices[c0][0], vertices[c1][0]], [vertices[c0][1], vertices[c1][1]], '-ok', linewidth = linewidth, markersize = 12, alpha = 0.4)

    # Axis
    al = d / 1.0 # axis length
    al2 = al/(np.sqrt(2.0)*1.5)
    plt.arrow(d, d, 0.0,  al,   shape='full', lw=3, color = 'k', length_includes_head = True, head_width = 0.02)
    plt.arrow(d, d, al,   0.0,  shape='full', lw=3, color = 'k', length_includes_head = True, head_width = 0.02)
    plt.arrow(d, d, -al2, -al2, shape='full', lw=3, color = 'k', length_includes_head = True, head_width = 0.02)
    ax.text(1.7*d, 1.05*d,'x')
    ax.text(0.8*d, 1.6*d, 'y')
    ax.text(0.8*d, 0.6*d, 'z')

    # Corners
    g = d / 8.0 # Gap
    ax.text(-g, 0,      '(i,j,k+1)',    verticalalignment = 'center', horizontalalignment = 'right') # 0
    ax.text(s+g,0,      '(i+1,j,k+1)',  verticalalignment = 'center', horizontalalignment = 'left')  # 1
    ax.text(-g, s,      '(i,j+1,k+1)',  verticalalignment = 'center', horizontalalignment = 'right') # 2
    ax.text(s+g,s,      '(i+1,j+1,k+1)',verticalalignment = 'center', horizontalalignment = 'left')  # 3
    ax.text(d-g, d,      '(i,j,k)',     verticalalignment = 'center', horizontalalignment = 'right') # 4
    ax.text(d+s+g,d,    '(i+1,j,k+1)',  verticalalignment = 'center', horizontalalignment = 'left')  # 5
    ax.text(d-g, s+d,   '(i,j+1,k)',    verticalalignment = 'center', horizontalalignment = 'right') # 6
    ax.text(s+d+g, s+d, '(i+1,j+1,k)',  verticalalignment = 'center', horizontalalignment = 'left')  # 7

# ******************************************************************************************
# Yee cell (FDTD)
fig_yee = plot.figure()
ax = fig_yee.add_subplot(1,1,1, aspect = 'equal', xticks=[], yticks=[], frameon=False)

plot_cell(ax)

# H arrows              x           y           dx      dy
ax.add_patch(plt.Arrow((s+d)/2.0,   d/2.0,      0.0,    d,      width = d/2.0, edgecolor='none', facecolor = 'blue'))
ax.add_patch(plt.Arrow(d/2.0,       (s+d)/2.0,  d,      0.0,    width = d/2.0, edgecolor='none', facecolor = 'blue'))
ax.add_patch(plt.Arrow(d+s/2.0,     d+s/2.0,    -d/3.0, -d/3.0, width = d/2.0, edgecolor='none', facecolor = 'blue'))
# H text
ax.text(d/2.0,       (s+d)/2.0,     '$H_{x(i,j+1/2,k+1/2)}$', color = 'blue', verticalalignment = 'top',    horizontalalignment = 'center')
ax.text((s+d)/2.0,   d/2.0,         '$H_{y(i+1/2,j,k+1/2)}$', color = 'blue', verticalalignment = 'top',    horizontalalignment = 'center')
ax.text(d+s/2.0,   d+s/2.0,         '$H_{z(i+1/2,j+1/2,k)}$', color = 'blue', verticalalignment = 'bottom', horizontalalignment = 'center')

# E arrows              x           y           dx      dy
ax.add_patch(plt.Arrow(d+s/2.0,     d,          d,      0.0,    width = d/2.0, edgecolor='none', facecolor = 'red'))
ax.add_patch(plt.Arrow(d,           d+s/2.0,    0,      d,      width = d/2.0, edgecolor='none', facecolor = 'red'))
ax.add_patch(plt.Arrow(d/2.0,       d/2.0,      -d/3.0, -d/3.0, width = d/2.0, edgecolor='none', facecolor = 'red'))
# E text
ax.text(d+s/2.0,    d,          '$E_{x(i+1/2,j,k)}$', color = 'red', verticalalignment = 'bottom', horizontalalignment = 'left')
ax.text(d,          d+s/2.0,    '$E_{y(i,j+1/2,k)}$', color = 'red', verticalalignment = 'bottom', horizontalalignment = 'right')
ax.text(d/3.0,      d/3.0,      '$E_{z(i,j,k+1/2)}$', color = 'red', verticalalignment = 'top',    horizontalalignment = 'left')

ax.set_xlim((-d/7.0, s+d+d/7.0))
ax.set_ylim((-d/7.0, s+d+d/7.0))

# ******************************************************************************************
# QFDTD cell
fig_qfdtd = plot.figure()
ax = fig_qfdtd.add_subplot(1,1,1, aspect = 'equal', xticks=[], yticks=[], frameon=False)

plot_cell(ax)

ax.plot([d], [d], 'o', markeredgecolor = 'green', markerfacecolor = 'green', markersize = 14)
ax.text(1.1*d, 0.9*d, '$\psi_{i,j,k}$', color = 'green', verticalalignment = 'top', horizontalalignment = 'left')

ax.set_xlim((-d/7.0, s+d+d/7.0))
ax.set_ylim((-d/7.0, s+d+d/7.0))

for ext in ['pdf', 'svg']:
    plot.savefig(['fdtd_cell_yee.' + ext, 'fdtd_cell_qfdtd.' +  ext])
plot.show()
