#!/usr/bin/env python2

import numpy as np
from matplotlib import pyplot as plt

import plot
import plot_params



np.random.seed(0)

rmax = 1.0
N = 300
max_level = 4.0

r       = rmax * np.random.randn(N,1)
theta   = 2.0 * np.pi * np.random.random((N, 1))

linewidth = 2.0

class Part:
    def __init__(self, x, y):
        self.x = np.asarray(x)
        self.y = np.asarray(y)
    def is_inside(self, cell):
        indices = np.nonzero(  (cell.mins[0] < self.x) & (self.x <= cell.maxs[0])
                             & (cell.mins[1] < self.y) & (self.y <= cell.maxs[1]))[0]
        return (len(indices) > 0)

class Cell:
    def __init__(self, x, y, sizes, level = 0):
        self.pos   = np.asarray([x, y])
        self.sizes = np.asarray(sizes)
        self.mins  = self.pos - self.sizes/2.0
        self.maxs  = self.pos + self.sizes/2.0
        self.level = level
    def draw_quadtree(self, ax, particles):
        rect = plt.Rectangle(self.mins, *self.sizes, edgecolor='black', facecolor='none', linewidth = linewidth)
        ax.add_patch(rect)
        ax.plot([self.mins[0], self.maxs[0]], self.pos[1]*np.ones(2), '-', color = 'black', linewidth = linewidth)
        ax.plot(self.pos[0]*np.ones(2), [self.mins[1], self.maxs[1]], '-', color = 'black', linewidth = linewidth)
        self.draw_daugthers(ax, particles)
    def draw_daugthers(self, ax, particles):
        if particles.is_inside(self):
            # Draw 4 daugthers
            daugthers = [Cell(self.pos[0] - self.sizes[0]/4.0, self.pos[1] - self.sizes[1]/4.0, self.sizes/2.0, level = self.level + 1),
                         Cell(self.pos[0] - self.sizes[0]/4.0, self.pos[1] + self.sizes[1]/4.0, self.sizes/2.0, level = self.level + 1),
                         Cell(self.pos[0] + self.sizes[0]/4.0, self.pos[1] - self.sizes[1]/4.0, self.sizes/2.0, level = self.level + 1),
                         Cell(self.pos[0] + self.sizes[0]/4.0, self.pos[1] + self.sizes[1]/4.0, self.sizes/2.0, level = self.level + 1)]
            for d, daugther in enumerate(daugthers):
                alpha = 1.0 - (self.level) / (max_level+2.0)
                R = daugther.level / (max_level+2)
                #print 'd', d, '  level =', self.level, '  alpha = ', alpha, '  red =', R
                if (particles.is_inside(daugther)):
                    ax.plot([daugther.mins[0], daugther.maxs[0]], daugther.pos[1]*np.ones(2), '-', color = (R, 0, 0), linewidth = linewidth, alpha = alpha)
                    ax.plot(daugther.pos[0]*np.ones(2), [daugther.mins[1], daugther.maxs[1]], '-', color = (R, 0, 0), linewidth = linewidth, alpha = alpha)
                if (self.level < max_level and particles.is_inside(daugther)):
                    daugther.draw_daugthers(ax, particles)


x = r * np.cos(theta)
y = r * np.sin(theta)

particles = Part(x, y)
cell = Cell(0.0, 0.0, (2.05*abs(x).max(), 2.05*abs(y).max()))

fig = plot.figure()
ax = fig.add_subplot(1,1,1, xticks=[], yticks=[])
#ax = fig.add_subplot(1,1,1)

ax.plot(x[0:N/2], y[0:N/2], 'ob', ms = 14)
ax.plot(x[N/2:], y[N/2:],   'om', ms = 7)
cell.draw_quadtree(ax, particles)
ax.set_ylim(cell.mins[0], cell.maxs[0])
ax.set_ylim(cell.mins[1], cell.maxs[1])

plot.savefig('quadtree.svg')
plot.savefig('quadtree.pdf')
plot.show()
