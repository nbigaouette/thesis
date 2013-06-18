#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import scipy.special
import plot

import plot_params
import constants as cst
import fleches

class photon:
    def __init__(self, ax, xy0, xy1, A):
        self.xy0 = np.asarray(xy0)
        self.xy1 = np.asarray(xy1)
        self.A  = A
        self.dx = self.xy1[0] - self.xy0[0]
        self.dy = self.xy1[1] - self.xy0[1]
        self.x = np.linspace(self.xy0[0], self.xy1[0], 1000)
        self.x0 = self.xy0[0] + self.dx / 2.0
        self.fwhm  = 0.25 * self.dx
        self.N = 10.0
        self.lamda = self.dx / self.N
        self.k = 2.0*np.pi / self.lamda
        self.sigma = self.fwhm / np.sqrt(2.0 * np.log(2.0))
        #print "A     =", self.A
        #print "xy0   =", self.xy0
        #print "xy1   =", self.xy1
        #print "dx    =", self.dx
        #print "dy    =", self.dy
        #print "x0    =", self.x0
        #print "fwhm  =", self.fwhm
        #print "N     =", self.N
        #print "lamda =", self.lamda
        #print "k     =", self.k
        #print "sigma =", self.sigma
        self.y = self.xy0[1] + self.A * np.exp(-(self.x - self.x0)**2 / self.sigma**2) * np.sin(self.k * (self.x - self.x0))
        ax.plot(self.x, self.y, '-r')
        ann = ax.annotate('', xy = self.xy1,
                          xycoords='data',
                          xytext = (self.xy1[0]-self.dx/10000.0, self.xy1[1]),
                          arrowprops=dict(arrowstyle = "-|>", color = 'r'),
                          annotation_clip=False)
        ann.arrow_patch.set_clip_box(ax.bbox)

Xe_Z0_Ip = 12.265625 * cst.eV_to_Eh

Umin = -1.2 * Xe_Z0_Ip
Umax = 0.25

GammaE = 1.25 * Xe_Z0_Ip

r = np.linspace(-20.0, 20.0, 1000)
cs = 1.0
U = -cs / abs(r)
U[U<Umin] = np.NaN

lr = r[-1] - r[0]
lU = Umax - Umin

r0 = abs(cs/Xe_Z0_Ip)

fig = plot.figure()

ax1 = fig.add_subplot(1,1,1)

# Ip
ax1.plot([-r0, r0], [-Xe_Z0_Ip, -Xe_Z0_Ip], '-m')

# U(r)
ax1.plot(r, U, '-k')

# Continuum
ax1.fill([r[0],r[-1],r[-1],r[0]],[0.0,0.0,Umax,Umax], fill=False, hatch='/', alpha = 0.5)

# Photon energy
ar1 = fleches.arrow('GammaE', [0.0, -Xe_Z0_Ip], [0.0, GammaE-Xe_Z0_Ip])
ar1.Plot(ax1, color = 'r', label = r'$E_\gamma$', horizontalalignment = 'left')
ax1.plot([0.0], [-Xe_Z0_Ip], 'ob', ms = 14, alpha = 0.6)

# New electron velocity
ar2 = fleches.arrow('N_e', [0.0, GammaE-Xe_Z0_Ip], [5.0, GammaE-Xe_Z0_Ip])
ar2.Plot(ax1, color = 'b', label = '$v_e$', horizontalalignment = 'center', verticalalignment = 'bottom')
ax1.plot([0.0], [GammaE-Xe_Z0_Ip], 'ob', ms = 14)

# Delta E
x = -1.0
ar3 = fleches.arrow('D_e', [x, 0.0], [x, GammaE-Xe_Z0_Ip])
ar3.Plot(ax1, color = 'g', label = '$\Delta E.$', horizontalalignment = 'right', bidirectional = True)

# Photon
#photon(ax1, [r[0]+lr*0.1, -Xe_Z0_Ip/2.0], [0.0, -Xe_Z0_Ip/2.0], Xe_Z0_Ip/5.0)
photon(ax1, [r[0]+lr*0.1, -Xe_Z0_Ip], [0.0, -Xe_Z0_Ip], Xe_Z0_Ip/5.0)

ax1.set_xlim((r[0], r[-1]))
ax1.set_ylim((Umin, Umax))

ax1.set_xlabel('r [bohr]')
ax1.set_ylabel('Energy [Hartree]')

for ext in ['pdf', 'svg']:
    plot.savefig('ionization_single.' + ext)
plot.show()

