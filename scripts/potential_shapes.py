#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import scipy.special
import plot

import plot_params
plot.matplotlib.rcParams['figure.subplot.hspace'] = 0.0

class Shape():
    def __init__(self, Z, r):
        self.Z = Z
        self.r = r
        # Long range Coulomb
        self.phi = Z / r
        self.E   = -Z / r**2
    def Plot(self, ax, color, label):
        ax[0].plot(self.r, self.phi, color, label=label, linewidth=3)
        ax[1].plot(self.r, self.E,   color, label=label, linewidth=3)

class Coulomb(Shape):
    def __init__(self, Z, r):
        Shape.__init__(self, Z, r)

class Harmonic(Shape):
    def __init__(self, Z, r, phi0):
        Shape.__init__(self, Z, r)

        self.phi0 = phi0
        self.R = 3.0 / (2.0 * self.phi0)
        #self.A = 1.0 / (2.0 * self.R**3)
        self.A = 4.0 * self.phi0**3 / 27.0
        self.ind = np.nonzero(self.r < self.R)

        self.phi[self.ind] = -self.A*self.r[self.ind]**2 + self.phi0
        self.E[self.ind]   = -2.0 * self.A * self.r[self.ind]


class SuperGaussian(Shape):
    def __init__(self, Z, r, phi0, m):
        Shape.__init__(self, Z, r)

        self.phi0 = phi0
        self.m    = m
        self.sigma = (self.Z*self.m**(1.0/(2.0*self.m)) / self.phi0) * np.exp(1.0 / (2.0 * self.m))
        self.R = (self.Z / self.phi0) * np.exp(1.0 / (2.0*self.m))
        self.ind = np.nonzero(self.r < self.R)

        self.phi[self.ind] = self.phi0 * np.exp(-0.5 * (self.r[self.ind]/self.sigma)**(2*self.m))
        self.E[self.ind]   = -(self.phi0 * self.m * np.exp(-0.5 * (self.r[self.ind]/self.sigma)**(2*self.m)) / self.r[self.ind]) * (self.r[self.ind]/self.sigma)**(2*m)


class ChargeDistribution(Shape):
    def __init__(self, Z, r, phi0):
        self.Z      = Z
        self.r      = r
        self.phi0   = phi0
        self.sigma  = np.sqrt(2.0 / np.pi) / self.phi0

        erf_r       = scipy.special.erf(self.r / (self.sigma * np.sqrt(2.0)))
        self.phi    = erf_r / self.r
        self.E      = -erf_r / (self.r**2) + np.sqrt(2.0/np.pi) * np.exp(-self.r**2 / (2.0*self.sigma)**2)/(self.sigma*self.r)


r = np.linspace(0.000001, 3.0, 1000.0)

Z    = 1
phi0 = 1.5 # Hartree

coulombic   = Coulomb(      Z, r)
harmonic    = Harmonic(     Z, r, phi0)
sg_one      = SuperGaussian(Z, r, phi0, 1)
sg_three    = SuperGaussian(Z, r, phi0, 3)
cs          = ChargeDistribution(Z, r, phi0)

fig = plot.figure(figsize = (12.0, 7.0))

ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2, sharex=ax1)
axes = [ax1, ax2]

coulombic.Plot( axes, ':k', 'Coulomb')
harmonic.Plot(  axes, '-b', 'Harmonic')
sg_one.Plot(    axes, '-r', 'Super-Gaussian (m=1)')
sg_three.Plot(  axes, '-c', 'Super-Gaussian (m=3)')
cs.Plot(        axes, '-m', 'Charge distribution')


ax1.grid()
ax1.set_ylabel('Potential [Eh]')
ax2.grid()
ax2.set_ylabel('Field [a.u.]')
ax2.set_xlabel('r [Bohr]')

ax1.set_ylim((0.0, 1.25*phi0))
ax2.set_ylim((-2.0, 0.05))
plot.setp(ax1.get_xticklabels(), visible=False)


leg = ax1.legend(loc = "best")
leg.get_frame().set_alpha(0.6)
#leg.set_zorder(100)
ax2.set_zorder(-100)

plot.savefig('potential_shapes.pdf')
plot.savefig('potential_shapes.svg')
plot.show()


