#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np

class photon:
    def __init__(self, ax, xy0, xy1, A, alpha = 1.0):
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
        ax.plot(self.x, self.y, '-r', alpha = alpha)
        ann = ax.annotate('', xy = self.xy1,
                          xycoords='data',
                          xytext = (self.xy1[0]-self.dx/10000.0, self.xy1[1]),
                          arrowprops=dict(arrowstyle = "-|>", color = 'r'),
                          annotation_clip=False, alpha = alpha)
        ann.arrow_patch.set_clip_box(ax.bbox)
