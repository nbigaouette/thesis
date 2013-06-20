#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np

import constants as cst


class arrow:
    def __init__(self, name, start_xy, stop_xy):
        self.name       = name
        self.start_xy   = np.asarray(start_xy)
        self.stop_xy    = np.asarray(stop_xy)
        self.length     = np.sqrt(((self.stop_xy - self.start_xy)**2).sum())
    def PrintHeader(self):
        print "Arrow name           Start         Stop         Length            Start          Stop         Length"
    def Print(self):
        print str('%-16s' % self.name) + " " + \
                    str('%8.4f' % self.start_xy) + " Eh   " + \
                    str('%8.4f' % self.stop_xy) + " Eh   " + \
                    str('%8.4f' % self.length) + " Eh   " + \
                "    " + \
                    str('%8.4f' % (self.start_xy * cst.Eh_to_eV)) + " eV   " + \
                    str('%8.4f' % (self.stop_xy * cst.Eh_to_eV)) + " eV   " + \
                    str('%8.4f' % (self.length * cst.Eh_to_eV)) + " eV"
    def Is_Null(self):
        if (self.length < 1.0e-5):
            return True
        else:
            return False
    def Plot(self, ax, color = 'k', label = None,
             horizontalalignment = 'center', verticalalignment = 'center',
             bidirectional = False,
             alpha = 1.0):
        if (self.Is_Null()):
            label_equal_zero = r" = 0"
        else:
            label_equal_zero = r""
            if (bidirectional):
                arrowstyle = "<|-|>"
            else:
                arrowstyle = "-|>"
            ann = ax.annotate('', xy     = self.stop_xy,  xycoords='data',
                                  xytext = self.start_xy, textcoords='data',
                                  arrowprops=dict(arrowstyle = arrowstyle, color = color),
                                  annotation_clip=False, alpha = alpha)
            ann.arrow_patch.set_clip_box(ax.bbox)
            #ax.plot([self.start_xy[0]], [self.start_xy[1]], 'o' + color)
        if (label != None):
            dxy = self.stop_xy + self.start_xy
            ax.text(dxy[0]/2.0, dxy[1]/2.0,
                    label + label_equal_zero,
                    horizontalalignment = horizontalalignment,
                    verticalalignment = verticalalignment,
                    color = color)
