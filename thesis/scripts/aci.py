#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#import matplotlib.pyplot as plt
#import matplotlib.transforms as mtransforms
#from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
import numpy as np


import plot
import plot_params

import colors_and_symbols
import constants as cst


plot.matplotlib.rcParams['figure.subplot.hspace'] = 0.0
plot.matplotlib.rcParams['figure.subplot.wspace'] = 0.0

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
             bidirectional = False):
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
                                  annotation_clip=False )
            ann.arrow_patch.set_clip_box(ax.bbox)
            #ax.plot([self.start_xy[0]], [self.start_xy[1]], 'o' + color)
        if (label != None):
            dxy = self.stop_xy + self.start_xy
            ax.text(dxy[0]/2.0, dxy[1]/2.0,
                    label + label_equal_zero,
                    horizontalalignment = horizontalalignment,
                    verticalalignment = verticalalignment,
                    color = color)


Xe_Z0_Ip = 12.265625 * cst.eV_to_Eh
# See ionization.git/src/data/ImpactCrossSections_Xe_GroundToExcited.cpp
# Eigenvalues
#Xe_Z0_es_E = np.asarray([-0.1235, -0.08198, -0.06377, -0.02002, -0.053875, -0.04072, -0.035485, -0.013905]) # [Eh]
#Xe_Z0_es_n = ['6s', '6p', '5d', '5f', '7s', '7p', '6d', '6f']
Xe_Z0_es_E = np.asarray([-0.1235, -0.08198, -0.06377, -0.02002]) # [Eh]
Xe_Z0_es_n = ['6s', '6p', '5d', '5f']
Xe_Z0_gses = Xe_Z0_Ip + Xe_Z0_es_E

#Umin = -1.2 * Xe_Z0_Ip
Umin = -1.8 * Xe_Z0_Ip

print "Xe_Z0_gses [Eh] =", Xe_Z0_gses
print "Xe_Z0_gses [eV] =", Xe_Z0_gses * cst.Eh_to_eV

Ke = 11.0 * cst.eV_to_Eh
Ke2 = 5.0 * cst.eV_to_Eh

r = np.linspace(-20.0, 20.0, 1000)
cs = 1.0
U = -cs / abs(r)
U[U<Umin] = np.NaN

#fig1 = plot.figure()
#fig2 = plot.figure()
#fig3 = plot.figure()
#fig4 = plot.figure()

#ax1, ax1b = plot.get_axis_two_scales(fig1,
                                    #scale_y = cst.Eh_to_eV,
                                    #ax2_ylabel = 'Energy [eV]')
#ax2, ax2b = plot.get_axis_two_scales(fig2,
                                    #scale_y = cst.Eh_to_eV,
                                    #ax2_ylabel = 'Energy [eV]')
#ax3, ax3b = plot.get_axis_two_scales(fig3,
                                    #scale_y = cst.Eh_to_eV,
                                    #ax2_ylabel = 'Energy [eV]')
#ax4, ax4b = plot.get_axis_two_scales(fig4,
                                    #scale_y = cst.Eh_to_eV,
                                    #ax2_ylabel = 'Energy [eV]')

fig = plot.figure()
ax1, ax1b = plot.get_axis_two_scales(fig,
                                     scale_y = cst.Eh_to_eV,
                                     ax2_ylabel = 'Energy [eV]',
                                     subplot = 221)
ax2, ax2b = plot.get_axis_two_scales(fig,
                                     scale_y = cst.Eh_to_eV,
                                     ax2_ylabel = 'Energy [eV]',
                                     sharex = ax1, sharey = ax1,
                                     subplot = 222)
ax3, ax3b = plot.get_axis_two_scales(fig,
                                     scale_y = cst.Eh_to_eV,
                                     ax2_ylabel = 'Energy [eV]',
                                     sharex = ax1, sharey = ax1,
                                     subplot = 223)
ax4, ax4b = plot.get_axis_two_scales(fig,
                                     scale_y = cst.Eh_to_eV,
                                     ax2_ylabel = 'Energy [eV]',
                                     sharex = ax1, sharey = ax1,
                                     subplot = 224)
ax1b.get_yaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)
ax3b.get_yaxis().set_visible(False)
ax4.get_yaxis().set_visible(False)
for ax in [ax1, ax1b, ax2, ax2b]:
    ax.get_xaxis().set_visible(False)

# ********************** Before ************************
ax1.plot(r, U)
ax1.axhline(-Xe_Z0_Ip, linestyle = '-', color = 'm', label = '5p')
for ei in xrange(len(Xe_Z0_es_E)):
    ax1.axhline(Xe_Z0_es_E[ei], linestyle = colors_and_symbols.symbol(ei+5),
                                color = colors_and_symbols.color(ei),
                                label = Xe_Z0_es_n[ei])

# Impacting electron
ar1Ke = arrow('K_e', [0.75*r[0], 0.0], [0.75*r[0], Ke])    # Vertical (Ke)
ar2Ke = arrow('K_e', [0.75*r[0], Ke],  [0.75*r[0]+10, Ke]) # Horizontal (v)
ar1Ke.Plot(ax1, color = 'b', label = '$K_e$', horizontalalignment = 'right', bidirectional = True)
ar2Ke.Plot(ax1, color = 'b', label = '$v_e$', verticalalignment = 'top')
ax1.plot([0.75*r[0]], [Ke], 'ob', ms = 14)

# Bound electron
ax1.plot([0.0], [-Xe_Z0_Ip], 'og', ms = 14)


# ********************** During ************************
ax2.plot(r, U)
ax2.axhline(-Xe_Z0_Ip, linestyle = '-', color = 'm', label = '5p')
for ei in xrange(len(Xe_Z0_es_E)):
    ax2.axhline(Xe_Z0_es_E[ei], linestyle = colors_and_symbols.symbol(ei+5),
                                color = colors_and_symbols.color(ei),
                                label = Xe_Z0_es_n[ei])
# Impacting electron
Kep = Ke - Xe_Z0_gses[0]
ar1Ke = arrow("K_e'", [0.0, 0.0], [0.0, Kep]) # Vertical (Ke)
ar2Ke = arrow("K_e'", [0.0, Kep], [5.0, Kep]) # Horizontal (v)
ar1Ke.Plot(ax2, color = 'b', label = r"$K_e'$   .", horizontalalignment = 'right', bidirectional = True)
ar2Ke.Plot(ax2, color = 'b', label = r"$v_e'$", verticalalignment = 'bottom')
ax2.plot([0.0], [Kep], 'ob', ms = 14)

# Bound electron
ax2.plot([0.0], [Xe_Z0_es_E[0]], 'og', ms = 14)
ax2.plot([0.0], [-Xe_Z0_Ip], 'og', ms = 14, alpha = 0.4)
arBe = arrow("", [0.0, -Xe_Z0_Ip], [0.0, Xe_Z0_es_E[0]]) # Vertical
arBe.Plot(ax2, color = 'g')


# ********************** During 2 ************************
ax3.plot(r, U)
ax3.axhline(-Xe_Z0_Ip, linestyle = '-', color = 'm', label = '5p')
for ei in xrange(len(Xe_Z0_es_E)):
    ax3.axhline(Xe_Z0_es_E[ei], linestyle = colors_and_symbols.symbol(ei+5),
                                color = colors_and_symbols.color(ei),
                                label = Xe_Z0_es_n[ei])
# Impacting electron
Kep = Ke - Xe_Z0_gses[0]
ar1Ke = arrow("K_e'", [0.75*r[-1], 0.0], [0.75*r[-1], Kep]) # Vertical (Ke)
ar2Ke = arrow("K_e'", [0.75*r[-1], Kep], [0.75*r[-1]+5.0, Kep]) # Horizontal (v)
ar1Ke.Plot(ax3, color = 'b', label = r"$K_e'$   .", horizontalalignment = 'right', bidirectional = True)
ar2Ke.Plot(ax3, color = 'b', label = r"$v_e'$", verticalalignment = 'bottom')
ax3.plot([0.75*r[-1]], [Kep], 'ob', ms = 14)

# New Impacting electron
ar1Ke = arrow('K_{e,2}', [0.75*r[0], 0.0], [0.75*r[0], Ke2])    # Vertical (Ke)
ar2Ke = arrow('K_{e,2}', [0.75*r[0], Ke2],  [0.75*r[0]+5, Ke2]) # Horizontal (v)
ar1Ke.Plot(ax3, color = 'r', label = '$K_{e,2}$', horizontalalignment = 'right', bidirectional = True)
ar2Ke.Plot(ax3, color = 'r', label = '$v_{e,2}$', verticalalignment = 'top')
ax3.plot([0.75*r[0]], [Ke2], 'or', ms = 14)

# Bound electron
ax3.plot([0.0], [Xe_Z0_es_E[0]], 'og', ms = 14)


# ********************** End/After ************************
ax4.plot(r, U)
ax4.axhline(-Xe_Z0_Ip, linestyle = '-', color = 'm', label = '5p')
for ei in xrange(len(Xe_Z0_es_E)):
    ax4.axhline(Xe_Z0_es_E[ei], linestyle = colors_and_symbols.symbol(ei+5),
                                color = colors_and_symbols.color(ei),
                                label = Xe_Z0_es_n[ei])

# New Impacting electron
#ar1Ke = arrow('K_{e,2}', [0.0, 0.0], [0.0, Ke2])    # Vertical (Ke)
ar2Ke = arrow('K_{e,2}', [0.0, Ke2],  [0.0+2.5, Ke2]) # Horizontal (v)
#ar1Ke.Plot(ax4, color = 'r', label = "$K'_{e,2}$", horizontalalignment = 'right', bidirectional = True)
ar2Ke.Plot(ax4, color = 'r', label = "$v'_{e,2}$", verticalalignment = 'top')
ax4.plot([0.0], [Ke2], 'or', ms = 14)

## Bound electron, now free
ax4.plot([0.0], [Xe_Z0_es_E[0]], 'og', ms = 14, alpha = 0.4)
ax4.plot([0.0], [Ke2+Xe_Z0_es_E[0]], 'og', ms = 14)
ar1Ke = arrow('', [0.0, Xe_Z0_es_E[0]], [0.0, Ke2+Xe_Z0_es_E[0]])    # Vertical (Ke)
ar1Ke.Plot(ax4, color = 'g')
ar2Ke = arrow('', [0.0, Ke2+Xe_Z0_es_E[0]],  [0.0+2.5, Ke2+Xe_Z0_es_E[0]]) # Horizontal (v)
ar2Ke.Plot(ax4, color = 'g', label = "$v_{e,3}$", verticalalignment = 'top')

for ax in [ax1, ax2, ax3, ax4]:
    ax.set_xlabel('r [bohr]')
    ax.set_ylabel('Energy [Hartree]')

    #ax.set_ylim((Umin, 0.05*abs(Umin)))
    ax.set_ylim((0.97*Umin, 1.05*Ke))
    #leg = ax.legend(loc = "lower right")
    ax.axhline(0.0, linestyle = ':', color = 'k')


ax1.text(0.9*r[0], 0.9*Umin, 'a)')
ax2.text(0.9*r[0], 0.9*Umin, 'b)')
ax3.text(0.9*r[0], 0.9*Umin, 'c)')
ax4.text(0.9*r[0], 0.9*Umin, 'd)')

#for ext in ['pdf', 'svg']:
   #plot.savefig(['aci_0.' + ext, 'aci_1.' +  ext, 'aci_2.' +  ext, 'aci_3.' +  ext])
for ext in ['pdf', 'svg']:
   plot.savefig('ionization_aci.' + ext)
plot.show()
