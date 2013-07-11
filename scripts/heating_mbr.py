#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#import matplotlib.pyplot as plt
#import matplotlib.transforms as mtransforms
#from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
import numpy as np


import plot
import plot_params

import fleches
import colors_and_symbols
import constants as cst


plot.matplotlib.rcParams['patch.linewidth'] = 3.0


Xe_Z0_Ip = 12.265625 * cst.eV_to_Eh
# See ionization.git/src/data/ImpactCrossSections_Xe_GroundToExcited.cpp
# Eigenvalues
#Xe_Z0_es_E = np.asarray([-0.1235, -0.08198, -0.06377, -0.02002, -0.053875, -0.04072, -0.035485, -0.013905]) # [Eh]
#Xe_Z0_es_n = ['6s', '6p', '5d', '5f', '7s', '7p', '6d', '6f']
Xe_Z0_es_E = np.asarray([-0.1235, -0.08198, -0.06377, -0.02002]) # [Eh]
Xe_Z0_es_n = ['6s', '6p', '5d', '5f']
Xe_Z0_gses = Xe_Z0_Ip + Xe_Z0_es_E

Umax = 0.5
Umin = -0.8 * Xe_Z0_Ip
dU = Umax - Umin

print "Xe_Z0_gses [Eh] =", Xe_Z0_gses
print "Xe_Z0_gses [eV] =", Xe_Z0_gses * cst.Eh_to_eV


r = np.linspace(-20.0, 20.0, 1000)
cs = 1.0
U = -cs / abs(r)
U[U<Umin] = np.NaN
r0 = abs(cs/Xe_Z0_Ip)


fig = plot.figure(figsize = (10.0, 7.0))

ax1, ax1b = plot.get_axis_two_scales(fig,
                                     scale_y = cst.Eh_to_eV,
                                     ax2_ylabel = 'Energy [eV]')

ax1.plot(r, U)
for ei in xrange(len(Xe_Z0_es_E)):
    r0e = min(abs(cs/Xe_Z0_es_E[ei]), r[-1])
    ax1.plot([-r0e, r0e], [Xe_Z0_es_E[ei], Xe_Z0_es_E[ei]],
                linestyle = colors_and_symbols.symbol(ei+5),
                color = colors_and_symbols.color(ei),
                label = Xe_Z0_es_n[ei])

# 3 electrons
al = 0.09 # arrow length
e1_xy = [0.25*r[0],  0.8*Umax]
e2_xy = [0.20*r[-1], 0.25*Umax]
e3_xy = [0.75*r[0],  0.1*Umax]
ar_e1 = fleches.arrow('e1', e1_xy, [e1_xy[0], e1_xy[1]    +al*dU])
ar_e2 = fleches.arrow('e2', e2_xy, [e2_xy[0], e2_xy[1]-2.0*al*dU])
ar_e3 = fleches.arrow('e3', e3_xy, [e3_xy[0], e3_xy[1]    +al*dU])
ar_e1.Plot(ax1, color = 'b')
ar_e2.Plot(ax1, color = 'b')
ar_e3.Plot(ax1, color = 'b')
ax1.plot([e1_xy[0]], [e1_xy[1]], 'ob', ms = 14)
ax1.plot([e2_xy[0]], [e2_xy[1]], 'ob', ms = 14)
ax1.plot([e3_xy[0]], [e3_xy[1]], 'ob', ms = 14)

ax1.set_xlabel('r [bohr]')
ax1.set_ylabel('Energy [Hartree]')

#leg = ax1.legend(loc = "lower right")
ax1.axhline(0.0, linestyle = ':', color = 'k')

ax1.set_ylim((0.98*Umin,1.02*Umax))

for ext in ['pdf', 'svg']:
   plot.savefig('heating_mbr.' + ext)
plot.show()
