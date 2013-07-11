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


plot.matplotlib.rcParams['figure.subplot.hspace'] = 0.0
plot.matplotlib.rcParams['figure.subplot.wspace'] = 0.0



Xe_Z0_Ip = 12.265625 * cst.eV_to_Eh

#Umin = -1.2 * Xe_Z0_Ip
Umin = -1.8 * Xe_Z0_Ip


Ke = 25.0 * cst.eV_to_Eh
Kep = (Ke - Xe_Z0_Ip) / 2.0

print "Ip  =", Xe_Z0_Ip, "Eh =", Xe_Z0_Ip*cst.Eh_to_eV, 'eV'
print "Ke  =", Ke, "Eh =", Ke*cst.Eh_to_eV, 'eV'
print "Kep =", Kep, "Eh =", Kep*cst.Eh_to_eV, 'eV'

Umax = 1.1*Ke

r = np.linspace(-20.0, 20.0, 1000)
cs = 1.0
U  = -cs / abs(r)
U[U<Umin] = np.NaN
r0 = abs(cs/Xe_Z0_Ip)

fig = plot.figure()
ax1, ax1b = plot.get_axis_two_scales(fig,
                                     scale_y = cst.Eh_to_eV,
                                     ax2_ylabel = 'Energy [eV]',
                                     subplot = 121)
ax2, ax2b = plot.get_axis_two_scales(fig,
                                     scale_y = cst.Eh_to_eV,
                                     ax2_ylabel = 'Energy [eV]',
                                     sharex = ax1, sharey = ax1,
                                     subplot = 122)
ax1b.get_yaxis().set_visible(False)
ax2.get_yaxis().set_visible(False)

# ********************** Before ************************

# Continuum
ax1.fill([r[0],r[-1],r[-1],r[0]],[0.0,0.0,Umax,Umax], fill=False, hatch='/', alpha = 0.5)

# Ip
#ax1.axhline(-Xe_Z0_Ip, linestyle = '-', color = 'm', label = '5p')
ax1.plot([-r0, r0], [-Xe_Z0_Ip, -Xe_Z0_Ip], '-m')

ax1.plot(r, U, '-k')

# Impacting electron
ar1Ke = fleches.arrow('K_e', [0.75*r[0], 0.0], [0.75*r[0], Ke])    # Vertical (Ke)
ar2Ke = fleches.arrow('K_e', [0.75*r[0], Ke],  [0.75*r[0]+10, Ke]) # Horizontal (v)
ar1Ke.Plot(ax1, color = 'b', label = '$K_e$', horizontalalignment = 'right', bidirectional = True)
ar2Ke.Plot(ax1, color = 'b', label = '$v_e$', verticalalignment = 'top')
ax1.plot([0.75*r[0]], [Ke], 'ob', ms = 14)

# Bound electron
ax1.plot([0.0], [-Xe_Z0_Ip], 'og', ms = 14)


# ********************** After ************************

# Continuum
ax2.fill([r[0],r[-1],r[-1],r[0]],[0.0,0.0,Umax,Umax], fill=False, hatch='/', alpha = 0.5)

# Ip
#ax2.axhline(-Xe_Z0_Ip, linestyle = '-', color = 'm', label = '5p')
ax2.plot([-r0, r0], [-Xe_Z0_Ip, -Xe_Z0_Ip], '-m')

ax2.plot(r, U, '-k')

# New electron
ar1Ke = fleches.arrow("K_e'", [0.0, 0.0], [0.0, Kep]) # Vertical (Ke)
ar2Ke = fleches.arrow("K_e'", [0.0, Kep*1.3], [5.0, Kep*1.3]) # Horizontal (v)
ar1Ke.Plot(ax2, color = 'b', label = r"$K_e'$   .", horizontalalignment = 'right', bidirectional = True)
ar2Ke.Plot(ax2, color = 'b', verticalalignment = 'bottom')
ax2.plot([0.0], [Kep*1.3], 'ob', ms = 14)

# Bound electron
ax2.plot([0.0], [-Xe_Z0_Ip], 'og', ms = 14, alpha = 0.4)
ax2.plot([0.0], [Kep], 'og', ms = 14)
ar2a = fleches.arrow("", [0.0, 0.0], [0.0, Kep]) # Vertical
ar2a.Plot(ax2, color = 'g', label = r".$K_e'$", horizontalalignment = 'left', bidirectional = True)
ar2b = fleches.arrow("K_e'", [0.0, Kep], [5.0, Kep]) # Horizontal (v)
ar2b.Plot(ax2, color = 'g', verticalalignment = 'bottom')

for ax in [ax1, ax2]:
    ax.set_xlabel('r [bohr]')
    ax.set_ylabel('Energy [Hartree]')

    #ax.set_ylim((Umin, 0.05*abs(Umin)))
    ax.set_ylim((0.97*Umin, 1.05*Ke))
    #leg = ax.legend(loc = "lower right")
    ax.axhline(0.0, linestyle = ':', color = 'k')


ax1.text(0.9*r[0], 0.9*Umin, 'a)')
ax2.text(0.9*r[0], 0.9*Umin, 'b)')

for ext in ['pdf', 'svg']:
   plot.savefig('ionization_impact.' + ext)
plot.show()
