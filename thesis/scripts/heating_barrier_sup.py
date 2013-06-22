#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np

import plot
import plot_params

from photon import *
import fleches
import colors_and_symbols
import constants as cst


plot.matplotlib.rcParams['figure.subplot.hspace'] = 0.0
plot.matplotlib.rcParams['figure.subplot.wspace'] = 0.0


Xe_Z0_Ip = 12.265625 * cst.eV_to_Eh
# See ionization.git/src/data/ImpactCrossSections_Xe_GroundToExcited.cpp
# Eigenvalues
#Xe_Z0_es_E = np.asarray([-0.1235, -0.08198, -0.06377, -0.02002, -0.053875, -0.04072, -0.035485, -0.013905]) # [Eh]
#Xe_Z0_es_n = ['6s', '6p', '5d', '5f', '7s', '7p', '6d', '6f']
Xe_Z0_es_E = np.asarray([-0.1235, -0.08198, -0.06377, -0.02002]) # [Eh]
Xe_Z0_es_n = ['6s', '6p', '5d', '5f']
Xe_Z0_gses = Xe_Z0_Ip + Xe_Z0_es_E

Umax = 0.1
Umin = -2.0 * Xe_Z0_Ip
dU = Umax - Umin

print "Xe_Z0_gses [Eh] =", Xe_Z0_gses
print "Xe_Z0_gses [eV] =", Xe_Z0_gses * cst.Eh_to_eV


r = np.linspace(-20.0, 20.0, 1000)
lr = r[-1] - r[0]
r01 = 0.0
r02 = -8.2
r03 = 0.0
r04 = +8.2
cs = 1.0
U1 = -cs/abs(r-r01)
U2 = -cs/abs(r-r02) - cs/abs(r-r03) - cs/abs(r-r04)
U1[U1<Umin] = np.NaN
U2[U2<Umin] = np.NaN
r0 = abs(cs/Xe_Z0_Ip)

fig = plot.figure(figsize = (10.0, 7.0))

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

# ******************************************************
ax1.plot(r, U1)
ax1.plot([-r0, r0], [-Xe_Z0_Ip, -Xe_Z0_Ip], '-m', label = '5p')

# Electron
ax1.plot([0.0], [-Xe_Z0_Ip], 'og', ms = 14)

# Photon
gamma = 0.5*Xe_Z0_Ip
ar_e1 = fleches.arrow('e1', [0.0, -Xe_Z0_Ip], [0.0, gamma-Xe_Z0_Ip])
ar_e1.Plot(ax1, color = 'b', alpha = 0.5)
photon(ax1, [r[0]+lr*0.1, -Xe_Z0_Ip], [0.0, -Xe_Z0_Ip], Xe_Z0_Ip/5.0, alpha = 0.5)
ax1.plot([0.0], [gamma/2.0 - Xe_Z0_Ip], 'xr', ms = 30, markeredgewidth = 3)



# ******************************************************
Uep = -1.5*Xe_Z0_Ip
ax2.plot(r, U2)
ax2.plot([-r0, r0], [Uep, Uep], '-m', label = '5p')

# Electron
ax2.plot([0.0], [Uep], 'og', ms = 14, alpha = 0.6)
ax2.plot([0.0], [Uep+gamma], 'og', ms = 14  )

# Photon
gamma = 0.5*Xe_Z0_Ip
ar_e1 = fleches.arrow('e1', [0.0, Uep], [0.0, Uep+gamma])
ar_e1.Plot(ax2, color = 'b', alpha = 0.5)
photon(ax2, [r[0]+lr*0.1, Uep], [0.0, Uep], Xe_Z0_Ip/5.0)


for ax in [ax1, ax2]:
    ax.set_xlabel('r [bohr]')
    ax.set_ylabel('Energy [Hartree]')

    #leg = ax.legend(loc = "lower right")
    ax.axhline(0.0, linestyle = ':', color = 'k')

ax1.set_ylim((0.98*Umin,1.02*Umax))

# Remove overlapping labels
ax1.get_xticklabels()[-1].set_visible(False)
ax2.get_xticklabels()[0].set_visible(False)

for ext in ['pdf', 'svg']:
   plot.savefig('heating_barrier_sup.' + ext)
plot.show()
