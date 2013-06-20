#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import scipy.special
import plot

import plot_params
import constants as cst
import fleches


Xe_Z0_Ip = 12.265625 * cst.eV_to_Eh

#Umin = -1.2 * Xe_Z0_Ip
Umin = -2.2 * Xe_Z0_Ip
Umax = 0.25

r = np.linspace(-20.0, 20.0, 1000)
cs = 1.0
U = -cs / abs(r)
U[U<Umin] = np.NaN

laser = -r / 30.0
Ubent = U + laser

lr = r[-1] - r[0]
lU = Umax - Umin

r0 = abs(cs/Xe_Z0_Ip)

fig = plot.figure()

ax1, ax1b = plot.get_axis_two_scales(fig,
                                     scale_y = cst.Eh_to_eV,
                                     ax2_ylabel = 'Energy [eV]')

# Ip
ax1.plot([-r0, r0], [-Xe_Z0_Ip, -Xe_Z0_Ip], '-m')

# Threshold
ax1.axhline(0.0, color = 'k', ls = '-', alpha = 0.5)

# U(r)
ax1.plot(r, U, '--k', alpha = 0.5, label = 'Unperturbed ion')

# Laser
ax1.plot(r, laser, '-r', label = 'Laser')

# U(r) + laser
ax1.plot(r, Ubent, '-k', label = 'Effective')

# Electron
ax1.plot([0.0], [-Xe_Z0_Ip], 'ob', ms = 14)
ar1 = fleches.arrow('Tunnel', [0.0, -Xe_Z0_Ip], [0.9*r[-1], -Xe_Z0_Ip])
ar1.Plot(ax1, color = 'g')

ax1.set_xlim((r[0], r[-1]))
ax1.set_ylim((0.95*Umin, Umax))

ax1.set_xlabel('r [bohr]')
ax1.set_ylabel('Energy [Hartree]')

leg = ax1.legend(loc = "best")
leg.get_frame().set_alpha(0.75)

for ext in ['pdf', 'svg']:
    plot.savefig('ionization_tunnel.' + ext)
plot.show()

