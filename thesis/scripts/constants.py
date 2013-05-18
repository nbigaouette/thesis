#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy

pi              = numpy.pi
twoPi           = 2.0 * pi
one_over_four_Pi= 1.0 / (4.0 * pi)
c               = 299792458.0           # Speed of ligth in vacuum [m/s]
planck          = 6.62606896e-34        # Planck constant [J . s]
me              = 9.10938188e-31        # Electron mass [kg]
e0              = 1.602176462e-19       # Unitary charge [C]
inv_c2          = 1.0 / (c * c)         # 1/c^2 [s^2/m^2]
mu0             = 4.0e-7 * pi           # Magnetic constant (vacuum) [N/A^2]
eps0            = inv_c2 / mu0          # Electric constant (vacuum) [F/m]

one_over_4Pieps0= one_over_four_Pi / eps0 # Coulomb's law constant [N.m^2/C^2] [C]

hbar            = planck / twoPi        # Planck cst / 2 pi [J . s]
a0              = hbar*hbar / (one_over_4Pieps0 * me * e0*e0)
Eh              = hbar*hbar/(me*a0*a0)  # Hartree [J]

nm_to_m         = 1.0e-9
m_to_nm         = 1.0 / nm_to_m
angstrom_to_m   = 1.0e-10
m_to_angstrom   = 1.0 / angstrom_to_m
mum_to_m        = 1.0e-6
m_to_mum        = 1.0 / mum_to_m

Eh_to_eV        = Eh / e0               # eV . Eh^-1
eV_to_Eh        = 1.0 / Eh_to_eV        # Eh . eV^-1
eV_to_J         = e0                    # J . eV^-1
J_to_eV         = 1.0 / eV_to_J         # eV . J^-1
J_to_Eh         = eV_to_Eh * J_to_eV    # Eh . J^-1
Eh_to_J         = 1.0 / J_to_Eh         # J . Eh^-1
