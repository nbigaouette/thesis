#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy

pi              = numpy.pi
twoPi           = 2.0 * pi
sqrt_Pi        = numpy.sqrt(pi)
sqrt_2         = numpy.sqrt(2.0)
sqrt_2_over_pi = numpy.sqrt(2.0 / pi)
one_over_Pi    = 1.0 / pi

two_over_three = 2.0 / 3.0
three_over_two = 3.0 / 2.0

one_over_three = 1.0 / 3.0
one_over_four  = 1.0 / 4.0
one_over_seven = 1.0 / 7.0

one_over_four_Pi= 1.0 / (4.0 * pi)
sqrt_one_over_four_Pi = numpy.sqrt(one_over_four_Pi)

four_over_twenty_seven = 4.0 / 27.0
four_over_three        = 4.0 / 3.0
eight_over_three       = 8.0 / 3.0
three_over_eight       = 3.0 / 8.0
three_over_sixteen     = 3.0 / 16.0



acos_pow_half_fourth   = numpy.arccos(0.5**0.25)
sin4_fwhm_to_period    = pi / acos_pow_half_fourth






# ***********************************************************
# ******* http://physics.nist.gov/cuu/Constants *************
# ***********************************************************


# ***********************************************************
# ******* SI

c               = 299792458.0           # Speed of ligth in vacuum [m/s]
me              = 9.10938188e-31        # Electron mass [kg]
e0              = 1.602176462e-19       # Unitary charge [C]
inv_c2          = 1.0 / (c * c)         # 1/c^2 [s^2/m^2]
mu0             = 4.0e-7 * pi           # Magnetic constant (vacuum) [N/A^2]
eps0            = inv_c2 / mu0          # Electric constant (vacuum) [F/m]

one_over_4Pieps0= one_over_four_Pi / eps0 # Coulomb's law constant [N.m^2/C^2] [C]

planck          = 6.62606896e-34        # Planck constant [J . s]
hbar            = planck / twoPi        # Planck cst / 2 pi [J . s]
a0              = hbar*hbar / (one_over_4Pieps0 * me * e0*e0)
                                        # Bohr radius [m]


Na             = 6.02214179e23        # Avogadro constant [mol^-1]
R              = 8.3144               # Gas constant [J . kg^-1 . mol^-1]
kB             = R / Na               # Boltzmann constant [J . K^-1]

# ***********************************************************
# ******* Atomic Units

Eh             = hbar*hbar/(me*a0*a0) # ...of energy (Hartree) [J]
au_energy      = Eh                   # ...of energy (Hartree) [J]
au_time        = hbar / Eh            # ...of time [s]
au_electric_field = Eh / (e0 * a0)    # ...of electric field [V . m^-1]
w_au           = 1.0 / au_time        # ...of frequency [???]
alpha          = 7.29735257e-3        # fine-structure constant

# ***********************************************************
# ******* Conversions

Eh_to_eV       = Eh / e0              # eV . Eh^-1
eV_to_Eh       = 1.0 / Eh_to_eV       # Eh . eV^-1
eV_to_J        = e0                   # J . eV^-1
J_to_eV        = 1.0 / eV_to_J        # eV . J^-1
J_to_Eh        = eV_to_Eh * J_to_eV   # Eh . J^-1
Eh_to_J        = 1.0 / J_to_Eh        # J . Eh^-1

m_to_angstrom  = 1.0e10               # Å . m^-1
angstrom_to_m  = 1.0 / m_to_angstrom  # m . Å^-1
bohr_to_m      = a0                   # m . bohr^-1
m_to_bohr      = 1.0 / bohr_to_m      # bohr . m^-1
m_to_cm        = 100.0                # cm . m^-1
cm_to_m        = 1.0 / m_to_cm        # cm . m^-1
m_to_nm        = 1.0e9                # nm . m^-1
nm_to_m        = 1.0 / m_to_nm        # m . nm^-1
bohr_to_nm     = bohr_to_m * m_to_nm  # nm . bohr^-1
nm_to_bohr     = 1.0 / bohr_to_nm     # bohr . nm^-1
bohr_to_angstrom  = bohr_to_m * m_to_angstrom   # Å . bohr^-1
angstrom_to_bohr  = 1.0 / bohr_to_angstrom      # bohr . Å^-1



fs_to_s        = 1.0e-15              # s . fs^-1
s_to_fs        = 1.0 / fs_to_s        # fs . s^-1
as_to_s        = 1.0e-18              # s . as^-1
s_to_as        = 1.0 / as_to_s        # as . s^-1
fs_to_as       = fs_to_s * s_to_as    # as . fs^-1
as_to_fs       = 1.0 / fs_to_as       # fs . as^-1

deg_to_rad     = 2.0 * pi / 360.0     # radians . degrees^-1
rad_to_deg     = 1.0 / deg_to_rad     # degrees . radians^-1


# ********** SI  <--> Atomic Units Conversions
au_to_si_mass  = me                   # kg . au^-1
si_to_au_mass  = 1.0 / au_to_si_mass  # au . kg^-1

au_to_si_length= a0                   # m . bohr^-1
si_to_au_length= 1.0 / au_to_si_length# bohr . m^-1

au_to_si_charge= e0                   # C . au^-1
si_to_au_charge= 1.0 / au_to_si_charge# au . C^-1

au_to_si_energy= Eh                   # J . au^-1
si_to_au_energy= 1.0 / au_to_si_energy# au . J^-1

au_to_si_k     = one_over_4Pieps0     # C^-2 . N . m^2 . au^-1
si_to_au_k     = 1.0 / au_to_si_energy# au . C^2 . N^-1 . m^-2

au_to_si_time  = hbar / Eh            # s . au^-1
si_to_au_time  = 1.0 / au_to_si_time  # au . s^-1

au_to_si_force = Eh / a0              # N . au^-1
si_to_au_force = 1.0 / au_to_si_force # au . N^-1

au_to_si_field = au_electric_field    # V . m^-1 . au^-1
si_to_au_field = 1.0 / au_to_si_field # au . m . V^-1

au_to_si_pot   = au_electric_field * a0 # V . au^-1
si_to_au_pot   = 1.0 / au_to_si_pot   # au . V^-1

au_to_si_vel   = Eh * a0 / hbar       # m.s^-1 . au^-1
si_to_au_vel   = 1.0 / au_to_si_vel   # au . m.s^-1

au_to_si_temp  = Eh / kB              # K . au^-1
si_to_au_temp  = 1.0 / au_to_si_temp  # au . K^-1


au_to_fs       = au_to_si_time * s_to_fs # fs . au^-1
fs_to_au       = 1.0 / au_to_fs       # au . fs^-1
au_to_as       = au_to_si_time * s_to_as # as . au^-1
as_to_au       = 1.0 / au_to_as       # au . as^-1
