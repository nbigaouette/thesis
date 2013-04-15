#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# ****************************************************************************************************************************************************
# http://docs.python.org/library/optparse.html
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-w", "--wavelength", type=float,         dest="wavelength",  default=None,       help="Wavelength [default: %default]")
parser.add_option("-e", "--energy",     type=float,         dest="energy",      default=None,       help="Energy [default: %default]")
parser.add_option(      "--ue",         type=str,           dest="unit_energy", default="eV",       help="Energy units [default: %default]")
parser.add_option(      "--uw",         type=str,           dest="unit_wl",     default="nm",       help="Wavelength units [default: %default]")
(options, args) = parser.parse_args()
# ****************************************************************************************************************************************************

import sys
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

Eh_to_eV        = Eh / e0               # eV . Eh^-1
eV_to_Eh        = 1.0 / Eh_to_eV        # Eh . eV^-1
eV_to_J         = e0                    # J . eV^-1
J_to_eV         = 1.0 / eV_to_J         # eV . J^-1
J_to_Eh         = eV_to_Eh * J_to_eV    # Eh . J^-1
Eh_to_J         = 1.0 / J_to_Eh         # J . Eh^-1


if (options.wavelength != None):
    wavelength = options.wavelength
    if (options.unit_wl == 'm'):
        pass
    elif (options.unit_wl == 'nm'):
        wavelength *= nm_to_m
    elif (options.unit_wl == 'A' or options.unit_wl == 'angstrom'):
        wavelength *= angstrom_to_m
    else:
        print 'Unknown units for wavelength (' + options.unit_wl + ')'
        sys.exit(0)

    k       = 2.0 * pi / wavelength
    omega   = k * c
    f       = omega / (2.0 * pi)
    gamma   = hbar * omega

    print "Wavelength to energy:"
    print "    {: 12g} m           {: 12g} J".format(wavelength, gamma)
    print "    {: 12g} nm   --->   {: 12g} eV".format(wavelength * m_to_nm, gamma * J_to_eV)
    print "    {: 12g} A           {: 12g} Eh".format(wavelength * m_to_angstrom, gamma * J_to_Eh)
    print "                 -->        f = {: 12g} THz".format(f     * 1.0e-12)
    print "                 -->    omega = {: 12g} THz".format(omega * 1.0e-12)


if (options.energy != None):
    energy = options.energy
    if (options.unit_energy == 'J'):
        pass
    elif (options.unit_energy == 'eV'):
        energy *= eV_to_J
    elif (options.unit_energy == 'Eh' or options.unit_energy == 'Hartree'):
        energy *= Eh_to_J
    else:
        print 'Unknown units for energy (' + options.unit_energy + ')'
        sys.exit(0)

    omega       = energy / hbar
    f       = omega / (2.0 * pi)
    k           = omega / c
    wavelength  = 2.0 * pi / k

    print "Energy to wavelength:"
    print "    {: 12g} J           {: 12g} m".format(energy, wavelength)
    print "    {: 12g} eV   --->   {: 12g} nm".format(energy * J_to_eV, wavelength * m_to_nm)
    print "    {: 12g} Eh          {: 12g} Angstrom".format(energy * J_to_Eh, wavelength * m_to_angstrom)
    print "                 -->        f = {: 12g} THz".format(f     * 1.0e-12)
    print "                 -->    omega = {: 12g} THz".format(omega * 1.0e-12)

