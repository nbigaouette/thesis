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

import constants as cst


if (options.wavelength != None):
    wavelength = options.wavelength
    if (options.unit_wl == 'm'):
        pass
    elif (options.unit_wl == 'nm'):
        wavelength *= cst.nm_to_m
    elif (options.unit_wl == 'A' or options.unit_wl == 'angstrom'):
        wavelength *= cst.angstrom_to_m
    else:
        print 'Unknown units for wavelength (' + options.unit_wl + ')'
        sys.exit(0)

    k       = 2.0 * cst.pi / wavelength
    omega   = k * cst.c
    f       = omega / (2.0 * pi)
    gamma   = cst.hbar * omega

    print "Wavelength to energy:"
    print "    {: 12g} m           {: 12g} J".format(wavelength, gamma)
    print "    {: 12g} nm   --->   {: 12g} eV".format(wavelength * cst.m_to_nm, gamma * cst.J_to_eV)
    print "    {: 12g} A           {: 12g} Eh".format(wavelength * cst.m_to_angstrom, gamma * cst.J_to_Eh)
    print "                 -->        f = {: 12g} THz".format(f     * 1.0e-12)
    print "                 -->    omega = {: 12g} THz".format(omega * 1.0e-12)


if (options.energy != None):
    energy = options.energy
    if (options.unit_energy == 'J'):
        pass
    elif (options.unit_energy == 'eV'):
        energy *= cst.eV_to_J
    elif (options.unit_energy == 'Eh' or options.unit_energy == 'Hartree'):
        energy *= cst.Eh_to_J
    else:
        print 'Unknown units for energy (' + options.unit_energy + ')'
        sys.exit(0)

    omega       = energy / cst.hbar
    f           = omega / (2.0 * cst.pi)
    k           = omega / cst.c
    wavelength  = 2.0 * cst.pi / k

    print "Energy to wavelength:"
    print "    {: 12g} J           {: 12g} m".format(energy, wavelength)
    print "    {: 12g} eV   --->   {: 12g} nm".format(energy * cst.J_to_eV, wavelength * cst.m_to_nm)
    print "    {: 12g} Eh          {: 12g} Angstrom".format(energy * cst.J_to_Eh, wavelength * cst.m_to_angstrom)
    print "                 -->        f = {: 12g} THz".format(f     * 1.0e-12)
    print "                 -->    omega = {: 12g} THz".format(omega * 1.0e-12)

