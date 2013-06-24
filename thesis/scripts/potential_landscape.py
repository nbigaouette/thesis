#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
import numpy as np
import scipy.special as special
import math
import copy
import sys

import plot
import plot_params

import colors_and_symbols
import constants as cst


# ****************************************************************************************************************************************************
# http://docs.python.org/library/optparse.html
from optparse import OptionParser
parser = OptionParser()
parser.add_option(      "--two",                    action="store_true",dest="twofigures",  default=False,              help="Two figures (no subplots). [default: %default]")
parser.add_option(      "--ionization",             action="store_true",dest="ionization",  default=False,              help="Plot ionization. [default: %default]")
parser.add_option(      "--recomb",                 action="store_true",dest="recomb",      default=False,              help="Plot recombination. [default: %default]")
parser.add_option("-i", "--ion",        type=str,   action="append",    dest="ions",        default=[],                 help="List of ions <pos,cs> [default: %default]")
parser.add_option("-d", "--depth",      type=float,                     dest="depth",       default=1.0,                help="Potential depth [default: %default Hartree]")
parser.add_option("-K", "--impe_K",     type=float,                     dest="impe_K",      default=None,               help="Impacting electron's kinetic energy [default: %default Hartree]")
parser.add_option("-E", "--impe_E",     type=float,                     dest="impe_E",      default=None,               help="Impacting electron's total energy [default: %default Hartree]")
parser.add_option("-r", "--impe_r",     type=float,                     dest="impe_r",      default=1.0,                help="Impacting electron's position [default: %default Bohr]")
parser.add_option(      "--ion_impe",   type=int,                       dest="ion_impe",    default=0,                  help="Ion impacted on [default: %default]")
parser.add_option(      "--rmin",       type=float,                     dest="rmin",        default=None,               help="Lower r to plot [default: %default Bohr]")
parser.add_option(      "--rmax",       type=float,                     dest="rmax",        default=None,               help="Upper r to plot [default: %default Bohr]")
parser.add_option(      "--vmin",       type=float,                     dest="vmin",        default=None,               help="Lower potential to plot [default: %default Volt]")
parser.add_option(      "--vmax",       type=float,                     dest="vmax",        default=None,               help="Upper potential to plot [default: %default Volt]")
parser.add_option(      "--umin",       type=float,                     dest="umin",        default=None,               help="Lower energy to plot [default: %default Hartree]")
parser.add_option(      "--umax",       type=float,                     dest="umax",        default=None,               help="Upper energy to plot [default: %default Hartree]")
parser.add_option(      "--Ub",         type=float,                     dest="Ub",          default=None,               help="Force a value for Ub [default: %default Hartree]")
parser.add_option(      "--plot",       type=str,                       dest="toplot",      default="all",              help="Which plot (all,potential,energy,ip). [default: %default]")
(options, args) = parser.parse_args()
# ****************************************************************************************************************************************************

if (options.ionization == options.recomb):
    print "ERROR: Please specify either --ionization or --recomb"
    sys.exit(0)

if ((options.impe_E == None and options.impe_K == None) or (options.impe_E != None and options.impe_K != None)):
    print "ERROR: Please provide either --impe_E OR --impe_K."
    sys.exit(0)

if (len(options.ions) == 0):
    options.ions.append("-5,1")
    options.ions.append("5,1")

plot_V  = False
plot_U  = False
plot_Ip = False
if (options.toplot == "all"):
    plot_V  = True
    plot_U  = True
    plot_Ip = True
elif (options.toplot == "potential" or options.toplot == "pot" or options.toplot == "V"):
    plot_V  = True
elif (options.toplot == "energy" or options.toplot == "U"):
    plot_U  = True
elif (options.toplot == "ip"):
    plot_Ip  = True
else:
    print "ERROR: Non-valid option for --plot:", options.toplot
    import sys
    sys.exit(0)

if (not options.twofigures and (plot_U or plot_V)):
    print "ERROR: Please provide --plot=all when not using --two"
    import sys
    sys.exit(0)


assert(options.depth > 0.0)

def print_V(name, V):
    print str("%-28s" % name), "=", str("%11.7g" % (V * cst.Eh_to_eV)), "Volt ( = " + str("%11.7g" % (V)) + " Eh)"
def print_E(name, E):
    print str("%-28s" % name), "=", str("%11.7g" % E), "Eh =", str("%11.7g" % (E * cst.Eh_to_eV)), "eV"

def find_nearest(array, value):
    if (type(value) == str or type(value) == np.string_):
        for i in xrange(len(array)):
            if (array[i] == value):
                return array[i], i
    else:
        # http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
        idx = (np.abs(array-value)).argmin()
        return array[idx], idx
#

def sigma(potential_depth):
    return (1.0/potential_depth) * math.sqrt(2.0 / math.pi)
def potential_depth(sigma):
    return (1.0/sigma) * math.sqrt(2.0 / math.pi)

class Particle:
    def __init__(self, cs, pos, depth, r, K = 0.0):
        self.cs     = int(cs)
        self.pos    = pos
        self.sigma  = sigma(depth)
        self.K      = K     # Kinetic energy
        self.V      = 0.0   # Potential
        # Find index in r where particle is located
        nothing, self.index  = find_nearest(r, self.pos)
        self.nearest= None
    def Charge(self):
        return self.cs * cst.e0
    def Set_nearest(self, _nearest):
        self.nearest= _nearest
    def Add_to_V(self, pot):
        self.Set_V(self.V + pot)
    def Set_V(self, pot):
        if (type(pot) == np.float64 or type(pot) == float):
            self.V = pot
        else:
            self.V = pot[0]
    def Get_V(self, _r):
        r = np.array(_r, ndmin = 1)
        dist = abs(r-self.pos)
        dist[dist == 0] = 1.0e-9
        return (self.cs / dist) * special.erf(dist / (self.sigma * math.sqrt(2.0)))
    def U(self):
        return self.cs * self.V
    def E(self):
        return self.U() + self.K
    def Print(self):
        print   "    cs =", self.cs
        print   "    r =", self.pos, "Bohr"
        print   "    sigma =", self.sigma, "Bohr"
        print   "    depth =", potential_depth(self.sigma), "Hartree"
        print   "    V =", self.V, "au =", self.V * cst.au_to_si_pot, "Volt"
        print   "    nearest =", self.nearest
        print_E("    U", self.U())
        print_E("    K", self.K)
        print_E("    E", self.E())

class arrow:
    def __init__(self, name, start, length, is_energy):
        self.name   = name
        self.start  = start
        self.length = length
        self.stop   = self.start + self.length
        self.is_energy = is_energy
    def PrintHeader(self):
        print "Arrow name           Start         Stop         Length            Start          Stop         Length"
    def Print(self):
        if (self.is_energy):
            print str('%-16s' % self.name) + " " + \
                        str('%8.4f' % self.start) + " Eh   " + \
                        str('%8.4f' % self.stop) + " Eh   " + \
                        str('%8.4f' % self.length) + " Eh   " + \
                  "    " + \
                        str('%8.4f' % (self.start * cst.Eh_to_eV)) + " eV   " + \
                        str('%8.4f' % (self.stop * cst.Eh_to_eV)) + " eV   " + \
                        str('%8.4f' % (self.length * cst.Eh_to_eV)) + " eV"
        else:
            print str('%-16s' % self.name) + " " + \
                        str('%8.4f' % (self.start * cst.Eh_to_eV)) + " Volt " + \
                        str('%8.4f' % (self.stop * cst.Eh_to_eV)) + " Volt " + \
                        str('%8.4f' % (self.length * cst.Eh_to_eV)) + " Volt " + \
                  "   (" + \
                        str('%8.4f' % (self.start)) + " Eh   " + \
                        str('%8.4f' % (self.stop)) + " Eh   " + \
                        str('%8.4f' % (self.length)) + " Eh)"
    def Is_Null(self):
        if (self.start == self.stop):
            return True
        elif (abs(self.start - self.stop) < 1.0e-5):
            return True
        else:
            return False
    def Plot(self, ax, x, color = 'k', label = '', horizontalalignment = 'center'):
        if (self.Is_Null()):
            label_equal_zero = " = 0"
        else:
            label_equal_zero = ""
            ann = ax.annotate('', xy     = (x, self.stop),  xycoords='data',
                                  xytext = (x, self.start), textcoords='data',
                                  arrowprops=dict(arrowstyle="-|>", color = color),
                                  annotation_clip=False )
            ann.arrow_patch.set_clip_box(ax.bbox)
            ax.plot([x], [self.start], 'o' + color)
        ax.text(x, (self.start + self.stop) / 2.0,
                label + label_equal_zero,
                horizontalalignment = horizontalalignment,
                verticalalignment = "center",
                color = color)



nb_ions = len(options.ions)

# Find range of r
ions_rmin = +1.0e9
ions_rmax = -1.0e9
for i in xrange(nb_ions):
    option = options.ions[i].split(",")
    pos = float(option[0])
    cs  = int(option[1])
    if (pos < ions_rmin):
        ions_rmin = pos
    if (pos > ions_rmax):
        ions_rmax = pos
# Check impacting electron too
if (options.impe_r < ions_rmin):
    ions_rmin = options.impe_r
if (options.impe_r > ions_rmax):
    ions_rmax = options.impe_r
# And also check --rmin & --rmax
if (options.rmin != None and options.rmin < ions_rmin):
    ions_rmin = options.rmin
if (options.rmax != None and options.rmax > ions_rmax):
    ions_rmax = options.rmax

rr = ions_rmax - ions_rmin
r = np.linspace(ions_rmin - 0.1*rr, ions_rmax + 0.1*rr, 10000)
rr = r[-1] - r[0]

if (options.impe_K != None):
    impact_electron_K = options.impe_K
else:
    impact_electron_K = 0.0
    impact_electron_E = options.impe_E

# Create arrays of ions and all_particles
ions = []
for i in xrange(nb_ions):
    option = options.ions[i].split(",")
    pos = float(option[0])
    cs  = int(option[1])
    ions.append(Particle(cs = cs, pos = pos, depth = options.depth, r = r))

impacting_electron = Particle(cs = -1.0, pos = options.impe_r, depth = options.depth, r = r, K = impact_electron_K)
del impact_electron_K # Make sure we use impacting_electron.K from now on.
impacting_electron.Set_nearest(options.ion_impe)

all_particles = ions + [impacting_electron]
nb_particles = len(all_particles)

# For Xenon
GSIS = cst.eV_to_Eh * np.array([12.265625, 21.390625, 31.296875, 41.859375, 55.015625, 66.906250, 93.375000, 106.953125, 179.343750, 204.453125, 230.375000, 257.046875])
Ip_ionization    = GSIS[ions[impacting_electron.nearest].cs]
Ip_recombination = GSIS[ions[impacting_electron.nearest].cs-1]
if (options.ionization):
    Ip = Ip_ionization
    Ip_text = "ion"
else:
    assert(options.recomb)
    Ip = Ip_recombination
    Ip_text = "recomb"

# Set potential at each particles position
for i in xrange(nb_particles):
    for j in xrange(nb_particles):
        if (i != j):
            all_particles[i].Add_to_V(all_particles[j].Get_V(all_particles[i].pos))

# Now that we have the potential at each particle's location, set the impacting electron's
# kinetic energy if --impe_E was passed.
if (options.impe_K == None):
    impacting_electron.K = options.impe_E - impacting_electron.U()

# Store the total potential due to all particles
Vtot = np.zeros(len(r))
for i in xrange(nb_particles):
    Vtot += all_particles[i].Get_V(r)

# Store the total potential due to all particles without nearby electrons
Vtot_NoE = copy.deepcopy(Vtot)
for i in xrange(nb_particles):
    if (i != impacting_electron.nearest):
        dr = abs(all_particles[i].pos - all_particles[impacting_electron.nearest].pos)
        if (all_particles[i].cs == -1 and dr <= 4.0):
            Vtot_NoE -= all_particles[i].Get_V(r)

# Find Vb (potential at ion's location, without nearby electrons)
# Equivalent to the electron's total energy if there was only one ion
Vb = 0.0
Vbr = np.zeros(len(r)) # # Vb == Vbr @ impacted-on ion's position
Ubr = np.zeros(len(r)) # # Ub == Ubr @ impacted-on ion's position
Vbr_label = r'$V_{b}(r) = V_{tot}'
Ubr_label = r'$U_{b}(r) = U_{tot}'
if (options.Ub):
    Ub = options.Ub
    Vb = -options.Ub
    Vbr = Vb * np.ones(len(r))
    Ubr = Ub * np.ones(len(r))
    Vbr_label = Vbr_label + r'$'
    Ubr_label = Ubr_label + r'$'
else:
    for i in xrange(nb_particles):
        # Distance between the two particles
        dr = abs(all_particles[i].pos - all_particles[impacting_electron.nearest].pos)
        # Remove impacted-on ion and electrons nearby (4 bohr)
        if (not ((i == impacting_electron.nearest) or (all_particles[i].cs == -1 and dr <= 4.0))):
            Vb  += all_particles[i].Get_V(all_particles[impacting_electron.nearest].pos)
            Vbr += all_particles[i].Get_V(r)
            Ubr -= all_particles[i].Get_V(r)
        else:
            Vbr_label = Vbr_label + ' - V_{' + str(i) + '}'
            Ubr_label = Ubr_label + ' - U_{' + str(i) + '}'
    Vbr_label = Vbr_label + '$'
    Ubr_label = Ubr_label + '$'
    Ub = -1.0 * Vb # Carefull with the Potential <--> Energy conversion here.

# Calculate kinetic energy of the impacting electron with respect to the threshold
# First, send the electron to infinity since cross-sections requires that value
V_WImpE_DIon = ions[impacting_electron.nearest].Get_V(impacting_electron.pos)
U_ImpE_Ion = -1.0 * V_WImpE_DIon
impe_K_inf_ion = max(0.0, impacting_electron.K + U_ImpE_Ion)

assert(Ip > 0.0)

# Then, adapt for Vb (the new threshold caused by the cluster environment)
K_thres = max(0.0, impacting_electron.K - (Ub - impacting_electron.U()))

# Difference between the threshold value "Ub" and the potential energy curve due to the cluster
DeltaU =  (impacting_electron.U() - U_ImpE_Ion) - Ub

Ee_recomb = impacting_electron.K + U_ImpE_Ion

arrow_ImpE_Uion     = arrow(name = "ImpE_Uion",     start = 0.0,                    length = U_ImpE_Ion,            is_energy = True)
arrow_ImpE_Kinf1    = arrow(name = "ImpE_Kinf (1)", start = U_ImpE_Ion,             length = impe_K_inf_ion,        is_energy = True)
arrow_ImpE_Kinf2    = arrow(name = "ImpE_Kinf (2)", start = 0.0,                    length = impe_K_inf_ion,        is_energy = True)
arrow_Vb            = arrow(name = "Vb",            start = 0.0,                    length = Vb,                    is_energy = False)
arrow_Ub            = arrow(name = "Ub",            start = 0.0,                    length = Ub,                    is_energy = True)
arrow_ImpE_K        = arrow(name = "ImpE_K",        start = impacting_electron.U(), length = impacting_electron.K,  is_energy = True)
arrow_ImpE_U        = arrow(name = "ImpE_U",        start = 0.0,                    length = impacting_electron.U(),is_energy = True)
arrow_ImpE_E        = arrow(name = "ImpE_E",        start = 0.0,                    length = impacting_electron.E(),is_energy = True)
arrow_ImpE_Kshifted = arrow(name = "ImpE_K_shifted",start = Ub-V_WImpE_DIon,        length = impacting_electron.K,  is_energy = True)
arrow_ImpE_Eshifted = arrow(name = "ImpE_E_shifted",start = Ub,                     length = Ee_recomb,             is_energy = True)
arrow_Ip            = arrow(name = "Ip",            start = -abs(Ip),               length = abs(Ip),               is_energy = True)
arrow_Ip_prime      = arrow(name = "Ip'",           start = Ub-abs(Ip),             length = abs(Ip),               is_energy = True)

arrow_ImpE_K_thresh = arrow(name = "ImpE_K_thres",  start = Ub,                     length = K_thres,               is_energy = True)
arrow_ImpE_U_thresh = arrow(name = "ImpE_U_thres",  start = Ub,                     length = -V_WImpE_DIon,         is_energy = True)

arrow_DeltaU        = arrow(name = "DeltaU",        start = Ub,                     length = DeltaU,               is_energy = True)

if (options.ionization):
    can_recombine = False
    if (K_thres >= Ip):
        can_ionize = True
        can_ionize_yes_no = "yes (K_threshold >= Ip)"
    else:
        can_ionize = False
        can_ionize_yes_no = "no (K_threshold < Ip)"
else:
    can_ionize = False
    can_ionize_yes_no = ""

    # Is recombination possible?
    if (Ee_recomb <= -Ip):
        can_recombine = True
    else:
        can_recombine = False

# ******************************************************************************
# Print values
# ******************************************************************************

print "potential_depth =", options.depth * cst.Eh_to_eV, "eV =", options.depth, "Hartree"
print "Impacting electron kinetic energy =", impacting_electron.K*cst.Eh_to_eV, "eV =", impacting_electron.K, "Hartree"

print "Number of ions:", nb_ions
print "Ions info:"
for i in xrange(nb_ions):
    print "  Ion", i
    ions[i].Print()
print ""
print "Impacting electron:"
impacting_electron.Print()

print ""
print   "can_ionize                   =", can_ionize
print   "can_recombine                =", can_recombine
print_E("Ke", impacting_electron.K)
print_E("Ip (recombination)", Ip_recombination)
print_E("Ip (ionization)", Ip_ionization)
print_E("Ip", Ip)
print_V("V_WImpE_DIon", V_WImpE_DIon)
print_E("U_ImpE_Ion",   U_ImpE_Ion)
print_V("Vb", Vb)
print_E("Ub", Ub)
print_E("impe_K_inf_ion", impe_K_inf_ion)
print_E("K_thres", K_thres)
print_E("K_thres-Ip", K_thres-Ip)
print_E("Ub - (U_ImpE - U_ImpE_Ion)", DeltaU)
print_E("Ee_recomb", Ee_recomb)
print ""

arrow_ImpE_K.PrintHeader()

arrow_Ip.Print()
arrow_Ip_prime.Print()
arrow_ImpE_Uion.Print()
arrow_ImpE_K.Print()
arrow_ImpE_U.Print()
arrow_ImpE_E.Print()
arrow_ImpE_Kshifted.Print()
arrow_ImpE_Eshifted.Print()
arrow_Vb.Print()
arrow_Ub.Print()
arrow_ImpE_Kinf1.Print()
arrow_ImpE_Kinf2.Print()
arrow_ImpE_K_thresh.Print()
arrow_ImpE_U_thresh.Print()
arrow_DeltaU.Print()

# ******************************************************************************
# Plots
# ******************************************************************************

if (options.twofigures):
    figs = []
    if (plot_V):
        fig1 = plot.figure()
        figs.append(fig1)
        plt.subplots_adjust(hspace = 0.0)
    if (plot_U):
        fig2 = plot.figure()
        figs.append(fig2)
        plt.subplots_adjust(hspace = 0.0)
    figs_N = 1
    fig1_y = 1
    fig2_y = 1
else:
    fig1 = plot.figure()
    fig2 = fig1
    figs = [fig1]
    figs_N = 2
    fig1_y = 1
    fig2_y = 2
    plt.subplots_adjust(hspace = 0.0)
axprops = dict()

if (plot_V):
    ax_V_Eh = SubplotHost(fig1, figs_N,1,fig1_y, **axprops)
    ax_V_Eh_to_Volt = mtransforms.Affine2D().scale(1.0, cst.eV_to_Eh)
    ax_V_Volt = ax_V_Eh.twin(ax_V_Eh_to_Volt)
    ax_V_Volt.set_viewlim_mode("transform")
    ax_V_Eh.grid(True)

    ax_V_Volt.set_ylabel("Potential (Volt)")
    ax_V_Eh.set_ylabel("Potential energy of a 1+ (Hartree)")

    axprops['sharex'] = ax_V_Volt

    plt.setp(ax_V_Volt.get_xticklabels(), visible=False)

    #ax_V_Eh.set_title(r"Potential")

if (plot_U):
    ax_U_Eh = SubplotHost(fig2, figs_N,1,fig2_y, **axprops)
    ax_U_Eh_to_eV  = mtransforms.Affine2D().scale(1.0, cst.eV_to_Eh)
    ax_U_eV = ax_U_Eh.twin(ax_U_Eh_to_eV)
    ax_U_eV.set_viewlim_mode("transform")
    ax_U_Eh.grid(True)

    ax_U_Eh.set_xlabel("Position (Bohr)")
    ax_U_Eh.set_ylabel("Potential energy of an electron (Hartree)")
    ax_U_eV.set_ylabel("Potential energy of an electron (eV)")

    #ax_U_Eh.set_title(r"Potential Energy")

    plt.setp(ax_U_eV.get_xticklabels(), visible=False)

if (not options.twofigures):
    plt.setp(ax_V_Eh.get_xticklabels(), visible=False)

if (plot_V):
    fig1.add_subplot(ax_V_Eh)
if (plot_U):
    fig2.add_subplot(ax_U_Eh)


# Plot symbols for particles
for i in xrange(nb_particles):
    if (all_particles[i].cs >= 0):
        str_elem = 'i'
        str_cs = '+'
        if (i == impacting_electron.nearest):
            symbol = 'rs'
        else:
            symbol = 'ro'
    else:
        str_elem = 'e'
        symbol = 'bo'
        str_cs = ''

    # Potential
    if (plot_V):
        ax_V_Eh.plot(all_particles[i].pos, all_particles[i].V, symbol)
        ax_V_Eh.text(all_particles[i].pos, all_particles[i].V, r' $' + str_elem + '_{' + str(i) + '}\ (' + str(all_particles[i].cs) + str_cs + ')$', horizontalalignment = "left")
    # Energy
    #if (plot_U):
    #    ax_U_Eh.plot(all_particles[i].pos, all_particles[i].U(), symbol)
    #    ax_U_Eh.text(all_particles[i].pos, all_particles[i].U(), r' $' + str_elem + '_{' + str(i) + '}\ (' + str(all_particles[i].cs) + str_cs + ')$', horizontalalignment = "center")


# ******************************************************************************
# Potential plot
# ******************************************************************************

if (plot_V):
    # Total potential
    ax_V_Eh.plot(r, Vtot, '-k', label=r'$V_{tot}$')

    # Potential due to each particle
    for i in xrange(nb_particles):
        ax_V_Eh.plot(r, all_particles[i].Get_V(r), '--' + colors_and_symbols.color(i),  label=r'$V_{' + str(i) + '}$')
        Vi = all_particles[i].V
        arrow_to_plot = arrow(name = "Vi", start = 0.0, length = Vi, is_energy = False)
        arrow_to_plot.Plot(ax_V_Eh, x = all_particles[i].pos, label = r"$V(" + str(i) + ")$", horizontalalignment="left")

    # Plot Vb
    # Total potential minus the ion we are looking at and nearby electrons (Vb function of r)
    ax_V_Eh.plot(r, Vbr, "--k", label=Vbr_label)
    if (nb_ions > 1):
        ax_V_Eh.plot([r[0], r[-1]], [Vb, Vb], '--k')
    arrow_Vb.Plot(ax_V_Eh, x = ions[impacting_electron.nearest].pos, color = 'r', label = r'$V_b$', horizontalalignment = 'left')


# ******************************************************************************
# Energy plot
# ******************************************************************************

if (plot_U):
    dr = rr*0.005

    # Potential energy landscape of an electron in the different potentials
    # All
    #ax_U_Eh.plot(r, impacting_electron.cs * Vtot, '-k', label=r'$U_{e,tot}$')
    # Without nearby electrons
    ax_U_Eh.plot(r, impacting_electron.cs * Vtot_NoE, '--k', label=r'$U_{e,tot} - \sum_{\rm nearby\ el} U_{e,i}$')


    # Potential energy of an electron due to each particle
    for i in xrange(nb_particles):
        if (all_particles[i].cs != -1):
            ax_U_Eh.plot(r, -all_particles[i].Get_V(r), '--' + colors_and_symbols.color(i),  label=r'$U_{e,' + str(i) + '}$')
        if (i == impacting_electron.nearest):
            ax_U_Eh.plot(r, -all_particles[i].Get_V(r) + Ub, ':' + colors_and_symbols.color(i),  label=r'$U_{e,' + str(i) + '} + U_b$')

    arrow_Ub.Plot(ax_U_Eh, x = all_particles[impacting_electron.nearest].pos, color = 'r', label = r'$U_b$', horizontalalignment = 'left')

    # ************************************
    # Impacting electron's...
    # ************************************

    # ************************************
    # ...U
    dr = 0.05
    arrow_ImpE_U.Plot(ax_U_Eh, x = impacting_electron.pos-dr, color = 'g', label = r'$U_e$', horizontalalignment = 'right')

    # ************************************
    # ...K
    arrow_ImpE_K.Plot(ax_U_Eh, x = impacting_electron.pos+dr, color = 'g', label = r'$K_e$', horizontalalignment = 'left')

    # ************************************
    # ...E
    arrow_ImpE_E.Plot(ax_U_Eh, x = impacting_electron.pos, color = 'g', label = r' $E_e$', horizontalalignment = 'left')

    # ************************************
    # U (with respect to ion)
    #arrow_ImpE_Uion.Plot(ax_U_Eh, x = impacting_electron.pos+2.0*dr, color = 'b', label = r'$U_e^{\rm wrt\ ion}$', horizontalalignment = 'left')
    #ax_U_Eh.plot([impacting_electron.pos+dr, impacting_electron.pos+2.0*dr], [arrow_ImpE_Uion.stop, arrow_ImpE_Uion.stop], '-b')

    # ************************************
    # K (with respect to ion at infinity)
    #if (options.ionization):
    #    label = r'   $K_{\infty}^{\rm wrt\ ion}$'
    #    #if (can_ionize):
    #        #label = label + r" >= $Ip'_{\rm ion}$"
    #    #else:
    #        #label = label + r" < $Ip_{\rm ion}$"
    #    arrow_ImpE_Kinf2.Plot(ax_U_Eh, x = impacting_electron.pos+2.0*dr, color = 'b', label = label, horizontalalignment = 'left')
    #    arrow_ImpE_Kinf1.Plot(ax_U_Eh, x = impacting_electron.pos+4.0*dr, color = 'b', label = label, horizontalalignment = 'left')

    arrow_DeltaU.Plot(ax_U_Eh, x = impacting_electron.pos+4.0*dr, color = 'c', label = r"$\Delta U$", horizontalalignment = 'left')
    ax_U_Eh.plot([impacting_electron.pos-dr, impacting_electron.pos+5.0*dr], [arrow_DeltaU.stop, arrow_DeltaU.stop], '-c')

    # ************************************
    # K (with respect to threshold)
    arrow_ImpE_K_thresh.Plot(ax_U_Eh, x = impacting_electron.pos-2.0*dr, color = 'r', label = r"$K_{\rm thresh}$", horizontalalignment = 'right')
    # U (with respect to threshold)
    arrow_ImpE_U_thresh.Plot(ax_U_Eh, x = impacting_electron.pos-2.0*dr, color = 'r', label = r"$U_{\rm " + str(impacting_electron.nearest) + "," + str(nb_particles-1) + " (thresh)}$", horizontalalignment = 'right')

    if (options.recomb):
        label = ""
        if (abs(impacting_electron.pos - ions[options.ion_impe].pos) > 5.0):
            label = label + r" > $-Ip_{\rm recomb}$ (no recomb: too far)"
        elif (can_recombine):
            label = label + r" < $-Ip_{\rm recomb}$ (recomb. possible)"
        else:
            label = label + r" > $-Ip_{\rm recomb}$ (no recomb.)"
        arrow_ImpE_Kshifted.Plot(ax_U_Eh, x = impacting_electron.pos, color = 'g', label = r'$K_e$', horizontalalignment = 'left')
        arrow_ImpE_Eshifted.Plot(ax_U_Eh, x = impacting_electron.pos, color = 'g', label = r"$E'_e$" + label, horizontalalignment = 'left')

    # ************************************
    # Plot some horizontal lines at start and stop of arrows
    ax_U_Eh.plot([impacting_electron.pos-dr, impacting_electron.pos+dr], [arrow_ImpE_U.start, arrow_ImpE_U.start], '-g')
    ax_U_Eh.plot([impacting_electron.pos-dr, impacting_electron.pos+dr], [arrow_ImpE_U.stop,  arrow_ImpE_U.stop],  '-g')
    #ax_U_Eh.plot([impacting_electron.pos-dr, impacting_electron.pos+dr], [imp_el_K_start, imp_el_K_start], '-g') # Same as imp_el_U_stop
    ax_U_Eh.plot([impacting_electron.pos-dr, impacting_electron.pos+dr], [arrow_ImpE_K.stop,  arrow_ImpE_K.stop],  '-g')
    #ax_U_Eh.plot([impacting_electron.pos-dr, impacting_electron.pos+dr], [imp_el_E_start, imp_el_E_start], '-g') # Same as imp_el_U_start
    ax_U_Eh.plot([impacting_electron.pos-dr, impacting_electron.pos+dr], [arrow_ImpE_E.stop,  arrow_ImpE_E.stop],  '-g')
    assert(arrow_ImpE_U.start == 0.0)
    assert(arrow_ImpE_U.stop  == arrow_ImpE_K.start)
    assert(arrow_ImpE_K.stop  == arrow_ImpE_E.stop)
    assert(arrow_ImpE_E.start == 0.0)
    #if (impe_K_inf_ion > 0.0):
    #    ax_U_Eh.plot([impacting_electron.pos-dr, impacting_electron.pos+4.0*dr], [arrow_ImpE_Kinf1.start, arrow_ImpE_Kinf1.start], '-b')
    #    ax_U_Eh.plot([impacting_electron.pos-dr, impacting_electron.pos+4.0*dr], [arrow_ImpE_Kinf1.stop,  arrow_ImpE_Kinf1.stop],  '-b')
    #    ax_U_Eh.plot([impacting_electron.pos-dr, impacting_electron.pos+2.0*dr], [arrow_ImpE_Kinf2.start, arrow_ImpE_Kinf2.start], '-b')
    #    ax_U_Eh.plot([impacting_electron.pos-dr, impacting_electron.pos+2.0*dr], [arrow_ImpE_Kinf2.stop,  arrow_ImpE_Kinf2.stop],  '-b')

    # ************************************
    # Ub(r)
    #ax_U_Eh.plot(r, Ubr, "-.m", label=Ubr_label)

    # ************************************
    # Ip
    arrows      = [arrow_Ip, arrow_Ip_prime]
    alignments  = ["right", "left"]
    which_Ip    = ["", "'"]
    Ip_drs      = [-0.015*rr,   0.015*rr]
    ll = 1.0
    for i in xrange(len(arrows)):
        arrows[i].Plot(ax_U_Eh, x = ions[impacting_electron.nearest].pos+Ip_drs[i], color = 'm', label = r'$Ip' + which_Ip[i] + r'_{\rm ' + Ip_text + '}$', horizontalalignment = alignments[i])
        ax_U_Eh.plot([ions[impacting_electron.nearest].pos - ll, ions[impacting_electron.nearest].pos + ll],
                     [arrows[i].start, arrows[i].start], '--m')
        ax_U_Eh.axhline(arrows[i].stop, linestyle=':', color='m')


# ******************************************************************************
# Create a third plot to compare K_threshold and Ip
if (plot_Ip):
    fig3 = plot.figure()
    figs.append(fig3)

    ax3_Eh = SubplotHost(fig3, 1,1,1)
    ax3_Eh_to_eV = mtransforms.Affine2D().scale(1.0, cst.eV_to_Eh)
    ax3_eV = ax3_Eh.twin(ax3_Eh_to_eV)
    ax3_eV.set_viewlim_mode("transform")

    ax3_Eh.grid(True)
    ax3_Eh.grid(True)

    plt.setp(ax3_Eh.get_xticklabels(), visible=False)
    plt.setp(ax3_eV.get_xticklabels(), visible=False)

    fig3.add_subplot(ax3_Eh)

    ax3_Eh.axhline(0.0, linestyle='--', color='k')

    start = 0.0

    # Ip
    Ip_dr = -0.1

    arrow_Ip = arrow(name = "Ip", start = start, length = Ip, is_energy = True)
    arrow_Ip.Plot(ax3_Eh, x = +1.0*Ip_dr, color = 'm', label = "$Ip$", horizontalalignment = "right")
    ax3_Eh.axhline(arrow_Ip.stop, linestyle='--', color='m', label = "$Ip$")

    # K
    arrow_K = arrow(name = "K", start = start, length = impacting_electron.K, is_energy = True)
    arrow_K.Plot(ax3_Eh, x = -1.0*Ip_dr, color = 'g', label = "$K$", horizontalalignment = "left")
    ax3_Eh.axhline(arrow_K.stop, linestyle='--', color='g', label = "Impacting electron K")

    # K_ion_inf
    arrow_Kii = arrow(name = "K_ino_inf", start = start, length = impe_K_inf_ion,  is_energy = True)
    arrow_Kii.Plot(ax3_Eh, x = -2.0*Ip_dr, color = 'b', label = "$K_{e,ion\ \infty}$", horizontalalignment = "left")
    ax3_Eh.axhline(arrow_Kii.stop, linestyle='--', color='b', label = "K w.r.t. parent ion at infinity")

    # K_threshold
    arrow_Kt = arrow(name = "K_threshold", start = start, length = K_thres,  is_energy = True)
    arrow_Kt.Plot(ax3_Eh, x = -3.0*Ip_dr, color = 'r', label = "$K_{threshold}$", horizontalalignment = "left")
    ax3_Eh.axhline(arrow_Kt.stop, linestyle='--', color='r', label = "K w.r.t. threshold ($V_b$)")

    # Can it ionize?
    ax3_Eh.set_title(r"Can ionize: " + can_ionize_yes_no)

    ax3_Eh.set_ylabel("Hartree")
    ax3_eV.set_ylabel("eV")

    ax3_Eh.set_xlim((-0.2, 0.4))
    ax3_Eh.set_ylim((-0.5, 1.1*max(arrow_Ip.stop, arrow_K.stop)))
    leg = ax3_Eh.legend(loc="best")
    leg.get_frame().set_alpha(0.4)

# ******************************************************************************
# Set plot properties
# ******************************************************************************

# Default the xlim be r[0] and r[-1]
if (plot_V):
    ax_V_Eh.set_xlim((r[0], r[-1]))
if (plot_U):
    ax_U_Eh.set_xlim((r[0], r[-1]))

# Make sure the ions' label are visible
if (plot_V):
    ylim = ax_V_Volt.get_ylim()
    ax_V_Volt.set_ylim((min(-0.1, ylim[0]), ylim[1]))
if (plot_U):
    ylim = ax_U_Eh.get_ylim()
    ax_U_Eh.set_ylim((ylim[0], max(0.1, ylim[1])))

if (options.rmin != None or options.rmax != None):
    if (plot_V or plot_U):
        if (plot_V):
            xlim = list(ax_V_Volt.get_xlim())
        elif (plot_U):
            xlim = list(ax_U_Eh.get_xlim())
        if (options.rmin != None):
            xlim[0] = options.rmin
        if (options.rmax != None):
            xlim[1] = options.rmax
        if (plot_V):
            ax_V_Eh.set_xlim(xlim)
        if (plot_U):
            ax_U_Eh.set_xlim(xlim)

if (options.vmin != None or options.vmax != None):
    if (plot_V):
        ylim = list(ax_V_Eh.get_ylim())
        if (options.vmin != None):
            ylim[0] = options.vmin
        if (options.vmax != None):
            ylim[1] = options.vmax
        ax_V_Eh.set_ylim(ylim)

if (options.umin != None or options.umax != None):
    if (plot_U):
        ylim = list(ax_U_eV.get_ylim())
        if (options.umin != None):
            ylim[0] = options.umin
        if (options.umax != None):
            ylim[1] = options.umax
        ax_U_Eh.set_ylim(ylim)

if (plot_V):
    leg = ax_V_Eh.legend(loc="best")
    leg.get_frame().set_alpha(0.4)
#if (plot_U):
    #leg = ax_U_Eh.legend(loc="best")
    #leg.get_frame().set_alpha(0.4)

#plot.savefig('potential_landscape.svg')
#plot.savefig('potential_landscape.pdf')
plot.show()
