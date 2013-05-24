#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import sys, os, glob, math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
import matplotlib.transforms as mtransforms

import netCDF4

import constants as cst
import on_key
import read_data
import read_xml
import get_simulations_parameters
import colors_and_symbols

# ****************************************************************************************************************************************************
# http://docs.python.org/library/optparse.html
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-i", "--input",                          type=str,   dest="folder",      default="output",   help="Folder contianing simulations data")
parser.add_option("-p", "--pattern",    action="append",    type=str,   dest="patterns",    default=None,       help="Input patterns")
parser.add_option("-r", "--reload",     action="store_true",            dest="reload",      default=False,      help="Reload data. [default: %default]")
parser.add_option("-m", "--max",                            type=int,   dest="max_plot",    default=10,         help="Maximum number of plot lines [default: %default]")
parser.add_option("-s", "--save",       action="store_true",            dest="save_figure", default=False,      help="Save figure [default: %default]")
(options, args) = parser.parse_args()
# ****************************************************************************************************************************************************

import savefigure
if (options.save_figure):
    savefigure.prepare_for_save()
#

if (options.patterns == None):
    options.patterns = ["*"]

def extract_data(input_folder, input_patterns, base_potentials, potential_shapes, dts):

    # Change in total energy before and after ionization event
    NumericalHeating= np.zeros((len(base_potentials), len(potential_shapes), len(dts)))
    E_t0            = np.zeros((len(base_potentials), len(potential_shapes), len(dts)))
    E_t1            = np.zeros((len(base_potentials), len(potential_shapes), len(dts)))
    E_delta         = np.zeros((len(base_potentials), len(potential_shapes), len(dts)))

    for input_pattern in input_patterns:
        for base_potential in base_potentials:
            for potential_shape in potential_shapes:
                for dt in dts:
                    pattern_string = "*" + input_pattern + "*_dt" + str("%08.4f" % dt) + "as*_pot" + str(potential_shape) + "*_b" + str("%05.2f" % float(base_potential)) + "*"
                    globber = glob.glob(os.path.join(input_folder, pattern_string))
                    if (len(globber) == 0):
                        print "Can't find", os.path.join(input_folder, pattern_string)
                        continue
                    else:
                        basename = globber[0]
                    print basename

                    nothing, index0 = read_data.find_nearest(base_potentials,  base_potential)
                    nothing, index1 = read_data.find_nearest(potential_shapes, potential_shape)
                    nothing, index2 = read_data.find_nearest(dts,              dt)

                    # Skip folder with non-existing files
                    if (not os.path.exists(os.path.join(basename, "dump_input.xml"))):
                        continue
                    if (not os.path.exists(os.path.join(basename, "energy.bin"))):
                        continue

                    xml_input_file = read_xml.xml_content(basename, print_filename = False)

                    photon_energy = cst.hbar * (2.0 * math.pi / (xml_input_file.laser_wavelength * cst.nm_to_m) * cst.co) * cst.J_to_eV # [eV]
                    Ip = xml_input_file.element.first_Ip # [eV]

                    # Read energy file
                    filename = os.path.join(basename, "energy.bin")
                    #print "Reading file \"" + filename + "\""
                    data = read_data.read_binary_file(filename, dtype_float = xml_input_file.dtype_float)

                    nb_columns = 5
                    data = data[0:int(math.floor(float(len(data)) / float(nb_columns)))*nb_columns]
                    nb_time_steps = len(data) / nb_columns
                    data = data.reshape((nb_time_steps,nb_columns))

                    # And assign columns to variables
                    t        = data[:,0]            # Time (fs)
                    Energy_K = data[:,1]            # Total Kinetic Energy (eV)
                    Energy_U = data[:,2]            # Total Potential Energy (eV)
                    Energy_D = -data[:,3]           # Energy delta due to ionization/recombination (eV)
                    #Energy_Heating = data[:,4]      # Numerical heating
                    Energy_T = Energy_U+Energy_K    # Total Energy

                    Energy_T0 = Energy_T[0]  + Energy_D[0]
                    # Average the last 10%
                    indices = np.where(t > t[-1]/2.0)
                    Energy_T1 = (Energy_T[indices] + Energy_D[indices]).mean()
                    Energy_Heating = Energy_T1 - Energy_T0

                    E_t0[index0, index1, index2]            = Energy_T0
                    E_t1[index0, index1, index2]            = Energy_T1
                    E_delta[index0, index1, index2]         = Energy_D[-1]
                    NumericalHeating[index0, index1, index2]= Energy_Heating

                    del t, data, Energy_K, Energy_U, Energy_D, Energy_Heating

    return E_t0, E_t1, E_delta, NumericalHeating

saved_data_folder = "saved_data"
saved_file    = os.path.join(saved_data_folder, "numerical_heating.npz")

if (options.reload):
    if (not os.path.exists(saved_file)):
        print "Can't find " + saved_file + ", Exiting."
        sys.exit(0)
    saved_data = np.load(saved_file)
    base_potentials = saved_data['base_potentials']
    potential_shapes= saved_data['potential_shapes']
    dts             = saved_data['dts']
    E_t0            = saved_data['E_t0']
    E_t1            = saved_data['E_t1']
    E_delta         = saved_data['E_delta']
    NumericalHeating= saved_data['NumericalHeating']
else:
    base_potentials = np.array(get_simulations_parameters.get_all_parameters(options.folder, options.patterns, "b"), dtype=float)  # Hartree
    potential_shapes=          get_simulations_parameters.get_all_parameters(options.folder, options.patterns, "pot")
    dts             = np.array(get_simulations_parameters.get_all_parameters(options.folder, options.patterns, "dt"), dtype=float) # as

    E_t0, E_t1, E_delta, NumericalHeating = extract_data(options.folder, options.patterns, base_potentials, potential_shapes, dts)

    print "Saving to", saved_file
    np.savez(saved_file,    base_potentials=base_potentials,
                            potential_shapes=potential_shapes,
                            dts=dts,
                            E_t0=E_t0,
                            E_t1=E_t1,
                            E_delta=E_delta,
                            NumericalHeating=NumericalHeating)

print "base_potentials =", base_potentials
print "potential_shapes =", potential_shapes
print "dts =", dts

fig1 = on_key.figure()
fig2 = on_key.figure()
figs = [fig1, fig2]

ax1_eV = SubplotHost(fig1, 1,1,1)
ax1_eV_to_Eh = mtransforms.Affine2D().scale(1.0, cst.Eh_to_eV)
ax1_Eh = ax1_eV.twin(ax1_eV_to_Eh)
ax1_Eh .set_viewlim_mode("transform")
fig1.add_subplot(ax1_eV)

ax2_eV = SubplotHost(fig2, 1,1,1)
ax2_eV_to_Eh = mtransforms.Affine2D().scale(cst.eV_to_Eh, cst.Eh_to_eV)
ax2_Eh = ax2_eV.twin(ax2_eV_to_Eh)
ax2_Eh .set_viewlim_mode("transform")
fig2.add_subplot(ax2_eV)


# Plot NumericalHeating as a function of dt for every potential depth
dbase_potentials = (base_potentials[-1] - base_potentials[0]) / float(options.max_plot-1) # -1 since we want the number of intervals
c = 0
for j in xrange(len(potential_shapes)):
    for base_potentials_close in np.arange(base_potentials[0], base_potentials[-1]+dbase_potentials/2.0, dbase_potentials):
        nothing, index0 = read_data.find_nearest(base_potentials,  base_potentials_close)
        nothing, index1 = read_data.find_nearest(potential_shapes, potential_shapes[j])
        #nothing, index2 = read_data.find_nearest(dts,              dts[i])
        ax1_eV.plot(dts, NumericalHeating[index0, index1, :],
                    colors_and_symbols.symb_col(c),
                    label = r"Potential depth: " + str(base_potentials[index0]) + " Eh (" + potential_shapes[j] + " )")
        # When the NumericalHeating is negative, we have cooling, which would not appear on the log scale plot.
        cooling_indices = np.where(NumericalHeating[index0, index1, :] <= 0.0)
        if (len(cooling_indices[0]) >= 1):
            ax1_eV.plot(dts[cooling_indices], abs(NumericalHeating[index0, index1, :][cooling_indices]),
                        '.' + colors_and_symbols.symb_col(c))
        c += 1

# Plot NumericalHeating as a function of potential depth for every dt
ddt = (dts[-1] - dts[0]) / float(options.max_plot-1) # -1 since we want the number of intervals
c = 0
for j in xrange(len(potential_shapes)):
    for dt_close in np.arange(dts[0], dts[-1]+ddt/2.0, ddt):
        nothing, index1 = read_data.find_nearest(potential_shapes, potential_shapes[j])
        nothing, index2 = read_data.find_nearest(dts,              dt_close)
        ax2_eV.plot(base_potentials, NumericalHeating[:, index1, index2],
                    colors_and_symbols.symb_col(c),
                    label = r"$\Delta t$ = " + str(dts[index2]) + " as (" + potential_shapes[j] + ")")
        # When the NumericalHeating is negative, we have cooling, which would not appear on the log scale plot.
        cooling_indices = np.where(NumericalHeating[:, index1, index2] <= 0.0)
        if (len(cooling_indices[0]) >= 1):
            ax2_eV.plot(base_potentials[cooling_indices], abs(NumericalHeating[:, index1, index2][cooling_indices]),
                        '.' + colors_and_symbols.symb_col(c))
        c += 1


ax1_eV.set_xlabel(r"$\Delta t$ [as]")
ax1_eV.set_xlim((dts[0], dts[-1]))
plt.setp(ax1_Eh.get_xticklabels(), visible=False)

ax2_eV.set_xlabel("Potential depth [Hartree]")
ax2_eV.set_xlim((base_potentials[0], base_potentials[-1]))

for ax in [ax1_eV, ax2_eV]:
    leg = ax.legend(loc = "best")
    leg.get_frame().set_alpha(0.4)
    ax.grid(True)
    ax.set_ylabel("Energy change [eV]")
    ax.set_xscale('log')
    ax.set_yscale('log')

for ax in [ax1_Eh, ax2_Eh]:
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylabel("Energy change [Hartree]")

savefigure.show(figs, figure_name = "numerical_heating")
