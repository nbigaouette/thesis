#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import sys, os, glob, math
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
import matplotlib.transforms as mtransforms

import netCDF4

import constants as cst
import plot
import colors_and_symbols

import plot_params

max_plot = 7
saved_file    = "numerical_heating.npz"

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

print "base_potentials =", base_potentials
print "potential_shapes =", potential_shapes
print "dts =", dts

fig1 = plot.figure()
fig2 = plot.figure()
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
dbase_potentials = (base_potentials[-1] - base_potentials[0]) / float(max_plot-1) # -1 since we want the number of intervals
c = 0
for j in xrange(len(potential_shapes)):
    for base_potentials_close in np.arange(base_potentials[0], base_potentials[-1]+dbase_potentials/2.0, dbase_potentials):
        nothing, index0 = find_nearest(base_potentials,  base_potentials_close)
        nothing, index1 = find_nearest(potential_shapes, potential_shapes[j])
        #nothing, index2 = find_nearest(dts,              dts[i])
        ax1_eV.plot(dts, NumericalHeating[index0, index1, :],
                    colors_and_symbols.symb_col(c),
                    label = r"" + str(base_potentials[index0]) + " Eh")
                    #label = r"Potential depth: " + str(base_potentials[index0]) + " Eh (" + potential_shapes[j] + " )")
        # When the NumericalHeating is negative, we have cooling, which would not appear on the log scale plot.
        cooling_indices = np.where(NumericalHeating[index0, index1, :] <= 0.0)
        if (len(cooling_indices[0]) >= 1):
            ax1_eV.plot(dts[cooling_indices], abs(NumericalHeating[index0, index1, :][cooling_indices]),
                        '.' + colors_and_symbols.symb_col(c))
        c += 1

# Plot NumericalHeating as a function of potential depth for every dt
ddt = (dts[-1] - dts[0]) / float(max_plot-1+1) # -1 since we want the number of intervals
                                               # +1 so the skipped dt = 0.005 as does not reduce the number of curves
c = 0
for j in xrange(len(potential_shapes)):
    for dt_close in np.arange(dts[0], dts[-1]+ddt/2.0, ddt):
        nothing, index1 = find_nearest(potential_shapes, potential_shapes[j])
        nothing, index2 = find_nearest(dts,              dt_close)
        if (dts[index2] <= 0.005):
           continue
        ax2_eV.plot(base_potentials, NumericalHeating[:, index1, index2],
                    colors_and_symbols.symb_col(c),
                    label = r"" + str(dts[index2]) + " as")
                    #label = r"$\Delta t$ = " + str(dts[index2]) + " as (" + potential_shapes[j] + ")")
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

plot.savefig('numerical_heating.svg')
plot.savefig('numerical_heating.pdf')
plot.show()
