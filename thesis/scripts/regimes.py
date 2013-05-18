#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
import plot
from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
import matplotlib.transforms as mtransforms
from matplotlib.patches import FancyBboxPatch

import plot_params

import constants as cst


def eV_to_mum(energy_eV):
    energy = energy_eV * cst.eV_to_J
    omega       = energy / cst.hbar
    f           = omega / (2.0 * cst.pi)
    k           = omega / cst.c
    wavelength  = 2.0 * cst.pi / k
    return wavelength * cst.m_to_mum
def eV_to_nm(energy_eV):
    return cst.mum_to_m * cst.m_to_nm * eV_to_mum(energy_eV)

def curve_from_loglog_pts(x0, x1, F0, F1, x):
    F = F0 * (x / x0)**( np.log(F1/F0) / np.log(x1/x0) )
    return F

def log_center(dmin, dmax):
    x = 10.0**((np.log10(dmax) + np.log10(dmin)) / 2.0)
    return x
def bbox_log_center(bb):
    return [log_center(bb.xmin, bb.xmax), log_center(bb.ymin, bb.ymax)]

def draw_text_box(ax, text, xmin, xmax, ymin, ymax, **kwargs):
    bb = mtransforms.Bbox([[xmin, ymin], [xmax, ymax]])
    ax.add_patch(FancyBboxPatch(xy          = (bb.xmin, bb.ymin),
                                width       = abs(bb.width),
                                height      = abs(bb.height),
                                boxstyle    = 'round',
                                **kwargs))
    ax.text(x = log_center(bb.xmin, bb.xmax),
            y = log_center(bb.ymin, bb.ymax),
            s = text,
            horizontalalignment = 'center',
            verticalalignment   = 'center')




eV_min = 0.1
eV_max = 1000.0
I_min  = 1.0e11
I_max  = 1.0e20
angle = 33


#fig = plot.figure(figsize = (12.0, 7.0))
fig = plot.figure()


ax_eV = fig.add_subplot(1,1,1)
ax_wl = ax_eV.twinx()
#ax_wl = ax_eV.twinx().twiny()

ax_eV.set_xlabel('Intensity [W/cm$^2$]')
ax_eV.set_ylabel('Photon energy [eV]')
#ax_wl.set_xlabel('Field strength [V/$\AA$]')
ax_wl.set_ylabel('Wavelength [nm]')

#for ax in [ax_eV]:
for ax in [ax_eV, ax_wl]:
    ax.set_yscale('log')
    ax.set_xscale('log')

I = 10.0**np.linspace(int(np.log10(I_min)), int(np.log10(I_max)), 10)

Up_1eV  = curve_from_loglog_pts(1e11, 6e18, 0.12, 1000, I)
Up_1meV = curve_from_loglog_pts(1e11,  7e15, 3.7, 1000, I)
Up_1keV = curve_from_loglog_pts(7e13,  1e20, 0.1, 120, I)


ax_eV.set_xlim((I_min, I_max))
ax_eV.set_ylim((eV_min, eV_max))
ax_wl.set_ylim((eV_to_nm(eV_min), eV_to_nm(eV_max)))




ax_eV.plot(I, Up_1eV,  '-k',  alpha = 0.5, label = '$U_p = 1 eV$')
ax_eV.plot(I, Up_1meV, '--r', alpha = 0.5, label = '$U_p = 1 meV$')
ax_eV.plot(I, Up_1keV, '--b', alpha = 0.5, label = '$U_p = 1 keV$')


ax_eV.annotate('$U_p = 1 meV$',
               fontsize = 16,
               xy=(4.8611e11, 8.44),
               xycoords='data',
               bbox=dict(boxstyle="round", fc="0.8"),
               rotation = angle, horizontalalignment = 'center', verticalalignment = 'center')
ax_eV.annotate('$U_p = 1 eV$',
               fontsize = 16,
               xy=(4.5e14, 8.2),
               xycoords='data',
               bbox=dict(boxstyle="round", fc="0.8"),
               rotation = angle, horizontalalignment = 'center', verticalalignment = 'center')
ax_eV.annotate('$U_p = 1 keV$',
               fontsize = 16,
               xy=(4.6e17, 8.11),
               xycoords='data',
               bbox=dict(boxstyle="round", fc="0.8"),
               rotation = angle, horizontalalignment = 'center', verticalalignment = 'center')

ax_eV.text(2.28e16, 124.52, 'Photon dominated regimes', rotation = angle,
           horizontalalignment = 'center', verticalalignment = 'center')
ax_eV.text(5.27e16, 40.7,   'Field dominated regimes',  rotation = angle,
           horizontalalignment = 'center', verticalalignment = 'center')



draw_text_box(ax_eV, "Relativistic regime",
              7.0e17, 4.0e19, 0.1, 1.5,
              facecolor  = (1., 1., 1.),
              edgecolor  = (1., 1., 1.))

draw_text_box(ax_eV, "Optical femtosecond lasers (IR)",
              1.0e11, 1.0e20, 1.5, 3.5,
              facecolor  = 'red',
              edgecolor  = 'red',
              alpha = 0.6)

draw_text_box(ax_eV, "Infrared FEL",
              1.0e11, 6.0e13, 0.1, 0.7,
              facecolor  = (1., .8, 1.),
              edgecolor  = (1., 0.5, 1.))

#draw_text_box(ax_eV, "VUV-FEL\nFLASH",
draw_text_box(ax_eV, "XUV\n \n \nVUV",
              1.0e11, 1.0e15, 11, 200,
              facecolor  = 'magenta',
              edgecolor  = 'magenta',
              alpha = 0.6)

draw_text_box(ax_eV, "XFEL",
              1.0e11, 1.0e17, 250, 1000,
              facecolor  = 'blue',
              edgecolor  = 'blue',
              alpha = 0.6)




plot.savefig('regimes.svg')
plot.savefig('regimes.pdf')
plot.show()

