#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import matplotlib
import matplotlib.pyplot
import math
import sys
import signal
import os

from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
import matplotlib.transforms as mtransforms

matplotlib.rcParams['figure.subplot.hspace'] = 0.6

def my_quit(signum, frame):
    sys.exit(0)

on_key_lines = []

def get_axis_two_scales(fig, scale_x, scale_y, \
                        ax2_xlabel = None, ax2_ylabel = None, \
                        subplot = 111,
                        sharex = None,
                        sharey = None):
    kargs = {}
    if (sharex != None):
        kargs['sharex'] = sharex
    if (sharey != None):
        kargs['sharey'] = sharey
    ax1 = SubplotHost(fig, subplot, **kargs)
    ax1_to_2 = mtransforms.Affine2D().scale(1.0/scale_x, 1.0/scale_y)
    ax2 = ax1.twin(ax1_to_2)
    ax2.set_viewlim_mode("transform")
    fig.add_subplot(ax1)
    if (ax2_xlabel != None):
        ax2.set_xlabel(ax2_xlabel)
    if (ax2_ylabel != None):
        ax2.set_ylabel(ax2_ylabel)
    if (scale_x == 1.0):
        ax2.get_xaxis().set_visible(False)
    if (scale_y == 1.0):
        ax2.get_yaxis().set_visible(False)
    return ax1, ax2

###################################################################
def on_key(event):
    global on_key_lines
    #print 'you pressed', event.key, event.xdata, event.ydata
    if   (event.key == 'q'):
        my_quit(0,0)
    #
    elif (event.key == 'w'):
        matplotlib.pyplot.close(event.canvas.figure)
    #
    elif (event.key == 'd'):
        print "##############################################################"
        print "Please click two points to get the distance and slope."
        from matplotlib.widgets import Cursor
        cursor = Cursor(event.inaxes, useblit=False, color='red', linewidth=1)
        points = matplotlib.pyplot.ginput(n=2, show_clicks=True, timeout=0)
        xs  = [points[0][0], points[1][0]]
        ys  = [points[0][1], points[1][1]]
        dx = xs[1] - xs[0]
        dy = ys[1] - ys[0]
        try:
            dy_log = math.log10(ys[0]) - math.log10(ys[1])
        except:
            dy_log = float('nan')
        try:
            dy_ln  = math.log(ys[0])   - math.log(ys[1])
        except:
            dy_ln  = float('nan')
        if (dx == 0.0):
            slope  = float('inf')
        else:
            slope  = dy/dx
        if (dy == 0.0):
            inv_slope = float('inf')
        else:
            inv_slope = dx/dy

        angle = math.atan2(dy, dx)
        print "points = ", points
        print "distance (x) =", dx
        print "distance (y) =", dy
        print "distance =", math.sqrt( dx**2 + dy**2 )
        print "Ratio: x0/x1 =", xs[0] / xs[1], "   x1/x0 =", xs[1] / xs[0], "  y0/y1 =", ys[0] / ys[1], "  y1/y0 =", ys[1] / ys[0]
        print "dy/y0 = ", dy/ys[0], "(", dy/ys[0]*100, "%)"
        print "dx/x0 = ", dx/xs[0], "(", dx/xs[0]*100, "%)"
        print "Angle: theta = atan2(y/x) =", angle, "rad =", angle*180.0/math.pi, "deg"
        print "Slope: ", slope, "  1/slope =", inv_slope
        print "Slope (log10 scale):", dy_log/dx
        print "Slope (ln scale):",    dy_ln/dx
        on_key_lines.append(event.inaxes.plot([points[0][0], points[1][0]], [points[0][1], points[1][1]], '--r', lw=1.0))
        on_key_lines.append(event.inaxes.plot([points[0][0], points[1][0]], [points[0][1], points[1][1]],  '+r', lw=1.0))
        matplotlib.pyplot.draw()
    #
    elif (event.key == 'a'):
        print "##############################################################"
        print "Please click a point to get the position."
        cursor = matplotlib.widgets.Cursor(event.inaxes, useblit=False, color='red', linewidth=1)
        point  = matplotlib.pyplot.ginput(n=1, show_clicks=False, timeout=0)
        print "    Point clicked:", point[0]
        on_key_lines.append(event.inaxes.plot(point[0][0], point[0][1],  '+r', lw=1.0))
        matplotlib.pyplot.draw()
    elif (event.key == 'x'):
        # Clear all annotations
        for line in on_key_lines:
           l = line.pop(0)
           l.remove()
           del line
        del on_key_lines
        on_key_lines = []
        matplotlib.pyplot.draw()
    elif (event.key == 'c'):
        # Switch current axis by changing their Z-order
        fig         = matplotlib.pyplot.gcf()
        current_ax  = matplotlib.pyplot.gca()
        axes        = [x for x in fig.canvas.figure.get_axes()]
        z_orders    = [x.get_zorder() for x in axes]
        nb_axes     = len(axes)
        # Find the current axis
        ci = 0
        for ax in axes:
            if (ax == current_ax):
                break
            ci += 1
        next_ci = (ci+1) % nb_axes

        # Compare z-orders and detect if they are all the same
        all_equal = True
        for i in xrange(nb_axes):
            if (z_orders[i] != z_orders[0]):
                all_equal = False
                break
        # If z-orders are all equal, set the current one on top
        if (all_equal):
            # Set a range around current values
            # Orders will start at -range/2 up to +range/2
            zorder_range = 0.1
            # Force the current one on top
            axes[ci].set_zorder(z_orders[ci] + zorder_range/2.0)
            # And set all others to under it
            k = 0
            for j in xrange(nb_axes):
                if (j != ci):
                    axes[j].set_zorder(z_orders[j] - zorder_range/2.0 + float(k)*zorder_range/float(nb_axes))
                    k += 1
            z_orders = [x.get_zorder() for x in axes]
        sorted_indices = [i[0] for i in sorted(enumerate(z_orders), key=lambda x:x[1])]

        # We push z-orders up the ladder, and send the first one in the sorted list (the smallest value)
        # to the end.
        # Store the smallest z-order
        smallest_zorder = z_orders[sorted_indices[0]]
        for i in xrange(nb_axes-1):
            j = sorted_indices[i]
            z_orders[sorted_indices[i]] = z_orders[sorted_indices[i+1]]
        # Set the smallest z-order to the last in the list
        z_orders[sorted_indices[-1]] = smallest_zorder

        # Now set the z-order of the axes
        for i in xrange(nb_axes):
            j = sorted_indices[i]
            axes[j].set_zorder(z_orders[j])

        matplotlib.pyplot.sca(axes[next_ci])
        matplotlib.pyplot.draw()

def on_pick(event):
    thisline = event.artist
    if (type(thisline) != matplotlib.legend.Legend):
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        print "Clicked on line \"" + thisline.get_label() + '\" on point',
        try:
            print (xdata[ind][0], ydata[ind][0])
        except IndexError:
            print "<no data>"

###################################################################

def figure(**kwargs):
    fig = matplotlib.pyplot.figure(**kwargs)
    fig.canvas.mpl_connect('key_press_event', on_key)
    fig.canvas.mpl_connect('pick_event',      on_pick)
    # Stupid Qt4Agg can't quit when ctrl+c (SIGINT) is send. Force it to quit.
    signal.signal(signal.SIGINT, my_quit)
    return fig

def setp(*args, **kwargs):
    matplotlib.pyplot.setp(*args, **kwargs)
def savefig(filenames):
    figures=[manager.canvas.figure
             for manager in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
    if (type(filenames) == str):
        filenames = [filenames]
    for i, figure in enumerate(figures):
        if (len(filenames) > 1):
            assert(len(filenames) == len(figures))
            filename = filenames[i]
        else:
            filename_base, filename_ext = os.path.splitext(filenames[0])
            filename = '{0}_{1}{2}'.format(filename_base, i, filename_ext)
        print "Saving to", filename
        figure.savefig(filename, transparent=True, bbox_inches='tight')


def show():
    fig         = matplotlib.pyplot.gcf()
    axes        = [x for x in fig.canvas.figure.get_axes()]
    # Set legends to be draggable
    for ax in axes:
        legend = ax.get_legend()
        if (ax.get_legend() != None):
            legend.draggable()
    # For every axes, get all lines plotted and set its picker value to 5 points
    for ax in axes:
        lines = ax.get_lines()
        for line in lines:
            line.set_picker(5)
    matplotlib.pyplot.show()

