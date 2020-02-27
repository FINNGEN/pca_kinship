from __future__ import absolute_import

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from verkko.misc import fig_utils as fu  # relative import

# This file contains small helper function to be used alongwith matplotlib


def setFigure(format='PNAS', fig_size=None, fig_ratio=None, ):
    '''
    Get a basic figure object
    Input:
        format : 'PNAS', 'PRE', 'Nature'

    '''
    if fig_ratio is None:
        fig_ratio = 1.0

    if format == 'PNAS' or 'Nature':
        font_size = {'label': 6, 'title': 6, 'text': 6, 'legend': 6, 'tick': 6}
        if fig_size is None:
            fig_size = 8.7
    if format == 'PRE':
        font_size = {'label': 7, 'title': 7, 'text': 7, 'legend': 7, 'tick': 7}
        if fig_size is None:
            fig_size = 8.6

    fig = plt.figure()
    params = fu.get_rcParams(fig_size, fig_ratio=fig_ratio, font_sizes=font_size)
    plt.rcParams.update(params)
    return fig


def setFigureGrid(fig, grid_geometry=(2, 2)):
    '''
    Set up the grid
    '''
    gs = gridspec.GridSpec(*grid_geometry)
    axarr = []
    for i in range(grid_geometry[0]):
        for j in range(grid_geometry[1]):
            axarr.append(fig.add_subplot(gs[i, j]))
    return gs, axarr


def setLegend(ax, lines=None, legends=None, loc='upper left'):
    '''
    Get the legend
    Input:
        ax: axes
        (optional arguments)
        lines: Lines
        legends: Legends
        loc: Location (default is upper left)
    Output:
        leg: The legend handle
    '''
    if lines is None and legends is None:
        lines, legends = ax.get_legend_handles_labels()
    elif lines is None:
        lines = ax.get_lines()
    elif legends is None:
        legends = ax.get_legend()

    if lines is None or legends is None:
        print("No labeled objects found.")
        return None
    else:
        leg = ax.legend(lines, legends, loc=loc, numpoints=1, fancybox=True, handletextpad=0.2, borderpad=0.15, handlelength=2.0, columnspacing=0.5, )
        leg.get_frame().set_alpha(0.5)
        leg.get_frame().set_lw(0.2)
        return leg


def setTickLabelsInvisible(axarr, whichaxis='x'):
    '''
    Make the x of y ticks invisible

    Parameters
    ----------
    ax : matplotlib axes obj
        axes or list of axes
    whichaxis : str
        x, y or xy
    '''
    if not isinstance(axarr, list):
        axarr = [axarr]
    if whichaxis == 'x' or whichaxis == 'xy':
        for ax in axarr:
            for label in ax.get_xticklabels():
                label.set_visible(False)
    if whichaxis == 'y' or whichaxis == 'xy':
        for ax in axarr:
            for label in ax.get_yticklabels():
                label.set_visible(False)


def setAxisLabels(axarr, labels, whichaxis='x'):
    '''
    Make the x of y labels
    Input:
        axarr: axes or list of axes
        labels : labels
        (optional arguments)
        ticklabels: x or y
    '''
    if not isinstance(axarr, list):
        axarr = [axarr]
        labels = [labels]
    if whichaxis == 'x':
        for ax, label in zip(axarr, labels):
            ax.set_xlabel(label)
    elif whichaxis == 'y':
        for ax, label in zip(axarr, labels):
            ax.set_ylabel(label)


def setAxisLimits(axarr, lims, whichaxis='x'):
    '''
    Set the x or y axis limits
    Input:
        axarr: axes or list of axes
        lims : limits
        (optional arguments)
        ticklabels: x or y
    '''
    if not isinstance(axarr, list):
        axarr = [axarr]
        lims = [lims]
    if whichaxis == 'x':
        for ax, lim in zip(axarr, lims):
            ax.set_xlim(*lim)
    elif whichaxis == 'y':
        for ax, lim in zip(axarr, lims):
            ax.set_ylim(*lim)


def setFigureIndex(axarr, letters=None, xy=(0, 1), xytext=(6, -2)):
    '''
    Put ABC to the figure index
    '''
    if letters is None:
        letters = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    if not isinstance(axarr, list):
        axarr = [axarr]

    for ax, letter in zip(axarr, letters):
        ax.annotate(letter, xy=xy, xytext=xytext, va='top', ha='center', xycoords='axes fraction', textcoords='offset points', size=8, weight='extra bold')
