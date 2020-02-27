import numpy as np
import matplotlib.pyplot as plt

from verkko.binner import bins


def plot_pdf(values, ax=None, xscale='lin', xParam=None, yscale='lin'):
    """
    Plots the probability density function of given values.

    Parameters
    ----------

    values : numpy ndarray
        the values for which the experimental pdf is computed
    ax : matplotlib axes object, optional
        axes to plot the figure in
    xscale : str
        'lin' or 'log', or ... (see binner.Bins for more details)
    yscale : str
        'lin' or 'log'
    xParam : different things, optional
        see binner.Bins for more details

    Returns
    -------
    fig : matplotlib Figure
        the parent figure of the axes
    ax : matplotlib Axes object
        the axes in which the pdf is plotted
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()

    indices = bins.get_reasonable_data_indices_for_binning(
        values, xscale=xscale)

    prop_vals = values[indices]

    if xParam is None:
        if xscale == 'log':
            xParam = np.sqrt(2)  # just getting a nice factor of two...
        if xscale == 'lin':
            xParam = 50
    xbins = bins.Bins(float, np.min(prop_vals),
                      np.max(prop_vals), xscale, xParam)
    ax.hist(values, bins=xbins.bin_limits, normed=True)

    if 'log' in xscale:
        ax.set_xscale('log')
    if 'log' in yscale:
        ax.set_yscale('log')

    return fig, ax


def plot_ccdf(values, ax=None, xscale='lin', xParam=None, yscale='log',
              threshold_data=False):
    """
    Plot the experimental 1-CDF of values.

    Parameters
    ----------
    See :py:func:`plot_pdf` for explanations of the parameters.

    bin_data : bool
        whether to use thresholds for drawing the plot
        (more efficient drawing if a lot of points present)

    Returns
    -------
    fig : matplotlib Figure
        the parent figure of the axes
    ax : matplotlib Axes object
        the axes in which the ccdf is plotted
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()

    if threshold_data:
        print "no thresholding implemented yet in plotting" + \
            " 1 - CDF, defaulting to basic stuff"

    xvals = np.sort(values)
    yvals = np.linspace(1, 1 / len(xvals), len(xvals))
    ax.plot(xvals, yvals)
    if 'log' in xscale:
        ax.set_xscale('log')
    if 'log' in yscale:
        ax.set_yscale('log')
    return fig, ax
