import scipy.stats
import numpy as np
from verkko.binner import bins
from matplotlib import pyplot as plt
from matplotlib import cm, colors
from matplotlib.ticker import LogFormatterMathtext


def plot_counts(x,
                y,
                cmap=cm.get_cmap('summer'),
                xscale='log',
                yscale='log',
                xParam=1.5,
                yParam=1.5,
                ax=None,
                use_gridspec=False
                ):
    """
    Plots the counts as a heatmap.

    Parameters
    ----------
    x : list-like (1D numpy array)
        values on x-axis
    y : list-like (1D numpy array)
        values on x-axis
        if xbins and ybins are both given
    cmap : matplotlib.cm
        colormap to use for the plot, defaulting to summer
    xscale : {"log", "linear", "linlog", ...}, optional
        see binner.binner or binner.bins for more info on the options
    yscale : {"log", "linear", "linlog", ...}, optional
        see binner.binner or binner.bins for more info on the options
    xParam : varies according to xscale
        if xscale == 'log'
            xParam equals to the multiplier
        if xscale == 'lin'
            xParam equals to the number of bins
    yParam : varies according to yscale
        see xParam
    ax : axis object
        give this if plotting is used as a subroutine for existing
        axis. Created if None.
    use_gridspec : bool
         set true if subroutine plotting to a subplot with gridspec
         fixes the colorbar to the correct subplot
    """

    fig, ax, cax = _get_fig_ax_and_colorbar_ax(xscale, yscale, ax=ax)

    X, Y, counts_matrix, bin_centers, means, Bins2D = \
        _get_count_data_to_plot(x, y, xscale, yscale, xParam, yParam)

    im = ax.pcolor(X, Y, counts_matrix,
                   cmap=cmap, norm=colors.LogNorm()
                   )

    ax.plot(bin_centers, means, "go-")

    cbar = fig.colorbar(im, cax, ax=ax, use_gridspec=use_gridspec,
                        orientation='vertical',
                        format=LogFormatterMathtext())
    return fig, ax, cbar, im


def plot_prob_density(x,
                      y,
                      cmap=cm.get_cmap('summer'),
                      xscale='log',
                      yscale='log',
                      xParam=1.5,
                      yParam=1.5,
                      ax=None,
                      use_gridspec=False
                      ):
    """
    Plots the normalized probability density.

    See :py:func:`plot_counts` for explanations of the input variables.
    """
    fig, ax, cax = _get_fig_ax_and_colorbar_ax(xscale, yscale, ax=ax)

    X, Y, counts_matrix, bin_centers, means, Bins2D = \
        _get_count_data_to_plot(x, y, xscale, yscale, xParam, yParam)

    prob_density = counts_matrix / (Bins2D.sizes.T * np.sum(counts_matrix))

    assert np.abs(np.sum(prob_density * Bins2D.sizes.T) - 1) < 0.000001
    im = ax.pcolor(X, Y, prob_density,
                   cmap=cmap, norm=colors.LogNorm()
                   )

    ax.plot(bin_centers, means, "go-")

    cbar = fig.colorbar(im, cax, ax=ax, use_gridspec=use_gridspec,
                        orientation='vertical',
                        format=LogFormatterMathtext())
    return fig, ax, cbar, im


def plot_conditional_prob_density(x,
                                  y,
                                  cmap=cm.get_cmap('summer'),
                                  xscale='log',
                                  yscale='log',
                                  xParam=np.sqrt(2),
                                  yParam=np.sqrt(2),
                                  ax=None,
                                  use_gridspec=False
                                  ):
    """
    Plots the conditional probability density. (P(y|x)).

    See plot_counts for explanations of the input variables.
    """
    fig, ax, cax = _get_fig_ax_and_colorbar_ax(xscale, yscale, ax=ax)

    X, Y, counts_matrix, bin_centers, means, Bins2D = \
        _get_count_data_to_plot(x, y, xscale, yscale, xParam, yParam)

    # count vals for each x-bin
    x_counts = np.sum(counts_matrix, axis=0)

    # normalize the counts for each x-bin by the total number
    norm_x_counts_mat = counts_matrix / x_counts[None, :]

    # to get the prob density, normalize by the bin widths
    cond_prob = norm_x_counts_mat / Bins2D.widths[1][:, None]

    im = ax.pcolor(X, Y, cond_prob,
                   cmap=cmap, norm=colors.LogNorm()
                   )

    ax.plot(bin_centers, means, "go-")

    cbar = fig.colorbar(im, cax, ax=ax, use_gridspec=use_gridspec,
                        orientation='vertical',
                        format=LogFormatterMathtext())
    return fig, ax, cbar, im


def _get_fig_ax_and_colorbar_ax(xscale, yscale, ax=None):
    # create figure and axis if not given
    if ax is None:
        fig = plt.figure()
        y_low = 0.07
        y_height = 0.85
        ax = fig.add_axes([0.1, y_low, 0.8, y_height])
        cax = fig.add_axes([0.92, y_low, 0.03, y_height])
    else:  # set colorbar axes to None so the
            # space is stolen from the given plotting axes
        cax = None
        fig = ax.get_figure()

    # set the scale
    if "log" in xscale:
        ax.set_xscale('log')
    if "log" in yscale:
        ax.set_yscale('log')
    return fig, ax, cax


def _get_count_data_to_plot(x, y, xscale, yscale, xParam, yParam):
    """
    See e.g. function plot_counts for interpreting the inner workings
    of this function.
    """

    if type(x) is not np.ndarray:
        x = np.array(x)
    if type(y) is not np.ndarray:
        y = np.array(y)

    xidx, yidx = bins.get_reasonable_data_indices_for_binning(
        x, y, xscale, yscale)
    min_x, max_x = min(x[xidx]), max(x[xidx])
    min_y, max_y = min(y[yidx]), max(y[yidx])

    bins2D = bins.Bins2D(
        float, min_x, max_x, xscale, xParam,
        float, min_y, max_y, yscale, yParam
    )

    counts_matrix, _, _ = np.histogram2d(
        x, y, bins=bins2D.bin_limits)

    # counts_matrix = counts_matrix.astype(np.float64)
    counts_matrix = counts_matrix.T
    counts_matrix = np.ma.masked_array(counts_matrix, counts_matrix == 0)

    X, Y = bins2D.edge_grids
    # compute bin average as a function of x:
    x_bin_means, _, _ = scipy.stats.binned_statistic(
        x, x, statistic='mean', bins=bins2D.xbin_lims
    )
    y_bin_means, _, _ = scipy.stats.binned_statistic(
        x, y, statistic='mean', bins=bins2D.xbin_lims
    )
    return X, Y, counts_matrix, x_bin_means, y_bin_means, bins2D
