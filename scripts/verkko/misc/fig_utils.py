# Different commands for making figures.

import sys
import os
import pylab
import matplotlib.ticker as ticker
import numpy as np
#import binner_noCoroutines as binner
from verkko.binner import binner

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def savefig(fig, name, extensions=None, verbose=False):
    """Save figure.

    Save matplotlib.figure object `fig` as `name`.EXT, where EXT are
    given in `extensions`. If only one save type is used, the full
    name (including the extension) can also be given as `name`.

    Note! Saving as 'svg' requires that the program 'pdf2svg' is
    installed on your system.
    """

    if len(name) == 0:
        raise ValueError("File name can not be empty.")

    if extensions is None:
        fields = name.split(".")
        if len(fields) == 1:
            raise ValueError("File name must contain an extension if"
                             " extensions are not given explicitely.")
        extensions = fields[-1]
        name = ".".join(fields[:-1])

    if isinstance(extensions, str):
        extensions = (extensions,)

    # Check if pdf should be generated (both eps and svg will be
    # created from pdf) and generate if necessary.
    pdf_generated = False
    pdf_tmp = "%s_tmp_%d.pdf" % (name, np.random.randint(100000))
    if set(['pdf', 'svg', 'eps']).intersection(extensions):
        fig.savefig(pdf_tmp)
        pdf_generated = True

    for ext in extensions:
        if not isinstance(ext, str):
            raise ValueError("'extensions' must be a list of strings.")
        if ext[0] == '.':
            ext = ext[1:]

        if ext == 'eps':
            pipe = os.popen("pdftops -eps %s %s.eps" % (pdf_tmp, name))
            exit_status = pipe.close()
            if exit_status:
                if os.WEXITSTATUS(exit_status) == 127:
                    sys.stderr.write("%s could not be created because program "
                                     "'pdftoeps' could not be found.\n" % (ext,))
                else:
                    sys.stderr.write("Problem saving '%s'.\n" % (ext,))
        elif ext == 'svg':
            pipe = os.popen("pdf2svg %s %s.svg" % (pdf_tmp, name))
            exit_status = pipe.close()
            if exit_status:
                if os.WEXITSTATUS(exit_status) == 127:
                    sys.stderr.write("%s could not be created because program "
                                     "'pdf2svg' could not be found.\n" % (ext,))
                else:
                    sys.stderr.write("Problem saving '%s'.\n" % (ext,))
        elif ext != 'pdf':
            # fig.savefig raises a ValueError if the extension is not identified.
            #
            # According to "http://stackoverflow.com/questions/4581504/how-to-set-
            # opacity-of-background-colour-of-graph-wit-matplotlib" it is necessary
            # to set the face- and edgecolor again!.
            fig.savefig(name + "." + ext, dpi=300)

    if pdf_generated:
        if 'pdf' in extensions:
            os.popen("mv %s %s.pdf" % (pdf_tmp, name))
        else:
            os.popen("rm %s" % (pdf_tmp,))


def get_rcParams(fig_width_cm, fig_ratio=0.8, font_sizes=None):
    """Set good parameters for LaTeX-figures.

    The idea idea is to set the figure width in centimeters to be the
    same as the final size in your LaTeX document. This way the font
    sizes will be correct also.

    Parameters
    ----------
    fig_width_cm: int or float
        The width of the final figure in centimeters.
    fig_ratio: float (between 0 and 1)
        The ratio height/width. < 1 is landscape, 1.0 is square and
        > 1.0 is portrait.
    font_sizes: dictionary
        The font sizes used in the figure. Default is size 10 for the
        title and 8 for everything else. Possible keys are 'default',
        'label', 'title', 'text', 'legend' and 'tick'.  'default' is
        used when the specific value is not defined, other keys should
        be self explanatory.
    """
    default_font_sizes = {
        'label': 8, 'title': 10, 'text': 8, 'legend': 8, 'tick': 8}
    font_sizes = (font_sizes or {})
    for k in default_font_sizes:
        if k not in font_sizes:
            font_sizes[k] = (
                font_sizes.get('default') or default_font_sizes[k])

    inches_per_cm = 1 / 2.54
    fig_width = 1.0 * fig_width_cm * inches_per_cm  # width in inches
    fig_height = 1.0 * fig_width * fig_ratio        # height in inches
    fig_size = [fig_width, fig_height]
    params = {'font.family': 'serif',
              'font.serif': 'Computer Modern Roman',
              'axes.labelsize': font_sizes['label'],
              'axes.titlesize': font_sizes['title'],
              'text.fontsize': font_sizes['text'],
              'font.size': font_sizes['text'],
              'legend.fontsize': font_sizes['legend'],
              'xtick.labelsize': font_sizes['tick'],
              'ytick.labelsize': font_sizes['tick'],
              'text.usetex': True,
              'figure.figsize': fig_size,
              'legend.labelspacing': 0.0,
              'lines.markersize': 3,
              'lines.linewidth': 0.5}
    return params


class EvenExpFormatter(ticker.LogFormatterMathtext):

    """Print labels only for even exponentials. Exponents given in
    'exclude' will also be skipped.
    """

    def __init__(self, base=10.0, labelOnlyBase=True, exclude=None):
        if exclude == None:
            self.exclude = []
        else:
            self.exclude = exclude
        ticker.LogFormatterMathtext.__init__(self, base, labelOnlyBase)

    def __call__(self, val, pos=None):
        fx = int(np.floor(np.log(abs(val)) / np.log(self._base) + 0.5))
        isDecade = self.is_decade(fx)
        if not isDecade and self.labelOnlyBase:
            return ''
        if (fx % 2) == 1 or (fx in self.exclude):  # odd, skip
            return ''

        return ticker.LogFormatterMathtext.__call__(self, val, pos)


class SkipMathFormatter(ticker.LogFormatterMathtext):

    """Skip exponents given in 'exclude'.
    """

    def __init__(self, base=10.0, labelOnlyBase=True, exclude=None):
        if exclude == None:
            self.exclude = []
        else:
            self.exclude = exclude
        ticker.LogFormatterMathtext.__init__(self, base, labelOnlyBase)

    def __call__(self, val, pos=None):
        fx = int(np.floor(np.log(abs(val)) / np.log(self._base) + 0.5))
        isDecade = self.is_decade(fx)
        if not isDecade and self.labelOnlyBase:
            return ''
        if fx in self.exclude:  # skip
            return ''

        return ticker.LogFormatterMathtext.__call__(self, val, pos)


class SelectiveScalarFormatter(ticker.ScalarFormatter):

    """Print only the ticks in the input list.
    """

    def __init__(self, printList=None, useOffset=True, useMathText=False):
        if printList is None:
            printList = []
        self.printList = printList
        ticker.ScalarFormatter.__init__(self, useOffset, useMathText)

    def __call__(self, val, pos=None):
        if val in self.printList:
            return ticker.ScalarFormatter.__call__(self, val, pos)
        else:
            return ''


class DivisorFormatter(ticker.FormatStrFormatter):

    """Divide all numbers by a constant."""

    def __init__(self, fmt, divisor=None):
        if divisor == None:
            self.divisor = 1
        else:
            self.divisor = divisor
        ticker.FormatStrFormatter.__init__(self, fmt)

    def __call__(self, val, pos=None):
        return ticker.FormatStrFormatter.__call__(self, int(1.0 * val / self.divisor), pos)


def distribution_2d(ax, bins, data, log_scaling=True, plot_average=True):
    """Plot a 2-dimensional distribution P(Y|X)

    Parameters
    ----------
    ax : pylab.axes
       The axis object where the distribution is drawn.
    bins : binner.Bins2D
       The 2-D bins to use for binning.
    data : (x,y)-points, iterable
       The data points to plot. If `plot_average` is true, this must
       be a sequence (iterable is not enough).
    log_scaling : bool
       If true, the colors that denote probabilities will be scaled
       logarithmically. (This is usually a good idea.)
    plot_average : bool
       If true, the average value will be drawn with black line (if
       x-values are floats) or with black dots (if x-values are
       integers).
    """

    # Get binned data for the distribution
    binned_data = bins.bin_count_divide(data)

    # Create coordinates and create an array from the data
    X, Y = bins.edge_grids
    Z = np.array(binned_data, float)
    Z = np.ma.masked_array(Z, Z == 0)

    # Normalize each weight bin separately.
    for j in range(len(Z)):
        row_sum = np.sum(Z[j])
        if row_sum:
            Z[j] = Z[j] / np.sum(Z[j])
    Z = np.transpose(Z)

    # Plot
    if log_scaling:
        norm = pylab.matplotlib.colors.LogNorm()
    else:
        norm = pylab.matplotlib.colors.Normalize()
    im = ax.pcolor(X, Y, Z, norm=norm)

    # Set scaling based on the bins used.
    if isinstance(bins.x_bin_finder, binner._logBinFinder):
        ax.set_xscale("log")
    if isinstance(bins.y_bin_finder, binner._logBinFinder):
        ax.set_yscale("log")

    # Set axis limits to correspond to bin limits
    ax.axis((bins.bin_limits[0][0], bins.bin_limits[0][-1],
             bins.bin_limits[1][0], bins.bin_limits[1][-1]))

    if plot_average:
        # Create bins for the average. The first line creates a dummy
        # binner; the intestines are then changed to those in `bins`
        # on the next two lines.
        avg_bins = binner.Bins(int, 0, 10, 'lin', 1)  # Dummy
        avg_bins.bin_limits = bins.bin_limits[0]
        avg_bins.bin_finder = bins.x_bin_finder
        binned_avg = avg_bins.bin_average(data)
        if avg_bins.bin_limits.dataType == 'int':
            plot_type = 'ko'
        else:
            plot_type = 'k-,'
        ax.plot(avg_bins.centers, binned_avg, plot_type, lw=1.0)

    return im


def scatter_2d(ax, bins, data, log_scaling=True, draw_diagonal=False):
    """Plot a binned scatter plot.

    Parameters
    ----------
    ax : pylab.axes
       The axis object where the distribution is drawn.
    bins : binner.Bins2D
       The 2-D bins to use for binning.
    data : (x,y)-points, iterable
       The data points to plot.
    log_scaling : bool
       If true, the colors that denote probabilities will be scaled
       logarithmically. (This is usually a good idea.)
    draw_diagonal : bool
       If true, the diagonal line is draws on top of the plot.
    """

    # Get binned data for the distribution
    binned_data = bins.bin_count_divide(data)

    # Create coordinates and create an array from the data
    X, Y = bins.edge_grids
    Z = np.array(binned_data, float)
    Z = np.ma.masked_array(Z, Z == 0)

    # Normalize the data as a whole.
    Z = Z / Z.sum()
    Z = np.transpose(Z)

    # Plot
    if log_scaling:
        norm = pylab.matplotlib.colors.LogNorm()
    else:
        norm = pylab.matplotlib.colors.Normalize()
    im = ax.pcolor(X, Y, Z, norm=norm)

    # Set scaling based on the bins used.
    if isinstance(bins.x_bin_finder, binner._logBinFinder):
        ax.set_xscale("log")
    if isinstance(bins.y_bin_finder, binner._logBinFinder):
        ax.set_yscale("log")

    # Set axis limits to correspond to bin limits
    ax.axis((bins.bin_limits[0][0], bins.bin_limits[0][-1],
             bins.bin_limits[1][0], bins.bin_limits[1][-1]))

    if draw_diagonal:
        v_ax = ax.axis()
        ax.plot(v_ax[0:2], v_ax[2:], 'k-', lw=1.0)

    return im


def pretty_number(x):
    """Turn integer into a pretty string.

    Make a pretty number by adding extra spaces between every three
    digits. The LaTeX environment must be active.

    Parameters
    ----------
    x : int
       The integer to print.

    Returns
    -------
    s : str
       String representation of the integer.
    """

    if isinstance(x, float):
        return r"$%.2f$" % x

    s_tmp = str(abs(x))
    s = ""
    while s_tmp:
        if len(s_tmp) > 3:
            s = r"\," + s_tmp[-3:] + s
            s_tmp = s_tmp[:-3]
        else:
            s = ("-" if x < 0 else "") + s_tmp + s
            break
    return r"$%s$" % s


def save_all_open_figs_to_pdf(fname):
    with PdfPages(fname) as pdf:
        fignums = plt.get_fignums()
        for i in fignums:
            fig = plt.figure(i)
            pdf.savefig(fig)
