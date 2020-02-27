#!/usr/bin/python2.6
import numpy as np
from math import ceil, floor, sqrt

# Based on the work of Lauri Kovanen, 2009-2011 (lauri.kovanen@gmail.com)
#
# Department of Biomedical Engineering and Computational Science
# Aalto University School of Science

"""Convenience classes describing 1- or 2-dimensional bins.

   Bins: D 1-dimensional bins.
   Bins2D: Bin data into 1-dimensional bins.

Example of giving data as input:
   bins = binner.Bins(int, 0, 200, 'lin', 10)

"""


class Bins(object):

    """ Class for describing 1-dimensional bins."""

    def __init__(self, dataType, minValue, maxValue, binType, param=None):
        """Initialize bins.

        Constructs the bins to be used in binning. You can select
        between several different types of bin limits, and if none of
        these seem fit, you can also supply your own bin limits.

        Parameters
        ----------
        dataType : type
            The type of the data, either int or float.  The data type
            affects both bin widths and cecnters. For example bin
            [1.1, 2.5] has width 1.4 and center 1.8 if `dataType`
            float, but width 1 and center 2 if `dataType` int.
        minValue : float or integer
            The minimum value that can be placed into the bins.
        maxValue : float or integer
            The maximum value that can be placed into the
            bins. maxValue must always be larger than minValue,
            otherwise BinLimitError is raised.
        binType : string
            The method used for constructing bin limits. More
            information below in section 'Bin types'.
        param : float, integer or list
            Additional parameter for controlling the construction of
            bin limits. The meaning of this parameter depends on bin
            type.

        Bin types
        ---------
        There are several methods for constructing the bin limits
        automatically for given parameters. The possible values for
        binType and the corresponding meaning of the extra parameter
        'param' are listed below.

        'lin' or 'linear': param = N_bin (integer, > 0)
            Construct N_bin equally sized bins minValue and maxValue, both
            inclusive. The first (last) bin is centered around minValue
            (maxValue).

        'log' or 'logarithmic': param = factor (float, > 1)
            Construct logarithmic bins between positive values minValue
            and maxValue, both inclusive, using factor to increase bin
            size. First bin will be centered at minValue and subsequent
            bins limits are calculated by multiplying previous limit by
            factor.

        'linlog': param = factor (float, > 1)
            Construct linear bins with width 1.0 from max(0, minValue)
            to min(10, maxValue) and logarithmic bins from thereon,
            using factor to increase bin size. Both minValue and
            maxValue must be integers whenever they are smaller than
            11, otherwise expection ValueError is raised.

        'maxlog': (no parameter)
            Construct as many logarithmic bins as possible between
            minValue-diff and maxValue+diff so that each bin contains
            at least one integer. minValue and maxValue must be
            integers.

        'linmaxlog': (no parameter)
            Same as 'maxlog', but the use linear bins for small values
            exactly as in 'linlog'.

        'custom': param = bin_limits (sequence)
            Use any strictly increasing sequence as bin limits. A
            value is put into bin i with
                bin_limits[i] <= value < bin_limits[i+1],
            with the exception of the last bin, which also includes
            bin_limits[-1]. In other words, both bin_limits[0] and
            bin_limits[-1] are inclusive. `minValue` and `maxValue`
            are ignored.
        """

        if not isinstance(dataType, type):
            raise ParameterError("dataType must be int or float.")
        self.bin_limits = _BinLimits(dataType, minValue, maxValue,
                                     binType, param)

    def __len__(self):
        """Return the number of bins."""
        return len(self.bin_limits) - 1

    # Create getter for bin centers.
    @property
    def centers(self):
        """Return bin centers as array."""
        try:
            return self._bin_centers
        except AttributeError:
            self._bin_centers = self.bin_limits.centers()
            return self._bin_centers

    # Create getter for bin widths.
    @property
    def widths(self):
        """Return bin widths as array."""
        try:
            return self._bin_widths
        except AttributeError:
            self._bin_widths = self.bin_limits.widths()
            return self._bin_widths

    def __check_data_element(self, elem, N):
        """Check one element of input data.

        Makes sure that elem is a sequence and elem[0] fits inside the
        bin limits. The required length N is _not_ checked (for
        performance reasons this is better done when actually needed)
        but is only shown in the error message if elem is not a sequence.
        """
        # Check bin limits and correct sequence type.
        try:
            if (elem[0] < self.bin_limits.minValue or
                    elem[0] > self.bin_limits.maxValue):

                raise BinLimitError("Value %g is not in the interval [%g, "
                                    "%g]." % (elem[0],
                                              self.bin_limits.minValue,
                                              self.bin_limits.maxValue))
        except (TypeError, IndexError):
            # TypeError occurs when data is a list and elem is
            # integer or a float. Rather surprisingly, numpy
            # raises an IndexError in the same situation; try for
            # instance creating a=numpy.array([1,2,3]) and then
            # call a[0][0].
            raise DataTypeError("Elements of input data must be sequences"
                                " with length at least %d." % (N,))


class Bins2D(object):

    """ Class for describing 2-dimensional bins."""

    def __init__(self, X_dataType, X_minValue, X_maxValue, X_binType, X_param,
                 Y_dataType, Y_minValue, Y_maxValue, Y_binType, Y_param
                 ):
        """Initialize bins.

        Constructs the bins to be used in binning. You can select
        between several different types of bin limits, and if none of
        these seem fit, you can also supply your own bin limits.

        Usage is identical to Bins, except that you have to supply
        necessary parameters for bins in both x- and y-directions. See
        the documentation of Bins for more information.
        """
        self.xbin_lims = _BinLimits(X_dataType, X_minValue, X_maxValue,
                                    X_binType, X_param)
        self.ybin_lims = _BinLimits(Y_dataType, Y_minValue, Y_maxValue,
                                    Y_binType, Y_param)
        self.bin_limits = [self.xbin_lims, self.ybin_lims]

    @property
    def shape(self):
        """Shape of bins."""
        try:
            return self._shape
        except AttributeError:
            """Return the number of bins in x- and y-directions."""
            self._shape = (
                len(self.bin_limits[0]) - 1, len(self.bin_limits[1]) - 1)
            return self._shape

    # Create getter for bin centers.
    @property
    def centers(self):
        """Return bin centers as 2 arrays (X,Y)."""
        try:
            return self._centers
        except AttributeError:
            self._centers = (self.bin_limits[0].centers(),
                             self.bin_limits[1].centers())
            return self._centers

    @property
    def center_grids(self):
        """Meshgrid of bin centers."""
        try:
            return self._center_grids
        except AttributeError:
            self._center_grids = np.meshgrid(self.bin_limits[0].centers(),
                                             self.bin_limits[1].centers())
            return self._center_grids

    @property
    def edge_grids(self):
        """Meshgrid of bin edges.

        The edge meshgrids should be used with the matplotlib.pcolor
        command.
        """
        try:
            return self._edge_grids
        except AttributeError:
            self._edge_grids = np.meshgrid(
                self.bin_limits[0], self.bin_limits[1])
            return self._edge_grids

    @property
    def widths(self):
        """Return bin widths as array."""
        try:
            return self._bin_widths
        except AttributeError:
            self._bin_widths = [
                self.bin_limits[0].widths(),
                self.bin_limits[1].widths()
            ]
            return self._bin_widths

    @property
    def sizes(self):
        """Bin sizes as 2-d array."""
        try:
            return self.bin_widths
        except AttributeError:
            self.bin_widths = np.outer(self.bin_limits[0].widths(),
                                       self.bin_limits[1].widths())
            return self.bin_widths


def get_reasonable_data_indices_for_binning(x,
                                            y=None,
                                            xscale='log',
                                            yscale='log'):
    """
    Parameters
    ----------
    x : 1D numpy array
    y : 1D numpy array, optional
    xscale : str, defaulting to 'log', optional
        only if xscale is 'log' zero values in x are omitted
    yscale : str, defaulting to 'log', optional
        only if yscale is 'log' zero values in y are omitted

    Returns
    -------
    xidx : numpy array of indices with binnable values
    yidx : numpy array of indices with binnable values, only returned if
           parameter y is given
    """
    xidx = (x != np.nan)
    if xscale == 'log':
        xidx *= (x != 0)
    if y is not None:
        yidx = (y != np.nan)
        if yscale == 'log':
            yidx = (y != 0)
        return xidx, yidx
    return xidx

# def normalize(x):
#     """Normalize a sequence.

#     Returns the sequence where each element is divided by the sum of
#     all elements. The return type is the same as type of
#     `x`. Guaranteed to work with tuples, lists, numpy arrays and
#     numpy.ma masked arrays; other types may also work but have not
#     been tested.
#     """
#     if hasattr(x, 'sum'):
#         x_sum = x.sum()
#     else:
#         x_sum = sum(x)
#     x_sum = float(x_sum)

#     if isinstance(x, (tuple, list)):
#         return type(x)(np.array(x, float) / x_sum)
#     else:
#         return x / float(x_sum)


class Error(Exception):

    """Base class for exceptions in this module."""
    pass


class ParameterError(Error):

    """Exception raised for errors in the parameter."""
    pass


class BinLimitError(Error):

    """Exception raised when bin limits are invalid."""
    pass


class DataTypeError(Error):

    """Exception raised when wrong kind of data is given."""
    pass


class _BinLimits(tuple):

    """Class that represents the bin limits.

    Implemented as a tuple with additional methods to facilitate
    construction and getting the centers and widths of the bins.
    """

    @classmethod
    def __generateLinearStart(cls, minValue, maxValue):
        """Generate the linear start for linlog bins."""
        return range(max(0, int(minValue)), min(12, int(maxValue) + 2))

    @classmethod
    def __generateLinbins(cls, left_edge, right_edge, N_bins):
        """Generate linear bins.

        Parameters
        ----------
        left_edge: int or float
            The left edge of the first bin
        right_edge: int or float
            The right edge of the last bin
        N_bins : integer
            Number of bins to generate

        Return
        ------
        limit_seq : list
            The bin limits. len(limit_seq) = N_bins + 1
        """
        return list(np.linspace(left_edge, right_edge, N_bins + 1))

    @classmethod
    def __generateLogbins(cls, left_edge, right_edge, factor, uselinear=True):
        """Generate logarithmic bins.

        The upper bound of each bin is created by multiplying the
        lower bound with factor. The lower bound of the first bin
        chosen so that left_edge is the mid point of the first
        bin. Bins are created until right_edge fits into the last
        bin.

        Parameters
        ----------
        left_edge: int or float (int if useLinear=True)
            The left edge of the first bin.
        right_edge: int or float (int if useLinear=True)
            The largest value that must fit into the last bin.
        factor : float
            The factor for increasing bin limits.
        uselinear : boolean (True)
            If True, linear bins will be used between
            min(1, left_edge) and max(10, right_edge). Each
            linear bin will include one integer.

        Return
        ------
        limit_seq : list
            The bin limits. len(limit_seq) = N_bins + 1
        """

        if uselinear:
            bins = cls.__generateLinearStart(left_edge, right_edge)
        else:
            bins = [left_edge]

        i = len(bins)
        while bins[i - 1] < right_edge:
            bins.append(bins[i - 1] * factor)
            i += 1

        return bins

    @classmethod
    def __generateMaxLogbins(cls, left_edge, right_edge):
        """Generate as many logarithmic bins as possible.

        Construct logarthmic bins from left_edge to right_edge
        so that each bin contains at least one integer.

        Parameters
        ----------
        left_edge: int
            The left edge of the first bin.
        right_edge: int
            The right edge of the last bin.

        Return
        ------
        limit_seq : list
            The bin limits. len(limit_seq) = N_bins + 1

        Notes
        -----
        The number of bins grows quickly when left_edge is increased
        beyond 1. It is a good idea to check that the resulting bins
        are still suitable for your purpose.
        """
        # Initial values
        if left_edge > right_edge:
            return []
        max_bins = right_edge - left_edge + 1
        if max_bins == 1:
            return [left_edge, right_edge]

        # Find the integer that diffs factor size.
        i, factor = 1, 0
        cmp_factor = sqrt(float(left_edge + 1) / left_edge)
        while factor < cmp_factor and i < max_bins:
            factor = cmp_factor
            i += 1
            cmp_factor = (float(left_edge + i) / left_edge) ** (1.0 / (i + 1))

        # Calculate the correct number of bins and the exact factor so
        # that this number of bins is reached.
        N_bin = int(np.log(float(right_edge) / left_edge) / np.log(factor))
        factor = (float(right_edge) / left_edge) ** (1.0 / N_bin)

        # Create bins.
        bins = [left_edge]
        for i in range(N_bin):
            bins.append(bins[i] * factor)
        bins[-1] = right_edge  # Make sure the last bin is exactly right_edge.
        return bins

    @classmethod
    def __check_parameters(cls, dataType, minValue, maxValue, binType, param):
        """Make sure the construction parameters are valid."""

        if dataType not in (int, float):
            raise ParameterError("'dataType' must be int or float")

        if minValue >= maxValue and binType != 'custom':
            raise BinLimitError("minValue must be larger than maxValue.")

        if binType in ('log', 'logarithmic', 'maxlog'):
            if minValue <= 0:
                # minValue must be strictly positive in logarithmic bins.
                raise BinLimitError("minValue must be strictly positive "
                                    "with '%s' bin type." % (binType,))

        if binType in ('linlog', 'linmaxlog'):
            # minValue must be non-negative in lin-log bins.
            if minValue < 0:
                raise BinLimitError("minValue must be non-negative "
                                    "with '%s' bin type." % (binType,))

            # The lower limit in linlog types must be integer, and
            # also the upper limit if it is below 10.5.
            if minValue < 11 and not isinstance(minValue, int):
                raise BinLimitError("When minValue < 11, minValue must be integer "
                                    "with '%s' bin type." % (binType,))

            if maxValue < 11 and not isinstance(maxValue, int):
                raise BinLimitError("When maxValue < 11, maxValue must be integer "
                                    "with '%s' bin type." % (binType,))

        if binType in ('log', 'logarithmic', 'linlog'):
            if param <= 1:
                # factor must be larger than 1.
                raise ParameterError("factor (param) must be larger than"
                                     " 1 with '%s' bin type." % (binType,))

        if binType in ('lin', 'linear'):
            if not isinstance(param, int) or param <= 0:
                raise ParameterError("Number of bins (param) must "
                                     "be a positive integer.")

        if binType == 'custom' and not (np.diff(param) > 0).all():
                raise ParameterError("Bin limits (param) must be a strictly "
                                     "increasing sequence.")

    @classmethod
    def __create_bins(cls, dataType, minValue, maxValue, binType, param):
        """Construct bins."""

        # Create bins
        if (binType == 'lin' or binType == 'linear'):
            limit_seq = cls.__generateLinbins(minValue, maxValue, param)

        elif (binType == 'log' or binType == 'logarithmic'):
            limit_seq = cls.__generateLogbins(minValue, maxValue, param, False)

        elif binType == 'linlog':
            limit_seq = cls.__generateLogbins(minValue, maxValue, param, True)

        elif binType == 'maxlog' and param is None:
            limit_seq = cls.__generateMaxLogbins(minValue, maxValue)

        elif binType == 'maxlog':
            limit_seq = cls.__generateMaxLogbins(minValue, maxValue)

        elif binType == 'linmaxlog':
            limit_seq = cls.__generateLinearStart(minValue, maxValue)
            if maxValue > limit_seq[-1]:
                lastLinValue = limit_seq[-1]
                limit_seq = limit_seq[:-1]
                limit_seq.extend(
                    cls.__generateMaxLogbins(lastLinValue, maxValue))

        elif binType == 'custom' and (len(param) > 0):
            limit_seq = param
            minValue = param[0]
            maxValue = param[-1]
        else:
            raise ParameterError("Unidentified parameter combination.")

        return minValue, maxValue, limit_seq

    def __new__(cls, dataType, minValue, maxValue, binType, param=None):
        """Initialize bin limits.

        The parameters are identical to those used in Bins.__init__,
        so see that method for explanation.

        If dataType is None, it is inferred from minValue and
        maxValue: if both are of type int, dataType is also int,
        otherwise float.
        """
        # Convert dataType to proper type if string given.
        if isinstance(dataType, str) and dataType in ('int', 'float'):
            dataType = eval(dataType)

        # Convert binType to lower case (just in case).
        binType = binType.lower()

        # Make sure the input parameters are valid.
        cls.__check_parameters(dataType, minValue, maxValue, binType, param)

        # Get the bin limits.
        minValue, maxValue, limit_seq = cls.__create_bins(
            dataType, minValue, maxValue, binType, param)

        # Initialize tuple with the bin limits.
        obj = super(_BinLimits, cls).__new__(cls, limit_seq)

        # Set object variables and return
        obj.minValue, obj.maxValue = minValue, maxValue
        obj.dataType = dataType
        return obj

    def centers(self):
        """Return bin centers as array."""
        bin_centers = np.zeros(len(self) - 1, float)
        if self.dataType == int:
            for i in range(len(self) - 2):
                bin_centers[i] = 0.5 * (ceil(self[i + 1]) - 1 + ceil(self[i]))
            bin_centers[-1] = 0.5 * (floor(self[-1]) + ceil(self[-2]))
        else:
            for i in range(len(self) - 1):
                bin_centers[i] = 0.5 * (self[i + 1] + self[i])
        return bin_centers

    def widths(self):
        """Return bin widths as array."""
        if self.dataType == int:
            bin_widths = np.zeros(len(self) - 1, int)
            for i in range(len(self) - 2):
                bin_widths[i] = (ceil(self[i + 1]) - ceil(self[i]))
            bin_widths[-1] = (floor(self[-1]) + 1 - ceil(self[-2]))
        else:
            bin_widths = np.zeros(len(self) - 1, float)
            for i in range(len(self) - 1):
                bin_widths[i] = float(self[i + 1] - self[i])
        return bin_widths
