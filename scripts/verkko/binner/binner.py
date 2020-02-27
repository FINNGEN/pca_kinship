#!/usr/bin/python2.6
#
# Lauri Kovanen, 2009-2011 (lauri.kovanen@gmail.com)
# Department of Biomedical Engineering and Computational Science
# Aalto University School of Science
"""Bin data into 1 or 2-dimensional bins.

   Bins: Bin data into 1-dimensional bins.
   Bins2D: Bin data into 1-dimensional bins.
   normalize(): Divide a sequence by its sum.


Note that all 1-d binning methods (Bins.bin_*) can be used either by
giving the data (or data iterator) as input, or as
coroutines. Coroutines are used when input data is not given, and
single data points are the given with the send()-method.

Example of giving data as input:
   bins = binner.Bins(int, 0, 200, 'lin', 10)
   binned_data = bins.bin_count(data)

Example of the coroutine interface:
   bins = binner.Bins(int, 0, 200, 'lin', 10)
   bin_counter = bins.bin_count()
   for x in data:
       binned_data = bin_counter.send(x)
   # Calling next() or send(None) returns the current result.
   binned_data = bin_counter.next()
"""

from __future__ import print_function
from math import ceil, floor, sqrt
from operator import itemgetter
import numpy as np
from verkko.misc import data_utils


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


def normalize(x):
    """Normalize a sequence.

    Returns the sequence where each element is divided by the sum of
    all elements. The return type is the same as type of
    `x`. Guaranteed to work with tuples, lists, numpy arrays and
    numpy.ma masked arrays; other types may also work but have not
    been tested.
    """
    if hasattr(x, 'sum'):
        x_sum = x.sum()
    else:
        x_sum = sum(x)
    x_sum = float(x_sum)

    if isinstance(x, (tuple, list)):
        return type(x)(np.array(x, float) / x_sum)
    else:
        return x / float(x_sum)


class _binFinder(object):

    """
    Find the correct bin for a given value by searching through all
    bins with binary search. The correct bin is bin i with
    bin_limits[i] <= value < bin_limits[i+1], except for the last bin,
    for which the right edge is also included.

    If value < bin_limits[0], the returned bin index is negative, and
    if value > bin_limits[-1], the returned bin index is greater than
    len(bins)-2.

    Finding the correct bin takes log(N) time. If it is possible to
    find the correct bin in constant time (as is the case with for
    instance linear and logarithmic bins), you should inherit this
    class and reimplement the __init__ and __call__ methods.

    Parameters
    ----------
    bin_limits : sequence
        The limits of the bins used in binning. If there are N bins,
        len(bin_limits) is N+1.
    value : integer or float
        The value for which the correct bin should be located.

    Returns
    -------
    bin_index : integer
        The index of the bin for value.

    """

    def __init__(self, bin_limits):
        self.right_edge = bin_limits[-1]
        self.last_bin = len(bin_limits) - 2
        # Add auxilary bins to front and back. This way we can return
        # -1 if the value is below the true first bin limit and N if
        # the value is beyond the actual last bin.
        self.bin_limits = ([bin_limits[0] - 1] + list(bin_limits) +
                           [bin_limits[-1] + 1])

    def __call__(self, value):
        if value == self.right_edge:
            return self.last_bin

        lo, hi = 0, len(self.bin_limits) - 1
        while (lo < hi - 1):
            mid = (lo + hi) / 2
            # print "bin_limits["+str(lo)+"] = ", str(self.bin_limits[lo]),
            # print "  bin_limits["+str(mid)+"] = ", str(self.bin_limits[mid]),
            # print "  bin_limits["+str(hi)+"] = ", str(self.bin_limits[hi])
            #raw_input("Press any key ...")
            if self.bin_limits[mid] <= value:
                lo = mid
            if self.bin_limits[mid] > value:
                hi = mid
        return lo - 1


class _linBinFinder(_binFinder):

    """
    Return the correct bin for a given value with linear bins.

    If value < bin_limits[0], the returned index is negative, and if
    value > bin_limits[-1], the returned index is > len(bins)-2. The
    right edge of the last bin is included in the last bin.

    Parameters
    ----------
    bin_limits : sequence
        The limits of the bins used in binning. If there are N bins,
        len(bin_limits) is N+1.
    value : integer or float
        The value for which the correct bin should be located.

    Returns
    -------
    bin_index : integer
        The index of the bin for value.
    """

    def __init__(self, bin_limits):
        self.right_edge = bin_limits[-1]
        self.last_bin = len(bin_limits) - 2
        self.left_edge = bin_limits[0]
        self.diff = float(bin_limits[-1] - self.left_edge) / (
            len(bin_limits) - 1)
        # print "lin self.diff:", repr(self.diff)

    def __call__(self, value):
        if value == self.right_edge:
            return self.last_bin
        # print "0:", repr(value), repr(self.left_edge)
        # print "1:", repr(value - self.left_edge)
        # print "2:", repr((value - self.left_edge)/self.diff)
        # print "3:", repr(floor( (value - self.left_edge)/self.diff ))
        return int(floor((value - self.left_edge) / self.diff))


class _logBinFinder(_binFinder):

    """
    Return the correct bin for a given value with logarithmic bins.

    If value < bin_limits[0], the returned index is negative, and if
    value > bin_limits[-1], the returned index is greater than
    len(bin_limits)-2.

    Parameters
    ----------
    bin_limits : sequence
        The limits of the bins used in binning. If there are N bins,
        len(bin_limits) is N+1.
    value : integer or float
        The value for which the correct bin should be located.

    Returns
    -------
    bin_index : integer
        The index of the bin for value.
    """

    def __init__(self, bin_limits):
        # Note that _logBinFinder is just _linBinFinder with logarithmic
        # bin_limits and values. However, that simple implementation
        # doesn't work because of floating point inaccuracies; with
        # floats there are cases when log(a)-log(b) != log(a/b).

        self.right_edge = bin_limits[-1]
        self.last_bin = len(bin_limits) - 2

        self.left_edge = bin_limits[0]
        self.N_bins = len(bin_limits) - 1
        self.ratio = np.log2(float(bin_limits[-1]) / bin_limits[0])
        # print "log2 self.ratio:", repr(self.ratio)

    def __call__(self, value):
        if value == self.right_edge:
            return self.last_bin
        try:
            return int(floor(self.N_bins * np.log2(float(value) / self.left_edge)
                             / self.ratio))
        except (OverflowError, ValueError):
            # value = 0 gives OverflowError, negative values give
            # ValueError.  Since these are always below the lowest bin
            # limit in a logarithmic bin, we return a negative index.
            return -1


class _linlogBinFinder(_binFinder):

    """
    Return the correct bin for a given value with linear-logarithmic
    bins: linear bins from left edge to 11, and logarithmic bins from
    thereon.

    Parameters
    ----------
    bin_limits : sequence
        The limits of the bins used in binning. If there are N bins,
        len(bin_limits) is N+1.
    value : integer or float
        The value for which the correct bin should be located.

    Returns
    -------
    bin_index : integer
        The index of the bin for value.
    """

    def __init__(self, bin_limits):
        self.bin_limit = np.min([11, bin_limits[-1]])
        self.N_lin_bins = int(self.bin_limit - bin_limits[0])
        # Create linear binner if linear bins exist.
        if self.N_lin_bins > 0:
            self.lin_bf = _linBinFinder(bin_limits[:self.N_lin_bins + 1])
        else:
            self.lin_bf = None
        # Create logarithmic binner if logarithmic bins exist.
        if bin_limits[-1] > 11:
            self.log_bf = _logBinFinder(bin_limits[self.N_lin_bins:])
        else:
            self.log_bf = None

    def __call__(self, value):
        if value < self.bin_limit or self.log_bf is None:
            return self.lin_bf(value)
        elif self.log_bf:
            return self.N_lin_bins + self.log_bf(value)
        else:
            return self.N_lin_bins


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
            bin_finder = _linBinFinder(limit_seq)

        elif (binType == 'log' or binType == 'logarithmic'):
            limit_seq = cls.__generateLogbins(minValue, maxValue, param, False)
            bin_finder = _logBinFinder(limit_seq)

        elif binType == 'linlog':
            limit_seq = cls.__generateLogbins(minValue, maxValue, param, True)
            bin_finder = _linlogBinFinder(limit_seq)

        elif binType == 'maxlog' and param is None:
            limit_seq = cls.__generateMaxLogbins(minValue, maxValue)
            bin_finder = _logBinFinder(limit_seq)

        elif binType == 'maxlog':
            limit_seq = cls.__generateMaxLogbins(minValue, maxValue)
            bin_finder = _logBinFinder(limit_seq)

        elif binType == 'linmaxlog':
            limit_seq = cls.__generateLinearStart(minValue, maxValue)
            if maxValue > limit_seq[-1]:
                lastLinValue = limit_seq[-1]
                limit_seq = limit_seq[:-1]
                limit_seq.extend(
                    cls.__generateMaxLogbins(lastLinValue, maxValue))
            bin_finder = _linlogBinFinder(limit_seq)

        elif binType == 'custom' and (len(param) > 0):
            limit_seq = param
            bin_finder = _binFinder(param)
            minValue = param[0]
            maxValue = param[-1]
        else:
            raise ParameterError("Unidentified parameter combination.")

        return minValue, maxValue, limit_seq, bin_finder

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
        minValue, maxValue, limit_seq, bin_finder = cls.__create_bins(
            dataType, minValue, maxValue, binType, param)

        # Initialize tuple with the bin limits.
        obj = super(_BinLimits, cls).__new__(cls, limit_seq)

        # Set object variables and return
        obj.minValue, obj.maxValue = minValue, maxValue
        obj.bin_finder = bin_finder
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


class Bins(object):

    """ Class for binning 1-dimensional data."""

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
        self.bin_finder = self.bin_limits.bin_finder

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

    def bin_count(self, coords=None):
        """Bin data and return the number of data points in each bin.

        Parameters
        ----------
        coords : iterable
            The data to be binned. An element coords[j] must be a
            comparable to bin limits, and is placed into the bin i
            with bins[i] <= coords[j] < bins[i+1]

        Returns
        -------
        binned_data : array
            The binned data, with length N, where binned_data[i] is
            the number of data points that fall into the bin.

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        """

        def bin_count_gen():
            binCounts = np.zeros(len(self), float)

            while True:
                # Return the count so far, get the next element.
                elem = (yield binCounts)

                # Return current result only when next() or send(None)
                # is called.
                if elem is None:
                    continue

                # Check element right here because with bin_count it is
                # not required to be a sequence.
                if (elem < self.bin_limits.minValue or
                    elem > self.bin_limits.maxValue):
                    raise BinLimitError("Value %g is not in the interval [%g, %g]."
                                        % (elem, self.bin_limits.minValue,
                                           self.bin_limits.maxValue))
                curr_bin = self.bin_finder(elem)
                try:
                    binCounts[curr_bin] += 1
                except IndexError:
                    print ("IndexError in bin_count.")
                    print ("   Value    %g" % (elem,))
                    print ("   Interval [%g, %g]" % (self.bin_limits.minValue,
                                                     self.bin_limits.maxValue))
                    print ("   curr_bin %d" % (curr_bin,))
                    raise

        # If coords is given, simply process it and return the
        # result. Otherwise return the generator iterator so that the
        # user may send the data. Also advance the generator to the
        # first yield-statement so that the user may start calling
        # send(data), and return the generator iterator.
        bingen = bin_count_gen()
        try:
            bingen.next()
        except:
            next(bingen)
        if coords is None:
            return bingen
        else:
            for elem in coords:
                bin_counts = bingen.send(elem)
            return bin_counts

    def bin_count_divide(self, coords=None):
        """
        Bin data and return the number of data points in each bin
        divided by the bin width.

        Parameters
        ----------
        coords : iterable
            The data to be binned. Each element coords[j] must be a
            comparable to bin limits, and is placed into the bin i
            with bins[i] <= coords[j] < bins[i+1]

        Returns
        -------
        binned_data : array
            The binned data, with length N, where binned_data[i] is
            the number of data points that fall into the bin.

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        """
        def bin_count_divide_gen():
            bc_gen = self.bin_count()
            elem = yield
            while True:
                elem = (yield bc_gen.send(elem) / self.widths)

        bcd_gen = bin_count_divide_gen()
        try:
            bcd_gen.next()
        except:
            next(bcd_gen)
        if coords is None:
            return bcd_gen
        else:
            # This would work too:
            #   "return self.bin_count(coords)/self.widths",
            # but just to make sure the generator iterator works as it
            # should, we use it also when `coords` is given:
            for elem in coords:
                bin_count_div = bcd_gen.send(elem)
            return bin_count_div

    def bin_sum(self, data=None):
        """
        Bin data and return the sum of data points in each bin.

        Parameters
        ----------
        data : iterable
            The data to be binned. Each element must be a pair (coord,
            value). The element is then placed into the bin i with
            bins[i] <= coord < bins[i+1]

        Returns
        -------
        binned_data : ma.masked_array
            The binned data, with length N, where binned_data[i] is
            the sum of all values that fall into the bin. The bins
            with no values are masked. To get a plain list, use
            binned_data.tolist() --- by default the missing values are
            replaced with None.

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        Raise DataTypeError if data does not consist of pairs.
        """
        def bin_sum_gen():
            binValues = np.zeros(len(self), float)
            binCounts = np.zeros(len(self), float)

            while True:
                elem = (yield np.ma.masked_array(binValues, binCounts == 0))

                # Return current result only when next() or send(None)
                # is called.
                if elem is None:
                    continue

                # Make sure the data is valid.
                self.__check_data_element(elem, 2)
                # Find the correct bin and increase count.
                curr_bin = self.bin_finder(elem[0])
                binCounts[curr_bin] += 1
                try:
                    binValues[curr_bin] += elem[1]
                except IndexError:
                    raise DataTypeError("Elements of input data must be sequences"
                                        " with length at least 2.")

        bs_gen = bin_sum_gen()
        try:
            bs_gen.next()
        except:
            next(bs_gen)
        if data is None:
            return bs_gen
        else:
            for elem in data:
                bin_sum = bs_gen.send(elem)
            return bin_sum

    def bin_sum_divide(self, data=None):
        """
        Bin data and return the sum of data points in each bin divided
        by the bin width.

        Parameters
        ----------
        data : iterable
            The data to be binned. Each element must be a pair (coord,
            value). The element is then placed into the bin i with
            bins[i] <= coord < bins[i+1]

        Returns
        -------
        binned_data : ma.masked_array
            The binned data, with length N, where binned_data[i] is
            the sum of all values that fall into the bin, divided by
            the width of the bin. The bins with no values are
            masked. To get a plain list, use binned_data.tolist().

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        Raise DataTypeError if data does not consist of pairs.
        """
        def bin_sum_divide_gen():
            bs_gen = self.bin_sum()
            elem = yield
            while True:
                elem = (yield bs_gen.send(elem) / self.widths)

        bsd_gen = bin_sum_divide_gen()
        try:
            bsd_gen.next()
        except:
            next(bsd_gen)
        if data is None:
            return bsd_gen
        else:
            for elem in data:
                bin_sum_div = bsd_gen.send(elem)
            return bin_sum_div

    def bin_average(self, data=None, variances=False):
        """
        Bin data and return the average of data points in each bin. If
        variances = True, return also the variance in each bin.

        Parameters
        ----------
        data : iterable
            The data to be binned. Each element must be a pair (coord,
            value). The element is then placed into the bin i with
            bins[i] <= coord < bins[i+1]

        Returns
        -------
        variances = False,
            binned_data : ma.masked_array
                The binned data, with length N, where binned_data[i]
                is the average of all values that fall into the
                bin. The bins with no values are masked. To get a
                plain list, use binned_data.tolist().

        variances = True,
            binned_data : tuple (ma.masked_array, ma.masked_array)
                The first element is the average in each bin as
                above, and the second element is the variance of the
                data in each bin.

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        Raise DataTypeError if data does not consist of pairs.
        """
        def bin_average_gen(variances):
            binValues = np.zeros(len(self), float)
            binSquares = np.zeros(len(self), float)
            binCounts = np.zeros(len(self), float)
            ma_averages = np.ma.masked_array(binCounts.copy(), binCounts == 0)
            ma_variances = np.ma.masked_array(binCounts.copy(), binCounts == 0)

            while True:
                if variances:
                    elem = (yield ma_averages, ma_variances)
                else:
                    elem = (yield ma_averages)

                # Return current result only when next() or send(None)
                # is called.
                if elem is None:
                    continue

                # Make sure the data is valid.
                self.__check_data_element(elem, 2)
                # Find the correct bin.
                curr_bin = self.bin_finder(elem[0])
                binCounts[curr_bin] += 1
                try:
                    binValues[curr_bin] += elem[1]
                    binSquares[curr_bin] += elem[1] ** 2
                except IndexError:
                    raise DataTypeError("Elements of input data must be sequences"
                                        " with length at least 2.")

                ma_averages[curr_bin] = binValues[
                    curr_bin] / binCounts[curr_bin]
                ma_averages.mask[curr_bin] = False
                if variances:
                    ma_variances[curr_bin] = (binSquares[curr_bin] / binCounts[curr_bin]
                                              - ma_averages[curr_bin] ** 2)
                    ma_variances.mask[curr_bin] = False

        ba_gen = bin_average_gen(variances)
        try:
            ba_gen.next()
        except:
            next(ba_gen)
        if data is None:
            return ba_gen
        else:
            for elem in data:
                ret_val = ba_gen.send(elem)
            return ret_val

    def bin_weighted_average(self, data, variances=False):
        """
        Bin data and return the weighted average of data points in
        each bin. If variances = True, return also the weighted
        variance in each bin.

        Parameters
        ----------
        data : iterable
            The data to be binned. Each element must be a triple
            (coord, value, weight). The element is then placed into
            the bin i with bins[i] <= coord < bins[i+1]

        Returns
        -------
        variances = False,
            binned_data : ma.masked_array
                The binned data, with length N, where binned_data[i]
                is the weighted average of all values that fall into
                the bin. The bins with no values are masked. To get a
                plain list, use binned_data.tolist().

        variances = True,
            binned_data : tuple (ma.masked_array, ma.masked_array)
                The first element is the weighted average in each bin
                as above, and the second element is the weighted
                variance of the data in each bin.

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        Raise DataTypeError if data does not consist of triples.
        """
        def bin_weighted_average_gen(variances):

            binValues = np.zeros(len(self), float)
            binSquares = np.zeros(len(self), float)
            binCounts = np.zeros(len(self), float)
            binWeights = np.zeros(len(self), float)

            while True:
                binValues_ma = np.ma.masked_array(binValues, binCounts == 0)
                binWeights_ma = np.ma.masked_array(binWeights, binWeights == 0)
                averages_ma = binValues_ma / binWeights_ma

                if variances:
                    elem = (yield averages_ma, binSquares / binWeights_ma - averages_ma ** 2)
                else:
                    elem = (yield averages_ma)

                # Return current result only when next() or send(None)
                # is called.
                if elem is None:
                    continue

                # Make sure the data is valid.
                self.__check_data_element(elem, 3)
                # Find the correct bin.
                curr_bin = self.bin_finder(elem[0])
                binCounts[curr_bin] += 1
                try:
                    binValues[curr_bin] += elem[1] * elem[2]
                    binSquares[curr_bin] += elem[1] ** 2 * elem[2]
                    binWeights[curr_bin] += elem[2]
                except IndexError:
                    raise DataTypeError("Elements of input data must be sequences"
                                        " with length at least 3.")

        bwa_gen = bin_weighted_average_gen(variances)
        bwa_gen.next()
        if data is None:
            return bwa_gen
        else:
            for elem in data:
                ret_val = bwa_gen.send(elem)
            return ret_val

    def bin_percentile(self, data=None, perc=(0.5,)):
        """
        Bin data and return the median in each bin.

        Parameters
        ----------
        data : iterable
            The data to be binned. Each element must be a pair (coord,
            value). The element is then placed into the bin i with
            bins[i] <= coord < bins[i+1]
        perc : sequence of floats
            The percentiles to calculate. The p:th percentile is the
            value below which one can find fraction p of the data.

        Returns
        -------
        binned_data : list of ma.masked_arrays
            The returned list will be of length len(perc), and there
            is one ma.masked_array for each p in perc. Percentiles are
            calculated and returned in increasing order. Each masked
            array contains the percentile of the data that falls into
            that bin.

            The bins with no values are masked. To get a plain list,
            use binned_data[i].tolist().

        Exceptions
        ----------
        Raise ParameterError if `perc` is empty.
        Raise ParameterError if any element of `perc` is not in (0, 1.0).
        Raise BinLimitError if any element does not fit into the bins.
        Raise DataTypeError if data does not consist of pairs.

        Notes
        -----
        Finding the percentiles of an arbitrary sequence requires
        first saving all elements. If the number of data points is
        very large, this method can take up a lot of memory.
        """
        if not perc:  # This covers both None and empty `perc`.
            raise ParameterError("No percentiles given.")
        perc = sorted(perc)

        def bin_percentile_gen(perc):
            binElements = [[] for i in range(len(self))]
            perc_arrays = [np.ma.masked_array(np.zeros(len(self), float),
                                              np.ones(len(self), float))
                           for j in range(len(perc))]

            while True:
                # Get next element and return percentiles.
                elem = (yield perc_arrays)

                if elem is None:
                    continue

                # Make sure the data is valid.
                self.__check_data_element(elem, 2)
                # Find the correct bin.
                curr_bin = self.bin_finder(elem[0])
                # Append the list of elements in the bin. BinLimitError
                # occurs if bin goes over the top, and DataTypeError
                # occurs if elem is not a sequence with length at least 2.
                try:
                    binElements[curr_bin].append(elem[1])
                except IndexError:
                    if curr_bin > len(self) - 1:
                        raise BinLimitError("Value %g is beyond the upper bin "
                                            "limit %g." %
                                            (elem, self.bin_limits.maxValue))
                    else:
                        raise DataTypeError("Elements of input data must be "
                                            "sequences with length at least 2.")

                # Sort the data in the bin that was altered.
                binElements[curr_bin] = sorted(binElements[curr_bin])

                # Find the percentiles in the bin that was
                # altered. Note that if a bin has no values,
                # data_utils.percentile returns None, which will be
                # turned into a nan by numpy.
                for i_p, p in enumerate(perc):
                    perc_arrays[i_p].mask[curr_bin] = False
                    perc_arrays[i_p][curr_bin] = data_utils.percentile(
                        binElements[curr_bin], p)

        bp_gen = bin_percentile_gen(perc)
        try:
            bp_gen.next()
        except:
            next(bp_gen)
        if data is None:
            return bp_gen
        else:
            for elem in data:
                ret_val = bp_gen.send(elem)
            return ret_val

    def bin_median(self, data=None):
        """
        Bin data and return the median in each bin.

        Parameters
        ----------
        data : iterable
            The data to be binned. Each element must be a pair (coord,
            value). The element is then placed into the bin i with
            bins[i] <= coord < bins[i+1]

        Returns
        -------
        binned_data : list of ma.masked_arrays
            Element binned_data[i] contains the median of all elements
            that fall into the bin.  The bins with no values are
            masked. To get a plain list, use binned_data.tolist().

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        Raise DataTypeError if data does not consist of pairs.

        Notes
        -----
        Finding the median of an arbitrary sequence requires first
        saving all elements. If the number of data points is very
        large, this method can take up a large amount of memory.
        """

        def bin_median_gen():
            bp_gen = self.bin_percentile(None, (0.5,))
            elem = yield
            while True:
                elem = (yield bp_gen.send(elem)[0])

        bm_gen = bin_median_gen()
        try:
            bm_gen.next()
        except:
            next(bm_gen)
        if data is None:
            return bm_gen
        else:
            for elem in data:
                ret_val = bm_gen.send(elem)
            return ret_val


class Bins2D(object):

    """ Class for binning 2-dimensional data."""

    def __init__(self, X_dataType, X_minValue, X_maxValue, X_binType, X_param,
                 Y_dataType, Y_minValue, Y_maxValue, Y_binType, Y_param):
        """Initialize bins.

        Constructs the bins to be used in binning. You can select
        between several different types of bin limits, and if none of
        these seem fit, you can also supply your own bin limits.

        Usage is identical to Bins, except that you have to supply
        necessary parameters for bins in both x- and y-directions. See
        the documentation of Bins for more information.
        """
        self.bin_limits = [_BinLimits(X_dataType, X_minValue, X_maxValue,
                                      X_binType, X_param),
                           _BinLimits(Y_dataType, Y_minValue, Y_maxValue,
                                      Y_binType, Y_param)]
        self.x_bin_finder = self.bin_limits[0].bin_finder
        self.y_bin_finder = self.bin_limits[1].bin_finder

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
    def sizes(self):
        """Bin sizes as 2-d array."""
        try:
            return self.bin_widths
        except AttributeError:
            self.bin_widths = np.outer(self.bin_limits[0].widths(),
                                       self.bin_limits[1].widths())
            return self.bin_widths

    def __check_data_element(self, elem, N):
        """Check one element of input data.

        Makes sure that len(elem) >= 2 and the it fits the bin
        limits. The required length of elem is N; this is not checked,
        but only shown in the error message if len(elem) < 2.
        """
        # Check bin limits and correct sequence type.
        try:
            if (elem[0] < self.bin_limits[0].minValue or
                elem[0] > self.bin_limits[0].maxValue):
                raise BinLimitError("X-coordinate %g is not in the interval [%g"
                                    ", %g]." % (elem[0], self.bin_limits[0].minValue,
                                                self.bin_limits[0].maxValue))
            elif (elem[1] < self.bin_limits[1].minValue or
                  elem[1] > self.bin_limits[1].maxValue):
                raise BinLimitError("Y-coordinate %g is not in the interval [%g"
                                    ", %g]." % (elem[1], self.bin_limits[1].minValue,
                                                self.bin_limits[1].maxValue))
        except TypeError:
            # TypeError occurs when data is a list and elem is
            # integer or a float. Rather surprisingly, numpy
            # raises an IndexError in the same situation; try for
            # instance creating a=numpy.array([1,2,3]) and then
            # call a[0][0].
            raise DataTypeError("Elements of input data must be sequences "
                                "with length at least %d." % (N,))

    def bin_count(self, coords):
        """
        Bin data and return the number of data points in each bin.

        Parameters
        ----------
        coords : iterable
            The data to be binned. An element coords[k] must be
            comparable to bin limits, and is placed into the bin (i,j)
            with x_bins[i] <= coords[k][0] < x_bins[i+1]
                 y_bins[j] <= coords[k][1] < x_bins[j+1]

        Returns
        -------
        binned_data : 2-d numpy array
            The binned data, where binned_data[i,j] is the number of
            data points that fall into the bin.

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        Raise DataTypeError if data does not consist of pairs.
        """
        binCounts = np.zeros(self.shape, float)

        for elem in coords:
            self.__check_data_element(elem, 2)
            x_bin = self.x_bin_finder(elem[0])
            y_bin = self.y_bin_finder(elem[1])
            binCounts[x_bin, y_bin] += 1

        return binCounts

    def bin_count_divide(self, coords):
        """
        Bin data and return the number of data points in each bin
        divided by the bin width.

        Parameters
        ----------
        coords : iterable
            The data to be binned. An element coords[k] must be
            comparable to bin limits, and is placed into the bin (i,j)
            with x_bins[i] <= coords[k][0] < x_bins[i+1]
                 y_bins[j] <= coords[k][1] < x_bins[j+1]

        Returns
        -------
        binned_data : 2-d numpy array
            The binned data, where binned_data[i,j] is the number of
            data points that fall into the bin divided by the bin
            size.

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        Raise DataTypeError if data does not consist of pairs.
        """
        return self.bin_count(coords) / self.sizes

    def bin_sum(self, data):
        """
        Bin data and return the sum of data points in each bin.

        Parameters
        ----------
        data : iterable
            The data to be binned. Each element must be a triple
            (x,y,value), where x and y are comparable to bin
            limits. value is placed into the bin (i,j) with
                x_bins[i] <= x < x_bins[i+1]
                y_bins[j] <= y < y_bins[j+1]

        Returns
        -------
        binned_data : ma.masked_array
            The binned data, with length N, where binned_data[i,j] is
            the sum of all values that fall into the bin. The bins
            with no values are masked. To get a plain list, use
            binned_data.tolist() --- by default the missing values are
            replaced with None.

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        Raise DataTypeError if data does not consist of triples.
        """
        binValues = np.zeros(self.shape, float)
        binCounts = np.zeros(self.shape, float)

        for elem in data:
            self.__check_data_element(elem, 3)
            x_bin = self.x_bin_finder(elem[0])
            y_bin = self.y_bin_finder(elem[1])
            binCounts[x_bin, y_bin] += 1
            try:
                binValues[x_bin, y_bin] += elem[2]
            except IndexError:
                raise DataTypeError("Elements of input data must be sequences"
                                    " with length at least 3.")

        return np.ma.masked_array(binValues, binCounts == 0)

    def bin_sum_divide(self, data):
        """
        Bin data and return the sum of data points in each bin divided
        by the bin width.

        Parameters
        ----------
        data : iterable
            The data to be binned. Each element must be a triple
            (x,y,value), where x and y are comparable to bin
            limits. value is placed into the bin (i,j) with
                x_bins[i] <= x < x_bins[i+1]
                y_bins[j] <= y < y_bins[j+1]

        Returns
        -------
        binned_data : ma.masked_array
            The binned data, with length N, where binned_data[i,j] is
            the sum of all values that fall into the bin, divided by
            the bin size. The bins with no values are masked. To get a
            plain list, use binned_data.tolist() --- by default the
            missing values are replaced with None.

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        Raise DataTypeError if data does not consist of triples.
        """
        return self.bin_sum(data) / self.sizes

    def bin_average(self, data, variances=False):
        """
        Bin data and return the average of data points in each bin. If
        variances = True, return also the variance in each bin.

        Parameters
        ----------
        data : iterable
            The data to be binned. Each element must be a triple
            (x,y,value), where x and y are comparable to bin
            limits. value is placed into the bin (i,j) with
                x_bins[i] <= x < x_bins[i+1]
                y_bins[j] <= y < y_bins[j+1]

        Returns
        -------
        variances = False,
            binned_data : ma.masked_array
                The binned data, with length N, where binned_data[i,j]
                is the average of all values that fall into the
                bin. The bins with no values are masked. To get a
                plain list, use binned_data.tolist() --- by default
                the missing values are replaced with None.

        variances = True,
            binned_data : tuple (ma.masked_array, ma.masked_array)
                The first element is the average in each bin as
                above, and the second element is the variance of the
                data in each bin.

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        Raise DataTypeError if data does not consist of triples.
        """

        binValues = np.zeros(self.shape, float)
        binSquares = np.zeros(self.shape, float)
        binCounts = np.zeros(self.shape, float)

        for elem in data:
            self.__check_data_element(elem, 3)
            x_bin = self.x_bin_finder(elem[0])
            y_bin = self.y_bin_finder(elem[1])
            binCounts[x_bin, y_bin] += 1
            try:
                binValues[x_bin, y_bin] += elem[2]
                binSquares[x_bin, y_bin] += elem[2] ** 2
            except IndexError:
                raise DataTypeError("Elements of input data must be sequences"
                                    " with length at least 3.")

        # Calculate averages as masked arrays.
        binCounts = np.ma.masked_array(binCounts, binCounts == 0)
        ma_averages = binValues / binCounts

        if variances:
            return ma_averages, binSquares / binCounts - ma_averages ** 2
        else:
            return ma_averages

    def bin_weighted_average(self, data, variances=False):
        """
        Bin data and return the weighted average of data points in
        each bin. If variances = True, return also the weighted
        variance in each bin.

        Parameters
        ----------
        data : iterable
            The data to be binned. Each element must be a tuple
            (x,y,value,weight), where x and y are comparable to bin
            limits. value is placed into the bin (i,j) with
                x_bins[i] <= x < x_bins[i+1]
                y_bins[j] <= y < y_bins[j+1]

        Returns
        -------
        variances = False,
            binned_data : ma.masked_array
                The binned data, with length N, where binned_data[i,j]
                is the weighted average of all values that fall into
                the bin. The bins with no values are masked. To get a
                plain list, use binned_data.tolist() --- by default
                the missing values are replaced with None.

        variances = True,
            binned_data : tuple (ma.masked_array, ma.masked_array)
                The first element is the weighted average in each bin
                as above, and the second element is the weighted
                variance of the data in each bin.

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        Raise DataTypeError if data does not consist of triples.
        """
        binValues = np.zeros(self.shape, float)
        binSquares = np.zeros(self.shape, float)
        binCounts = np.zeros(self.shape, float)
        binWeights = np.zeros(self.shape, float)

        for elem in data:
            self.__check_data_element(elem, 4)
            x_bin = self.x_bin_finder(elem[0])
            y_bin = self.y_bin_finder(elem[1])
            binCounts[x_bin, y_bin] += 1
            try:
                binValues[x_bin, y_bin] += elem[2] * elem[3]
                binSquares[x_bin, y_bin] += elem[2] ** 2 * elem[3]
                binWeights[x_bin, y_bin] += elem[3]
            except IndexError:
                raise DataTypeError("Elements of input data must be sequences"
                                    " with length at least 4.")

        # Calculate weighted average.
        binValues = np.ma.masked_array(binValues, binCounts == 0)
        binWeights[binWeights == 0] = 1.0
        ma_wAverages = binValues / binWeights

        if variances:
            return ma_wAverages, binSquares / binWeights - ma_wAverages ** 2
        else:
            return ma_wAverages

    def bin_median(self, data):
        """
        Bin data and return the median in each bin.

        Parameters
        ----------
        data : iterable
            The data to be binned. Each element must be a triple
            (x,y,value), where x and y are comparable to bin
            limits. value is placed into the bin (i,j) with
                x_bins[i] <= x < x_bins[i+1]
                y_bins[j] <= y < y_bins[j+1]

        Returns
        -------
        binned_data : ma.masked_array
            The binned data, with length N, where binned_data[i,j] is
            the median of all values that fall into the bin. The bins
            with no values are masked. To get a plain list, use
            binned_data.tolist() --- by default the missing values are
            replaced with None.

        Exceptions
        ----------
        Raise BinLimitError if any element does not fit into the bins.
        Raise DataTypeError if data does not consist of triples.

        Notes
        -----
        Finding the median of an arbitrary sequence requires first
        saving all elements. If the number of data points is very
        large, this method can take up a large amount of memory.
        """

        binElements = [[[] for i in range(self.shape[1])]
                       for j in range(self.shape[0])]

        for elem in data:
            self.__check_data_element(elem, 3)
            x_bin = self.x_bin_finder(elem[0])
            y_bin = self.y_bin_finder(elem[1])
            try:
                binElements[x_bin][y_bin].append(elem[2])
            except IndexError:
                raise DataTypeError("Elements of input data must be sequences"
                                    " with length at least 3.")

        binMedians = np.zeros(self.shape, float)

        # Find medians for each bin. If a bin has no values, np.median
        # returns nan.
        for i, i_elements in enumerate(binElements):
            for j, j_elements in enumerate(i_elements):
                binMedians[i, j] = np.median(j_elements)

        return np.ma.masked_array(binMedians, np.isnan(binMedians))


if __name__ == '__main__':
    """Run unit tests if called."""
    from tests.test_binner import *
    unittest.main()
