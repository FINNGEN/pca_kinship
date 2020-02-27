import unittest
from operator import itemgetter
import binner
import numpy # Used only for testing with numpy arrays.
import warnings


class TestBins(unittest.TestCase):
    def setUp(self):
        # Construct good data
        self.weights = [7.5,   2,   4,  3, 100, 0.1, 0.9,  2,  3]
        self.coords =  [ 20,   5,   5, 10,  16,   7,   7, 11, 12]
        self.values =  [ 10, 0.5, 3.5,  1, 1.6, 100,   2, 10, 20]
        self.data = zip(self.coords, self.values)
        self.weighted_data = zip(self.coords, self.values, self.weights)

        # Construct bad data
        self.bad_coords_A = [20,   -1,   5, 10,  16,   7,   7, 11, 12]
        self.bad_coords_B = [20,   2,   5, 10,  30,   7,   7, 11, 12]
        self.bad_data_A = zip(self.bad_coords_A, self.values)
        self.bad_data_B = zip(self.bad_coords_B, self.values)
        self.bad_wdata_A = zip(self.bad_coords_A, self.values, self.weights)
        self.bad_wdata_B = zip(self.bad_coords_B, self.values, self.weights)
        
        self.bins = binner.Bins(int, 0, 20, 'lin', 10)
        
    def __test_0printer(self):
        """Doesn't test anything, just prints info.

        NOTE! REMOVE THE UNDERSCORES BEFORE 'test' IN THE NAME TO RUN
        THIS METHOD DURING FULL TEST.
        """
        deco = [(c,d,w) for c,d,w in zip(self.coords, self.values, self.weights)]
        deco.sort()
        print "Coords : " + reduce(lambda x,y: str(x)+"\t"+str(y), map(itemgetter(0), deco))
        print "Values : " + reduce(lambda x,y: str(x)+"\t"+str(y), map(itemgetter(1), deco))
        print "Weights: " + reduce(lambda x,y: str(x)+"\t"+str(y), map(itemgetter(2), deco))
        print
        print "Bin limits: ", self.bins.bin_limits
        print "Bin centers:", self.bins.centers
        print "Bin widths: ", self.bins.widths
        print

    def test_binFinder(self):
        bins = [1, 5, 9, 13, 17, 21, 25, 29, 33]
        bf = binner._binFinder(bins)
        self.assertEqual(bf(1), 0)
        self.assertEqual(bf(2), 0)
        self.assertEqual(bf(21), 5)
        self.assertEqual(bf(32), 7)
        self.assert_(bf(33), len(bins)-2)
        
        bins = [1, 100, 1000, 10000, 10010, 10020, 10100, 10200]
        bf = binner._binFinder(bins)
        self.assertEqual(bf(1.001), 0)
        self.assertEqual(bf(99.99), 0)
        self.assertEqual(bf(10009.9), 3)
        self.assertEqual(bf(10010), 4)
        self.assertEqual(bf(10011), 4)
        self.assertEqual(bf(10200), len(bins)-2)

        bins = [-5.87, -1.74, 1.0/3, 1.1/3, 78.0001, 100124]
        bf = binner._binFinder(bins)
        self.assertEqual(bf(-5.86), 0)
        self.assertEqual(bf(-1.74), 1)
        self.assertEqual(bf(0.3333), 1)
        self.assertEqual(bf(1.0/3), 2)
        self.assertEqual(bf(1.05/3), 2)
        self.assert_(bf(100124) == len(bins)-2)
        
        self.assert_(bf(-5.88) < 0)
        self.assert_(bf(-78) < 0)
        self.assert_(bf(200000) > len(bins)-2)

    def test_linBinFinder(self):
        bins = [0.5,  1.5,  2.5,  3.5,  4.5,  5.5]
        bf = binner._linBinFinder(bins)
        self.assertEqual(bf(0.5), 0)
        self.assertEqual(bf(1), 0)
        self.assertEqual(bf(2.5), 2)
        self.assertEqual(bf(5), 4)
        self.assertEqual(bf(5.0), 4)
        self.assertEqual(bf(5.5), len(bins)-2)

        # Test out of bounds values
        self.assert_(bf(0.499) < 0)
        self.assert_(bf(0) < 0)
        self.assert_(bf(-100) < 0)
        self.assert_(bf(6) > len(bins)-2)
        self.assert_(bf(111) > len(bins)-2)

        bins = [2.5,   7. ,  11.5,  16. ,  20.5]
        bf = binner._linBinFinder(bins)
        self.assertEqual(bf(2.5), 0)
        self.assertEqual(bf(7.6), 1)
        self.assertEqual(bf(18), 3)
        self.assertEqual(bf(20.49), 3)
        self.assertEqual(bf(20.5), len(bins)-2)

    def test_logBinFinder(self):
        bins = [ 1, 2, 4, 8, 16, 32, 64 ]
        bf = binner._logBinFinder(bins)
        self.assertEqual(bf(1.5), 0)
        self.assertEqual(bf(2), 1)
        self.assertEqual(bf(2.1), 1)
        self.assertEqual(bf(10), 3)
        self.assertEqual(bf(63.429151), 5)
        self.assertEqual(bf(64), len(bins)-2)        

        # Test out of bounds values
        self.assert_(bf(0.999) < 0)
        oldErrs = numpy.seterr(all="ignore")
        self.assert_(bf(0) < 0)
        self.assert_(bf(-100) < 0)
        numpy.seterr(**oldErrs)
        self.assert_(bf(111) > len(bins)-2)

        bins = [ 5, 7.5, 11.25, 16.875, 25.3125, 37.96875 ]
        bf = binner._logBinFinder(bins)
        self.assertEqual(bf(5.5), 0)
        self.assertEqual(bf(25.3124), 3)
        self.assertEqual(bf(25.3125), 4)
        self.assertEqual(bf(37.96875), len(bins)-2)

    def test_linlogBinFinder(self):
        # Test typical case.
        bins = range(1, 12) + [22, 44, 88]
        bf = binner._linlogBinFinder(bins)
        self.assertEqual(bf(1), 0)
        self.assertEqual(bf(7), 6)
        self.assertEqual(bf(10), 9)
        self.assertEqual(bf(11), 10)
        self.assertEqual(bf(20), 10)
        self.assertEqual(bf(80), 12)
        self.assertEqual(bf(88), len(bins)-2)
        
        # Test out of bounds values
        self.assert_(bf(0.499) < 0)
        self.assert_(bf(0) < 0)
        self.assert_(bf(-100) < 0)
        self.assert_(bf(111) > len(bins)-2)

        # Smallest bin value > 1.
        bins = range(6,12) + [44, 176, 704, 2816]
        bf = binner._linlogBinFinder(bins)
        self.assertEqual(bf(6), 0)
        self.assertEqual(bf(10), 4)
        self.assertEqual(bf(44), 6)
        self.assertEqual(bf(100), 6)
        self.assertEqual(bf(176), 7)
        self.assertEqual(bf(2816), len(bins)-2)

        # Largest bin value < 11
        bins = [2, 3, 4, 5, 6]
        bf = binner._linlogBinFinder(bins)
        self.assert_(bf(1) < 0)
        self.assertEqual(bf(2), 0)
        self.assertEqual(bf(4), 2)
        self.assertEqual(bf(6), 3)
        self.assert_(bf(7) > 3)
        self.assert_(bf(100) > 3)

    def test_conformity(self):
        """ The general binary search in the base class binFinder
        should always give exactly the same result as any
        implementation for any special case.
        """
        bins = [6, 7, 8, 9, 10, 11, 44, 176, 704, 2816]
        bf_A = binner._linlogBinFinder(bins)
        bf_B = binner._binFinder(bins)
        for value in [6, 11, 100, 176, 704, 2815]:
            self.assertEqual(bf_A(value), bf_B(value))

    def test_Bins_errors(self):
        """Check that input parameters are correct when constructing bins."""
        # minValue == maxValue
        self.assertRaises(binner.BinLimitError, binner.Bins, int, 7, 7, 'lin', 10)
        # Logarithmic bins with minValue <= 0
        self.assertRaises(binner.BinLimitError, binner.Bins, int, 0, 1, 'log', 1.5)
        # Logarithmic bins with factor <= 1
        self.assertRaises(binner.ParameterError, binner.Bins, int, 1, 10, 'linlog', 1)
        self.assertRaises(binner.ParameterError, binner.Bins, int, 1, 10, 'logarithmic', 0.9)
        # Linear-logarithmic bins with non-integer bin limits
        self.assertRaises(binner.BinLimitError, binner.Bins, int, 1.5, 20, 'linlog', 2)
        self.assertRaises(binner.BinLimitError, binner.Bins, int, 1, 10.1, 'linmaxlog')
        # Linear bins with <= 0 bins or with a floating point number of bins.
        self.assertRaises(binner.ParameterError, binner.Bins, int, 1, 10, 'linear', 0)
        self.assertRaises(binner.ParameterError, binner.Bins, int, 1, 10, 'linear', 10.5)
        # Custom bins with non-increasing sequence
        bin_lim = [0, 1, 1.5, 3, 2]
        self.assertRaises(binner.ParameterError, binner.Bins, int, 1, 10, 'custom', bin_lim)
        bin_lim = [0, 1, 1.5, 3, 3]
        self.assertRaises(binner.ParameterError, binner.Bins, int, 1, 10, 'custom', bin_lim)
        
    def test_Linbins_int(self):
        N_bin = 3
        b = binner.Bins(int, 1, 5, 'lin', N_bin)
        self.assertEqual(len(b), N_bin)
        self.assertEqual(b.widths.sum(), 5)

        N_bin = 10
        b = binner.Bins(int, 1, 10, 'lin', N_bin)
        self.assertEqual(len(b), N_bin)
        self.assertEqual(b.widths.tolist(), [1]*N_bin)
        self.assertEqual(b.centers.tolist(), range(1,11))

        N_bin = 25
        b = binner.Bins(int, 26, 75, 'lin', N_bin)
        self.assertEqual(len(b), N_bin)
        self.assertEqual(b.widths.tolist(), [2]*N_bin)
        self.assertEqual(b.centers.tolist(), [x+0.5 for x in range(26,75,2)])

    def test_Linbins_float(self):
        N_bin = 3
        b = binner.Bins(float, 1, 5, 'lin', N_bin)
        self.assertEqual(len(b), N_bin)
        self.assertEqual(b.widths.sum(), 4.0)
        exp_centers = [1+(5-1)/6.0, 1+3*(5-1)/6.0, 1+5*(5-1)/6.0]
        for c, ce in zip(exp_centers, b.centers):
            self.assertEqual("%.8f" % c, "%.8f" % ce)
        for w, we in zip([(5-1)/3.0]*N_bin, b.widths):
            self.assertEqual("%.8f" % w, "%.8f" % we)

    def test_Logbins_int(self):
        b = binner.Bins(int, 1, 5, 'log', 2)
        expected_bins = (1, 2, 4, 8)
        self.assertEqual(b.bin_limits, expected_bins)
        expected_centers = [1, 2.5, 6]
        self.assertEqual(b.centers.tolist(), expected_centers)
        expected_widths = [1, 2, 5]
        self.assertEqual(b.widths.tolist(), expected_widths)

    def test_Linlogbins_int(self):
        b = binner.Bins(int, 7, 12, 'linlog', 2)
        exp_bins = (7, 8, 9, 10, 11, 22)
	self.assertEqual(b.bin_limits, exp_bins)
        exp_widths = [1,1,1,1,12]
        self.assertEqual(b.widths.tolist(), exp_widths)
        exp_centers = [7, 8, 9, 10, (22+11)/2.0]
        self.assertEqual(b.centers.tolist(), exp_centers)
        
        b = binner.Bins(float, 0, 11, 'linlog', 2)
        exp_widths = [1.0]*len(b)
        self.assertEqual(b.widths.tolist(), exp_widths)

    def test_Linlogbinfinder_zeroBehaviour(self):
        b = binner.Bins(float, 0, 20, 'linlog', 2)
        data = [(1,7), (0,8)]
        exp_result = [8.0, 7.0] + [None]*10
        self.assertEqual(b.bin_average(data).tolist(), exp_result)
        data.append((-1,2))
        self.assertRaises(binner.BinLimitError, b.bin_average, data)

        b = binner.Bins(float, 1, 20, 'linlog', 2)
        data = [(1,7), (0,8)]
        self.assertRaises(binner.BinLimitError, b.bin_average, data)

    def test_Linlogbins_float(self):
        b = binner.Bins(float, 7, 12, 'linlog', 2)
        exp_bins = (7, 8, 9, 10, 11, 22)
	self.assertEqual(b.bin_limits, exp_bins)
        exp_widths = [1,1,1,1,11]
        self.assertEqual(b.widths.tolist(), exp_widths)
        exp_centers = [7.5, 8.5, 9.5, 10.5, (22+11)/2.0]
        self.assertEqual(b.centers.tolist(), exp_centers)

        b = binner.Bins(float, 0, 11, 'linlog', 2)
        exp_widths = [1.0]*len(b)
        self.assertEqual(b.widths.tolist(), exp_widths)

    def test_Maxlogbins(self):
        minValue = 1
        maxValue = 20
        b = binner.Bins(int, minValue, maxValue, 'maxlog')
        self.assert_(b.bin_limits[0] == minValue)
        self.assert_(b.bin_limits[-1] == maxValue)
        for w in b.widths:
            self.assert_(w >= 1)

    def test_custom_bins(self):
        """Test arbitrary bins."""
        bins = binner.Bins(int, 0, 0, 'custom', [1, 2.5, 4.5])
        self.assertEqual(bins.widths.tolist(), [2, 2])
        self.assertEqual(bins.centers.tolist(), [1.5, 3.5])

        bins = binner.Bins(int, 0, 0, 'custom', [-1.5, 0, 6.5, 6.6])
        self.assertEqual(bins.widths.tolist(), [1, 7, 0])

        bins = binner.Bins(int, 0, 0, 'custom', [-10, -4, -1.5])
        self.assertEqual(bins.widths.tolist(), [6, 3])
        
    def test_Count(self):
        # Check correct result
        binned_data = self.bins.bin_count(self.coords)
        expected_result = [0,0,2,2,0,2,1,0,1,1]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.BinLimitError, self.bins.bin_count, self.bad_coords_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_count, self.bad_coords_B)

    def test_CountDiv(self):
        # Check correct result
        binned_data = self.bins.bin_count_divide(self.coords)
        expected_result = [0,0,1,1,0,1,0.5,0,0.5,1.0/3]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.BinLimitError, self.bins.bin_count_divide, self.bad_coords_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_count_divide, self.bad_coords_B)

    def test_Sum(self):
        # Check correct result
        binned_data = self.bins.bin_sum(self.data)
        expected_result = [None,None,4,102,None,11,20,None,1.6,10]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_sum, self.coords)
        self.assertRaises(binner.BinLimitError, self.bins.bin_sum, self.bad_data_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_sum, self.bad_data_B)

        # Make sure numpy arrays work as promised.
        binned_data = self.bins.bin_sum(numpy.array(self.data))
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_sum, numpy.array(self.coords))
        self.assertRaises(binner.BinLimitError, self.bins.bin_sum, numpy.array(self.bad_data_A))
        self.assertRaises(binner.BinLimitError, self.bins.bin_sum, numpy.array(self.bad_data_B))

    def test_SumDivide(self):
        # Check correct result
        binned_data = self.bins.bin_sum_divide(self.data)
        expected_result = [None,None,2,51,None,5.5,10,None,0.8,10.0/3]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_sum_divide, self.coords)
        self.assertRaises(binner.BinLimitError, self.bins.bin_sum_divide, self.bad_data_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_sum_divide, self.bad_data_B)

    def test_Average(self):
        # Check correct result
        binned_data = self.bins.bin_average(self.data)
        expected_result = [None,None,2,51,None,5.5,20,None,1.6,10]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_average, self.coords)
        self.assertRaises(binner.BinLimitError, self.bins.bin_average, self.bad_data_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_average, self.bad_data_B)

    def test_Average_Variance(self):
        # Check correct result
        binned_data = self.bins.bin_average(self.data, True)
        expected_average = [None,None,2,51,None,5.5,20,None,1.6,10]
        expected_variance = [None,None,2.25,2401.0,None,20.25,0.0,None,0.0,0.0]
        self.assertEqual(binned_data[0].tolist(), expected_average)
        self.assertEqual(binned_data[1].tolist(), expected_variance)

    def test_WeightedAverage(self):
        # Check correct result
        binned_data = self.bins.bin_weighted_average(self.weighted_data)
        expected_result = [None,None,2.5,11.8,None,23./5,20,None,1.6,10]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_weighted_average, self.coords)
        self.assertRaises(binner.DataTypeError, self.bins.bin_weighted_average, self.data)
        self.assertRaises(binner.BinLimitError, self.bins.bin_weighted_average, self.bad_wdata_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_weighted_average, self.bad_wdata_B)

    def test_WeightedAverage_Variance(self):
        # Check correct result
        oldErrs = numpy.seterr(invalid="ignore")
        binned_data = self.bins.bin_weighted_average(self.weighted_data, True)
        numpy.seterr(**oldErrs)

        expected_average = [None,None,2.5,11.8,None,23./5,20,None,1.6,10]
        expected_variance = [None,None,2.0, 0.9*4+1000-11.8**2,None,
                             (3+200)/5.0-(23./5)**2,0.0,None,0.0,0.0]
        self.assertEqual(binned_data[0].tolist(), expected_average)
        self.assertEqual(binned_data[1].tolist(), expected_variance)

    def test_Median(self):
        # Check correct result
        new_data = zip(self.coords + [4,11], self.values + [1,30])
        binned_data = self.bins.bin_median(new_data)
        expected_result = [None,None,1,51,None,10,20,None,1.6,10]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_median, self.coords)
        self.assertRaises(binner.BinLimitError, self.bins.bin_median, self.bad_data_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_median, self.bad_data_B)

    def test_Percentile(self):
        # Check correct result
        data = [(1, 1), (1, 2), (1, 3), (1, 4), (1, 5),
                (6, 10), (6, 30), (7, 100), (7, 300),
                (14, 42),
                (20, 1.1), (20, 1.5), (20, 1.5), (20, 1.9)]
        perc = (0.2, 0.25, 0.75, 0.8)
        expected_res = ([0.2*1+0.8*2,None,None,0.4*10+0.6*30,None,None,None,42,None,0.4*1.1+0.6*1.5],
                        [2.0,None,None,0.25*10+0.75*30,None,None,None,42,None,0.25*1.1+0.75*1.5],
                        [4.0,None,None,0.75*100+0.25*300,None,None,None,42,None,0.75*1.5+0.25*1.9],
                        [0.8*4+0.2*5,None,None,0.6*100+0.4*300,None,None,None,42,None,0.6*1.5+0.4*1.9])

        binned_data = self.bins.bin_percentile(data, perc)
        for res, exp_res in zip(binned_data, expected_res):
            for r, er in zip(res, exp_res):
                if r is not None and er is not None:
                    self.assertEqual("%.8f" % r, "%.8f" % er)

        # Check exceptions
        self.assertRaises(binner.ParameterError, self.bins.bin_percentile, data, None)
        self.assertRaises(binner.DataTypeError, self.bins.bin_percentile, self.coords, (0.1, 0.2))
        self.assertRaises(binner.BinLimitError, self.bins.bin_percentile, self.bad_data_A, (0.1, 0.2))
        self.assertRaises(binner.BinLimitError, self.bins.bin_percentile, self.bad_data_B, (0.1, 0.2))


class TestBins2D(unittest.TestCase):
    def setUp(self):
        self.x_coords = [1, 1, 2, 4, 4, 3, 1, 2]
        self.y_coords = [1, 4, 3, 5, 6, 1, 8, 6]
        values = [1, 3, 7, 5, 6, 3, 1, 10]
        weights = [1, 3, 1, 0, 100, 0, 1, 10]

        self.coords = zip(self.x_coords, self.y_coords)
        self.data = zip(self.x_coords, self.y_coords, values)
        self.weighted_data = zip(self.x_coords, self.y_coords, values, weights)

        bad_x_coords = [1, 1, 2, 4, 4, 5, 1, 2]
        bad_y_coords = [1, 4, 0, 5, 6, 1, 8, 6]
        self.bad_coords_A = zip(bad_x_coords, self.y_coords)
        self.bad_coords_B = zip(self.x_coords, bad_y_coords)
        self.bad_data_A = zip(bad_x_coords, self.y_coords, values)
        self.bad_data_B = zip(self.x_coords, bad_y_coords, values) 
        self.bad_wdata_A = zip(bad_x_coords, self.y_coords, values, weights)
        self.bad_wdata_B = zip(self.x_coords, bad_y_coords, values, weights) 

        self.bins = binner.Bins2D(int, 0, 0, 'custom', [1, 1.5, 4],
                                  int, 0, 0, 'custom', [1, 4.5, 5.5, 6.5, 10])

    def test_Count(self):
        binned_data = self.bins.bin_count(self.coords)
        expected_result = [[2,0,0,1],[2,1,2,0]]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_count, self.x_coords)
        self.assertRaises(binner.BinLimitError, self.bins.bin_count, self.bad_coords_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_count, self.bad_coords_B)

    def test_Count_with_data(self):
        binned_data = self.bins.bin_count(self.data)
        expected_result = [[2,0,0,1],[2,1,2,0]]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.BinLimitError, self.bins.bin_count, self.bad_data_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_count, self.bad_data_B)

    def test_Count_divide(self):
        binned_data = self.bins.bin_count_divide(self.coords)
        expected_result = [[2.0/4,0,0,1.0/4],
                           [2.0/12,1.0/3,2.0/3,0]]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_count_divide, self.x_coords)
        self.assertRaises(binner.BinLimitError, self.bins.bin_count_divide, self.bad_coords_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_count_divide, self.bad_coords_B)

    def test_Average(self):
        oldErrs = numpy.seterr(invalid='ignore')
        binned_data = self.bins.bin_average(self.data)
        numpy.seterr(**oldErrs)
        expected_result = [[2,None,None,1],[5,5,8,None]]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_average, self.x_coords)
        self.assertRaises(binner.DataTypeError, self.bins.bin_average, self.coords)
        self.assertRaises(binner.BinLimitError, self.bins.bin_average, self.bad_data_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_average, self.bad_data_B)

    def test_Average_Variance(self):
        oldErrs = numpy.seterr(invalid='ignore')
        binned_data = self.bins.bin_average(self.data, True)
        numpy.seterr(**oldErrs)
        expected_average = [[2,None,None,1],[5,5,8,None]]
        expected_variance = [[1.0,None,None,0.0],[4.0,0.0,4.0,None]]
        self.assertEqual(binned_data[0].tolist(), expected_average)
        self.assertEqual(binned_data[1].tolist(), expected_variance)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_average, self.x_coords)
        self.assertRaises(binner.DataTypeError, self.bins.bin_average, self.coords)
        self.assertRaises(binner.BinLimitError, self.bins.bin_average, self.bad_data_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_average, self.bad_data_B)

    def test_Sum(self):
        binned_data = self.bins.bin_sum(self.data)
        expected_result = [[4,None,None,1], [10,5,16,None]]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_sum, self.x_coords)
        self.assertRaises(binner.DataTypeError, self.bins.bin_sum, self.coords)
        self.assertRaises(binner.BinLimitError, self.bins.bin_sum, self.bad_data_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_sum, self.bad_data_B)

    def test_Sum_divide(self):
        binned_data = self.bins.bin_sum_divide(self.data)
        expected_result = [[4.0/4,None,None,1.0/4],
                           [10.0/12,5.0/3,16.0/3,None]]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_sum_divide, self.x_coords)
        self.assertRaises(binner.DataTypeError, self.bins.bin_sum_divide, self.coords)
        self.assertRaises(binner.BinLimitError, self.bins.bin_sum_divide, self.bad_data_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_sum_divide, self.bad_data_B)

    def test_Weighted_average(self):
        binned_data = self.bins.bin_weighted_average(self.weighted_data)
        expected_result = [[10./4,None,None,1],[7,0,700./110,None]]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_weighted_average, self.x_coords)
        self.assertRaises(binner.DataTypeError, self.bins.bin_weighted_average, self.coords)
        self.assertRaises(binner.DataTypeError, self.bins.bin_weighted_average, self.data)
        self.assertRaises(binner.BinLimitError, self.bins.bin_weighted_average, self.bad_wdata_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_weighted_average, self.bad_wdata_B)

    def test_Weighted_average_Variance(self):
        binned_data = self.bins.bin_weighted_average(self.weighted_data, True)
        expected_average = [[10./4,None,None,1],[7,0,700./110,None]]
        expected_variance = [[0.75,None,None,0.0],[0.0,0.0,4600/110.0-(700./110)**2,None]]
        self.assertEqual(binned_data[0].tolist(), expected_average)
        self.assertEqual(binned_data[1].tolist(), expected_variance)


    def test_Median(self):
        # to ignore python and numpy core warnings: 
        warnings.filterwarnings("ignore", category=RuntimeWarning) 
        binned_data = self.bins.bin_median(self.data)
        warnings.resetwarnings()
        expected_result = [[2,None,None,1],[5,5,8,None]]
        self.assertEqual(binned_data.tolist(), expected_result)

        # Check exceptions
        self.assertRaises(binner.DataTypeError, self.bins.bin_median, self.x_coords)
        self.assertRaises(binner.DataTypeError, self.bins.bin_median, self.coords)
        self.assertRaises(binner.BinLimitError, self.bins.bin_median, self.bad_data_A)
        self.assertRaises(binner.BinLimitError, self.bins.bin_median, self.bad_data_B)

    def test_normalize(self):
        """Function normalize(x)"""
        # Test with a list.
        input = [1,2,3,4]
        exp_output = [1./10, 2./10, 3./10, 4./10]
        act_output = binner.normalize(input)
        self.assert_(isinstance(act_output, type(input)))
        self.assertEqual(act_output, exp_output)

        # Test with a tuple.
        input = (1,2,3,4)
        exp_output = (1./10, 2./10, 3./10, 4./10)
        act_output = binner.normalize(input)
        self.assert_(isinstance(act_output, type(input)))
        self.assertEqual(act_output, exp_output)

        # Test with a numpy array.
        input = numpy.array([1,2,3,4])
        exp_output = numpy.array([1./10, 2./10, 3./10, 4./10])
        act_output = binner.normalize(input)
        self.assert_(isinstance(act_output, type(input)))
        self.assert_((act_output == exp_output).all())

        # Test with a masked array.
        x = numpy.array([1,2,3,4])
        input = numpy.ma.masked_array(x, x < 2.5)
        y = numpy.array([0, 0, 3./7, 4./7])
        exp_output = numpy.ma.masked_array(y, y < 2.5/7)
        act_output = binner.normalize(input)
        self.assert_(isinstance(act_output, type(input)))
        self.assert_(numpy.ma.equal(act_output, exp_output).all())
        self.assert_((act_output.mask == exp_output.mask).all())

if __name__ == '__main__':
    # The if-clause below exists only for debugging; it makes it
    # easier to run only one test instead of all of them. All test
    # will be run when binner.py is executed.
    if True:
        # Run only one test.
        suite = unittest.TestSuite()
        suite.addTest(TestBins("test_Average"))
        #suite.addTest(TestBins2D("test_Average"))
        unittest.TextTestRunner().run(suite)
    else:
        # Run all tests.
        unittest.main()
