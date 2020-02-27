import bootstrap as bs
import numpy as np

import unittest


class Test(unittest.TestCase):

    def test_boostrap_conf_intervals(self):
        data = np.zeros(100)
        val1, val2 = bs.mean_conf_interval(data, 100)
        assert val1 == 0.0
        assert val2 == 0.0

        data = np.arange(10)
        percentiles = bs.mean_conf_interval(data, 100)
        assert percentiles[0] > 0.1
        assert percentiles[1] < 8.9

        # if coverage == 0, the percentiles should be the same:
        percentiles = bs.mean_conf_interval(data, 100, 0)
        assert percentiles[0] == percentiles[1]

        # if coverage = 100, the min and max values should be obtained:
        # (if enough iterations)
        data = np.arange(2)
        percentiles = bs.mean_conf_interval(data, 1000, 100)
        # print probability of this test failing is pretty low
        assert percentiles[0] == np.min(data)
        assert percentiles[1] == np.max(data)

        # simple test with two-d data:
        data = np.array([[1, 2, 1, 2],
                         [2, 3, 2, 3],
                         [3, 4, 3, 4]])
        percentiles = bs.mean_conf_interval(data, 1000, 100)
        print percentiles.shape == (3, 2)
        assert (percentiles[0] < percentiles[1]).all()
        assert percentiles[0, 0] < 1.3
        assert percentiles[2, 1] > 3.7

    def test_mean_groupwise_conf_intervals_from_sim_matrix(self):
        mat = [[0, 1, 1, 0, 0, 0, 0],
              [0, 0, 1, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 2, 2, 2],
              [0, 0, 0, 0, 0, 2, 2],
              [0, 0, 0, 0, 0, 0, 2],
              [0, 0, 0, 0, 0, 0, 0]]
        mat = np.array(mat)
        means, intervals = bs.mean_groupwise_conf_intervals_from_sim_matrix(
            mat, 3, 10, 100)
        assert (means == np.array([1, 2])).all()
        assert intervals[0][0] == 1
        assert intervals[0][1] == 1
        assert intervals[1][0] == 2
        assert intervals[1][1] == 2

        # twist the matrix values a little:
        mat = [[0, 1, 0.5, 0, 0, 0, 0],
              [0, 0, 1.5, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 2, 2, 2.5],
              [0, 0, 0, 0, 0, 1.5, 2],
              [0, 0, 0, 0, 0, 0, 2],
              [0, 0, 0, 0, 0, 0, 0]]
        mat = np.array(mat)
        means, intervals = bs.mean_groupwise_conf_intervals_from_sim_matrix(
            mat, 3, 100, 75)
        assert (means == np.array([1, 2])).all()
        assert intervals[0][0] < 1
        assert intervals[0][1] > 1
        assert intervals[1][0] < 2
        assert intervals[1][1] > 2
