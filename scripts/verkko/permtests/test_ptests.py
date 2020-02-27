import ptests
import measures
import unittest
import numpy as np


class TestPtests(unittest.TestCase):

    def setUp(self):
        self.as_arr_alm_eq = np.testing.assert_almost_equal
        self.simple_data_n1 = 3
        self.simple_data_n2 = 3
        self.simple_data = np.array([-2, 1, 2,
                                     2, 3, 6])
        # multiply by two as we report two-sided p-values
        self.simple_data_p_paired = 1./2**3*2

        self.simple_data_p_unpaired = 2./(6*5*4/(3*2*1))*2
        # 2 (number of combinations for such a low value for md or t-value)
        # /(6*5*4/(3*2*1)) tot number of combinations
        # *2 two-sided p-value

        self.simple_data_meandiff = -10/3.
        self.simple_data_t_value_paired = measures.paired_t_value(
            self.simple_data,
            self.simple_data_n1
        )
        self.simple_data_t_value_unpaired = measures.unpaired_t_value(
            self.simple_data,
            self.simple_data_n1
        )

        # self.sim_mat = ...

    def test_simple_one_dim(self):
        # Testing t-value and meandiff
        # First with paired setup
        paired_study = True
        stats = ptests.mean_difference_permtest(self.simple_data,
                                               self.simple_data_n1,
                                               self.simple_data_n2,
                                               paired_study,
                                               'all',
                                               seed=1,
                                               n_cpus=1
                                               )
        self.assertEqual(stats[0], self.simple_data_p_paired)
        self.assertAlmostEqual(stats[1], self.simple_data_meandiff)
        stats = ptests.t_value_permtest(self.simple_data,
                                       self.simple_data_n1,
                                       self.simple_data_n2,
                                       paired_study,
                                       'all',
                                       seed=1,
                                       n_cpus=1
                                       )
        self.assertEqual(stats[0], self.simple_data_p_paired)
        self.assertAlmostEqual(stats[1], self.simple_data_t_value_paired)

        # Same stuff with a unpaired setup
        paired_study = False
        stats = ptests.mean_difference_permtest(self.simple_data,
                                               self.simple_data_n1,
                                               self.simple_data_n2,
                                               paired_study,
                                               10**4.,
                                               n_cpus=1,
                                               seed=10)
        self.assertTrue(np.abs(stats[0]-self.simple_data_p_unpaired) < 0.01)
        self.assertAlmostEqual(stats[1], self.simple_data_meandiff)

        stats = ptests.t_value_permtest(self.simple_data,
                                       self.simple_data_n1,
                                       self.simple_data_n2,
                                       paired_study,
                                       10**4.,
                                       seed=1,
                                       n_cpus=1
                                       )
        self.assertTrue(np.abs(stats[0]-self.simple_data_p_unpaired) < 0.01)
        self.assertAlmostEqual(stats[1], self.simple_data_t_value_unpaired)

    def test_two_dim_and_parallel(self):
        #test two dim data with just the basic mean difference
        sdn = 10
        two_dim_data = np.array([self.simple_data for i in range(sdn)]).T
        #note the transpose!!!

        two_dim_pvals = np.ones((sdn, ))*self.simple_data_p_paired
        two_dim_meandiffs = np.ones((sdn, ))*self.simple_data_meandiff
        paired_study = True
        stats = ptests.mean_difference_permtest(two_dim_data,
                                               self.simple_data_n1,
                                               self.simple_data_n2,
                                               paired_study,
                                               'all',
                                               seed=1,
                                               n_cpus=1
                                               )
        self.as_arr_alm_eq(stats[0], two_dim_pvals)
        self.as_arr_alm_eq(stats[1], two_dim_meandiffs)

        stats = ptests.mean_difference_permtest(two_dim_data,
                                               self.simple_data_n1,
                                               self.simple_data_n2,
                                               paired_study,
                                               'all',
                                               seed=1,
                                               n_cpus=5
                                               )
        self.as_arr_alm_eq(stats[0], two_dim_pvals)
        self.as_arr_alm_eq(stats[1], two_dim_meandiffs)

    def test_sim_mat_permtests(self):
        """Testing permtests for sim matrices"""
        sim_mat = np.array([
            [10, 1, 1, 5, 1.5, 1.5],
            [1, 10, 1, 1.5, 5, 1.5],
            [1, 1, 10, 1.5, 1.5, 5],
            [5, 1.5, 1.5, 20, 2, 2],
            [1.5, 5, 1.5, 2, 20, 2],
            [1.5, 1.5, 5, 2, 2, 20]
        ])
        assert (sim_mat == sim_mat.T).all()
        n1 = n2 = 3
        paired_p_min_within_group_means = 1./2**3*2
        mean_diff_within_group_diff = -1.  # only upper triangles
        paired_study = True
        # is the first group more consistent than the other?
        stats = ptests.sim_matrix_within_group_mean_diff_permtests(
            sim_mat,
            n1,
            n2,
            paired_study,
            'all',
            seed=0
        )

        self.assertAlmostEqual(paired_p_min_within_group_means, stats[0][2])
        self.assertAlmostEqual(mean_diff_within_group_diff, stats[1][2])

        unpaired_p_min_within_group_means = 1./(6*5*4/(3*2)) * 2
        paired_study = False
        stats = ptests.sim_matrix_within_group_mean_diff_permtests(
            sim_mat,
            n1,
            n2,
            paired_study,
            10**4,
            seed=0
        )

        self.assertTrue(
            np.abs(stats[0][2]-unpaired_p_min_within_group_means) < 0.01
        )
        self.as_arr_alm_eq(stats[1][2], mean_diff_within_group_diff)

        sim_mat = np.array([
            [10, 1, 1, 5, 0, 0],
            [1, 10, 1, 0, 5, 0],
            [1, 1, 10, 0, 0, 5],
            [5, 0, 0, 20, 2, 2],
            [0, 5, 0, 2, 20, 2],
            [0, 0, 5, 2, 2, 20]
        ])
        assert (sim_mat == sim_mat.T).all()
        n1 = n2 = 3

        #def sim_matrix_inter_group_means_permtest(matrix, nIt=1e6, seed=None):

        distance_mat = np.array([
            [0, 0, 0, 0, 1, 1, 1, 1],
            [0, 0, 0, 0, 1, 1, 1, 1],
            [0, 0, 0, 0, 1, 1, 1, 1],
            [0, 0, 0, 0, 1, 1, 0.9, 1.1],
            [1, 1, 1, 1, 0, 0, 0, 0],
            [1, 1, 1, 1, 0, 0, 0, 0],
            [1, 1, 1, 0.9, 0, 0, 0, 0],
            [1, 1, 1, 1.1, 0, 0, 0, 0]
        ])
        paired_study = True
        paired_p_val = 2./(2**4)*2
        p_val, inter_mean = \
            ptests.sim_matrix_group_distance_permtest(distance_mat,
                                                     4,
                                                     4,
                                                     paired_study,
                                                     'all',
                                                     seed=5)
        self.assertTrue(np.abs(p_val-paired_p_val) < 0.01)
        self.as_arr_alm_eq(inter_mean, 1)

        distance_mat = np.array([
            [0, 0, 0, 0, 3, 1, 1, 1],
            [0, 0, 0, 0, 1, 3, 1, 1],
            [0, 0, 0, 0, 1, 1, 3, 1],
            [0, 0, 0, 0, 1, 1, 1, 3],
            [3, 1, 1, 1, 0, 0, 0, 0],
            [1, 3, 1, 1, 0, 0, 0, 0],
            [1, 1, 3, 1, 0, 0, 0, 0],
            [1, 1, 1, 3, 0, 0, 0, 0]
        ])

        #test pairedness
        func_to_test = ptests.sim_matrix_semidiag_vs_inter_group_permtest
        pvals, dists = func_to_test(distance_mat, 10**4, seed=10)
        self.assertEqual(dists[0], 1)  # inter_mean (without semidiag)
        self.assertEqual(dists[1], 3)  # semidiag mean
        self.assertEqual(dists[2], 3-1)  # semidiag - inter_mean
        pval_should_be = 1.0/(4*3*2*1)*2
        self.assertTrue(np.abs(pvals[-1] - pval_should_be) < 0.01)
