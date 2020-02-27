import measures
import numpy as np
import unittest
from nose.tools import assert_almost_equal


class TestMeasures(unittest.TestCase):

    def setUp(self):
        self.as_arr_alm_eq = np.testing.assert_almost_equal
        self.one_dim_data_simple = np.array(
            [0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3]
        )
        self.one_dim_data_paired = np.array(
            [0, 0, 0, 1, 1, 1, 3, 3, 3, 2, 2, 2]
        )
        self.meandiff = -2
        self.n1 = len(self.one_dim_data_simple) / 2
        self.sim_mat_1 = np.array(
            [[1, 2, 3, 4],
             [2, 3, 4, 5],
             [3, 4, 5, 6],
             [4, 5, 6, 7]]
        )
        self.sim_mat_1_n1 = 2
        self.sim_mat_2 = np.array(
            [[1, 2, 3, 4, 5],
             [2, 3, 4, 5, 6],
             [3, 4, 5, 6, 7],
             [4, 5, 6, 7, 8],
             [5, 6, 7, 8, 9]]
        )
        self.sim_mat_2_n1 = 2
        self.sim_mat_3 = np.array(
            [[1, 2, 3, 10, 5, 6],
             [2, 3, 4, 5, 10, 7],
             [3, 4, 5, 6, 7, 10],
             [10, 5, 6, 7, 8, 9],
             [5, 10, 7, 8, 9, 10],
             [6, 7, 10, 9, 10, 11]]
        )
        self.sim_mat_3_n1 = 3

    def test_paired_t_value(self):
        """
        See documentation from here:
        http://en.wikipedia.org/wiki/Student's_t-test
        for the t-values
        """
        d_sdevs = 1 ** 2 * np.sqrt(self.n1 / float(self.n1 - 1))
        paired_t_val_should_be = -2 / (d_sdevs / np.sqrt(self.n1))
        computed_paired_t_val = measures.paired_t_value(
            self.one_dim_data_paired, self.n1)
        assert_almost_equal(paired_t_val_should_be, computed_paired_t_val)

    def test_mean_difference(self):
        md1 = measures.mean_difference(self.one_dim_data_paired, self.n1)
        md2 = measures.mean_difference(self.one_dim_data_simple, self.n1)
        assert_almost_equal(self.meandiff, md1)
        assert_almost_equal(md1, md2)

    def test_two_dim_paired_t_value(self):
        self.two_dim_helper(measures.paired_t_value)

    def test_two_dim_t_value(self):
        self.two_dim_helper(measures.unpaired_t_value)

    def test_two_dim_mean_difference(self):
        self.two_dim_helper(measures.mean_difference)

    def two_dim_helper(self, func):
        one_dim_val = func(self.one_dim_data_paired, self.n1)
        true_vals = []
        two_dim_data = []
        for i in range(3 * self.n1):  # arbitrary no of loops
            two_dim_data.append(self.one_dim_data_paired)
            true_vals.append(one_dim_val)

        compt_vals = func(np.array(two_dim_data).T, self.n1)
        self.as_arr_alm_eq(compt_vals, np.array(true_vals))

    def test_sim_matrix_within_groups_means(self):
        group_means_1 = measures.sim_matrix_within_group_means(
            self.sim_mat_1, self.sim_mat_1_n1)
        self.as_arr_alm_eq(
            np.array(group_means_1), np.array([2., 6, -4]))
        group_means_2 = measures.sim_matrix_within_group_means(
            self.sim_mat_2, self.sim_mat_2_n1)
        self.as_arr_alm_eq(
            np.array(group_means_2), np.array([2., 7, -5]))
        group_means_2b = measures.sim_matrix_within_group_means(
            self.sim_mat_2, len(self.sim_mat_2) - self.sim_mat_2_n1)
        self.as_arr_alm_eq(
            np.array(group_means_2b), np.array([3., 8, -5]))
        group_means_3 = measures.sim_matrix_within_group_means(
            self.sim_mat_3, self.sim_mat_3_n1)
        self.as_arr_alm_eq(
            np.array(group_means_3), np.array([3., 9, -6])
        )

    def test_paired_sim_matrix_inter_group_means(self):
        inter_group_means_paired = \
            measures.paired_sim_matrix_inter_group_means(self.sim_mat_3)
        self.as_arr_alm_eq(
            np.array(inter_group_means_paired), np.array([6., 10, 4.])
        )

    def test_sim_matrix_within_group_means_minus_inter_group_mean(self):
        wmi = measures.sim_matrix_within_groups_mean_minus_inter_group_mean(
            self.sim_mat_3, True
        )
        assert_almost_equal(wmi, 0)
        wmi = measures.sim_matrix_within_groups_mean_minus_inter_group_mean(
            self.sim_mat_3, False, n1=3

        )
        assert_almost_equal(wmi, (3 + 9) / 2. - 66 / 9.)

    def test_sim_matrix_mean_inter_group_similarity(self):
        val = measures.sim_matrix_mean_inter_group_similarity(self.sim_mat_2,
                                                              self.sim_mat_2_n1
                                                              )
        assert_almost_equal(val, 27/6.)
        val = measures.sim_matrix_mean_inter_group_similarity(self.sim_mat_3,
                                                              self.sim_mat_3_n1
                                                              )
        assert_almost_equal(val, 66/9.)
