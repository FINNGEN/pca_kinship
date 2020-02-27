import permute
import unittest
import numpy as np


class TestPermutation(unittest.TestCase):

    def setUp(self):
        self.as_arr_alm_eq = np.testing.assert_almost_equal
        self.data = np.array([1., 2, 3, 4, 5, 6])
        self.perm1 = [0, 1, 2, 3, 4, 5]
        self.perm2 = [1, 2, 3, 4, 5, 0]

        self.two_dim_data = np.array([
                                     [1, 1],
                                     [2, 2],
                                     [3, 3],
                                     [4, 4],
                                     [5, 5],
                                     [6, 6]
                                     ])

        self.matrix = np.array([[11, 12, 13, 14],
                                [21, 22, 23, 24],
                                [31, 32, 33, 34],
                                [41, 42, 43, 44]
                                ])
        self.mat_perm_1 = np.array([3, 1, 2, 0])
        self.mat_perm_2 = np.array([1, 2, 3, 0])
        self.rng = permute.get_random_state(seed=123456)

    def test_permute_array(self):
        newdata1 = permute.permute_array(self.data, self.perm1)
        self.as_arr_alm_eq(newdata1, np.array([1, 2, 3, 4, 5, 6]))
        newdata2 = permute.permute_array(self.data, self.perm2)
        self.as_arr_alm_eq(newdata2, np.array([2, 3, 4, 5, 6, 1]))

        two_dim_data_1 = permute.permute_array(self.two_dim_data, self.perm1)
        result_should_be = np.array([
                                    [1, 1],
                                    [2, 2],
                                    [3, 3],
                                    [4, 4],
                                    [5, 5],
                                    [6, 6]
                                    ])
        self.as_arr_alm_eq(two_dim_data_1, result_should_be)

        two_dim_data_2 = permute.permute_array(self.two_dim_data, self.perm2)
        result_should_be = np.array([
                                    [2, 2],
                                    [3, 3],
                                    [4, 4],
                                    [5, 5],
                                    [6, 6],
                                    [1, 1]
                                    ])
        self.as_arr_alm_eq(two_dim_data_2, result_should_be)

    def test_permute_matrix(self):
        perm_mat = permute.permute_matrix(self.matrix, self.mat_perm_1)
        result_should_be = np.matrix([
            [44, 42, 43, 41],
            [24, 22, 23, 21],
            [34, 32, 33, 31],
            [14, 12, 13, 11]]
        )
        self.as_arr_alm_eq(perm_mat, result_should_be)

        perm_mat_2 = permute.permute_matrix(self.matrix, self.mat_perm_2)
        result_should_be = np.matrix([
            [22, 23, 24, 21],
            [32, 33, 34, 31],
            [42, 43, 44, 41],
            [12, 13, 14, 11]]
        )
        self.as_arr_alm_eq(perm_mat_2, result_should_be)

    def test_half_permute_paired_matrix(self):

        initmat = np.matrix([
            [100, 1, 10, 0],
            [1, 100, 0, 10],
            [10, 0, 900, 2],
            [0, 10, 2, 800]
        ])
        half_perm_mat = permute.half_permute_paired_matrix(initmat,
                                                           self.mat_perm_1
                                                           )
        result_should_be = np.matrix([
            [100, 1, 0, 10],
            [1, 100, 10, 0],
            [0, 10, 800, 2],
            [10, 0, 2, 900]
        ])
        self.as_arr_alm_eq(half_perm_mat, result_should_be)
        result_should_be = \
            initmat[np.array([0, 1, 3, 2]), :][:, np.array([0, 1, 3, 2])]
        self.as_arr_alm_eq(half_perm_mat, result_should_be)

    def assert_paired_permutation_validity(self, perm, n1):
        """
        Test that the permutation is a *paired* permutation
        """
        index_sums = []
        for i in range(n1):
            index_sums.append(perm[i]+perm[i+n1])
        # one requirement for a paired permutation:
        for i in range(n1-1):
            assert index_sums[i] == index_sums[i+1] - 2
        # more explicit requirement:
        for i in range(2*n1):
            pair = perm[i]
            assert perm[pair] == i
        self.assert_permutation_validity(perm, 2*n1)

    def assert_permutation_validity(self, perm, n):
        self.as_arr_alm_eq(np.sort(np.unique(perm)), np.arange(n))

    def test_get_random_paired_permutation(self):
        paired = True
        n1 = n2 = 4
        #to test for different permutations
        for _ in range(10):
            perm = permute.get_permutation(paired, n1, n2, self.rng)
            self.assert_paired_permutation_validity(perm, n1)

    def test_get_paired_permutation(self):
        paired = True
        perms = []
        n1 = n2 = 4
        for i in range(2**n1):
            perm = permute.get_permutation(paired, n1, n2, None, i)
            self.assert_paired_permutation_validity(perm, n1)
            for other_perm in perms:
                assert not (other_perm == perm).all()
            perms.append(perm)

    def test_get_normal_permutation(self):
        paired = False
        n1 = n2 = 4
        for i in range(10):
            p = permute.get_permutation(paired, n1, n2, rng=self.rng)
            self.assert_permutation_validity(p, n1+n2)
