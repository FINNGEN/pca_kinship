import numpy as np
import parallel_run_helper as prh
import unittest


def worker_multiplier(vals):
    return vals[0] * vals[1]


class Test(unittest.TestCase):

    def setUp(self):
        pass

    def test_run_in_parallel(self):
        arg_list = zip([1, 2, 3], [0, 2, 1])
        answers = np.array([0, 4, 3])
        work_func = worker_multiplier

        for n_cpus in [1, 2]:
            comps = prh.run_in_parallel(work_func, arg_list, n_cpus)
            comps = np.array(comps)
            self.assertTrue(np.all(answers == comps))
