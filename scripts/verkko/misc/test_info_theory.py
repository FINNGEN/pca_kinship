import unittest
import numpy as np
import info_theory as it


# All verification results have been calculated with Matlab library
# "Mutual Information computation" by Hanchuan Peng
# (http://www.mathworks.co.uk/matlabcentral/fileexchange/14888)

class Test_IT(unittest.TestCase):
    def setUp(self):
        # First data (in format suitable for the Matlab library)
        a = [1, 2, 1, 2, 1]
        b = [2, 1, 2, 1, 1]
        c = [2, 1, 2, 2, 1]
        # Construct probability matrix
        self.P = np.zeros((max(a)-min(a)+1, \
                           max(b)-min(b)+1, \
                           max(c)-min(c)+1), float)
        for i,j,k in zip(a,b,c):
            self.P[(i-1,j-1,k-1)] += 1
        self.P = self.P/np.sum(self.P)

        # Second data 
        a = [3, 1, 2, 3, 1, 2, 1, 3, 2, 3]
        b = [2, 3, 1, 2, 3, 1, 1, 3, 2, 3]
        c = [3, 2, 1, 2, 2, 1, 3, 1, 3, 1]
        # Construct probability matrix
        self.Q = np.zeros((max(a)-min(a)+1, \
                           max(b)-min(b)+1, \
                           max(c)-min(c)+1), float)
        for i,j,k in zip(a,b,c):
            self.Q[(i-1,j-1,k-1)] += 1
        self.Q = self.Q/np.sum(self.Q)

    def test_entropy_X_Y(self):
        old_errs = np.seterr(under="ignore", divide="ignore")
        ans_str = "%.10f" % it.entropy_X_Y(np.sum(self.P, 2))
        corr_ans = "%.10f" % 0.550977500432694
        self.assertEqual(ans_str, corr_ans)

        ans_str = "%.10f" % it.entropy_X_Y(np.sum(self.P, 1))
        corr_ans = "%.10f" % 0.950977500432694
        self.assertEqual(ans_str, corr_ans)

        ans_str = "%.10f" % it.entropy_X_Y(np.sum(self.P, 0))
        corr_ans = "%.10f" % 0.550977500432694
        self.assertEqual(ans_str, corr_ans)
        np.seterr(**old_errs)

    def test_mutual_info_XY_Y(self):
        old_errs = np.seterr(under="ignore", divide="ignore")
        ans_str = "%.10f" % it.mutual_info_XY_Z(self.P)
        corr_ans = "%.10f" % 0.550977500432694
        self.assertEqual(ans_str, corr_ans)

        ans_str = "%.10f" % it.mutual_info_XY_Z(self.Q)
        corr_ans = "%.10f" % 0.950977500432694
        self.assertEqual(ans_str, corr_ans)
        np.seterr(**old_errs)

if __name__ == '__main__':
    #suite = unittest.TestSuite()
    #suite.addTest(Test_IT("test_entropy_X_Y"))
    #unittest.TextTestRunner().run(suite)
    unittest.main()
