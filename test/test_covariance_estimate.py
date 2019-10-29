import sys
sys.path.insert (0, '..')

import unittest
from covariance_estimate import calc_covariance_diagonal

class TestCovarianceMatrix (unittest.TestCase):

    def test_covariance (self):
        """ Tests if the module can correctly calculate the estimator
                Qn = sum_{i = 1}^n {(x_i - xbar)(x_i - xbar)^T}  / n """
        sample = [[1, 0, 1],
                  [0, 1, 0],
                  [1, 1, 0],
                  [0, 0, 1]]
        cov = calc_covariance_diagonal (sample)
        self.assertEqual (cov[0, 0], 1 / 4)
        self.assertEqual (cov[0, 1], 0)
        self.assertEqual (cov[0, 2], 0)
        self.assertEqual (cov[1, 0], 0)
        self.assertEqual (cov[1, 1], 1 / 4)
        self.assertEqual (cov[1, 2], -.25)
        self.assertEqual (cov[2, 0], 0)
        self.assertEqual (cov[2, 1], -.25)
        self.assertEqual (cov[2, 2], 1 / 4)

