import sys
sys.path.insert(0, '../src/distributions/')

import unittest
import numpy as np
from MultivariateLognormal import MultivariateLognormal

class TestMultivariateLognormal (unittest.TestCase):

    def test_get_mean (self):
        """ Tests if one can get the mean of the random variable. """
        mu = [1, 2, 3, 4]
        S = np.array ([[2, 0, 0, 0], 
             [0, 1, 0, 0], 
             [0, 0, 4, 1],
             [0, 0, 1, 2]])
        X = MultivariateLognormal (mu, S)
        assert all (abs (X.mean () - np.exp (np.array (mu) \
                + S.diagonal () / 2)) < 1e-4)


    def test_get_random_value (self):
        """ Tests if one can get the a random value from the random
            variable. """
        mu = [2, 1]
        s = [[1, 0], [0, 1]]
        X = MultivariateLognormal (mu, s)
        X_0 = X.rvs ()
        assert all (X_0 > 0)


    def test_get_pdf (self):
        """ Tests if we can get a point of the probability density
            function of this random variable. """
        mu = [0, 0, 0]
        S = [[1, 0, 0], 
             [0, 1, 0], 
             [0, 0, 1]]
        X = MultivariateLognormal (mu, S)
        x = [1, 1, 1]
        analytic = 1 / np.sqrt ((2 * np.pi) ** 3 )
        assert (abs (X.pdf (x) - analytic) < 1e-4)
