import sys
sys.path.insert(0, '../src/distributions/')

import unittest
import numpy as np
from Gamma import Gamma

class TestGamma (unittest.TestCase):

    def test_get_mean (self):
        """ Tests if one can get the mean of the random variable. """
        a = 2
        b = 2
        X = Gamma (a, b)
        assert (abs (X.mean () - a * b) < 1e-4)


    def test_get_random_value (self):
        """ Tests if one can get the a random value from the random
            variable. """
        a = 2
        b = 100
        X = Gamma (a, b)
        assert all (np.array (X.rvs (100)) > 0)

    def test_convergence_to_mean (self):
        """ Tests if random values has mean as Gamma should have. 
            WARNING: this test can fail with some small probability. """
        a = 2
        b = 3
        X = Gamma (a, b)
        N = 10000
        mean = 0
        for i in range (N):
            mean += X.rvs ()
        mean /= N
        assert (abs (mean - 6) < 1e-1)



    def test_get_pdf (self):
        """ Tests if we can get a point of the probability density
            function of this random variable. """
        a = 2
        b = 1
        X = Gamma (a, b)
        x = 1
        analytic = np.exp (-1)
        assert (abs (X.pdf (x) - analytic) < 1e-4)
