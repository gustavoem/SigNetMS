import sys
sys.path.insert(0, '..')

import unittest
import numpy as np
import statistics 
from distributions.Gamma import Gamma

class TestGamma (unittest.TestCase):

    def test_get_mean (self):
        """ Tests if one can get the mean of the random variable. """
        a = 2
        b = 2
        X = Gamma (a, b)
        assert (abs (X.mean () - a * b) / (a * b) < 1e-2)


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
        N = 3000
        sample = [X.rvs () for _ in range (N)]
        mean = statistics.mean (sample)
        assert (abs (mean - 6) / 6 < 5e-2)



    def test_get_pdf (self):
        """ Tests if we can get a point of the probability density
            function of this random variable. """
        a = 2
        b = 1
        X = Gamma (a, b)
        x = 1
        analytic = np.exp (-1)
        assert (abs (X.pdf (x) - analytic) / analytic < 1e-2)
