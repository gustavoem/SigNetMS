import sys
sys.path.insert(0, '../src/distributions/')

import unittest
import numpy as np
from Lognormal import Lognormal

class TestLognormal (unittest.TestCase):

    def test_get_mean (self):
        """ Tests if one can get the mean of the random variable. """
        mu = 2
        s = 2
        X = Lognormal (mu, s)
        assert (abs (X.mean () - np.exp (mu + s * s / 2)) < 1e-4)


    def test_get_random_value (self):
        """ Tests if one can get the a random value from the random
            variable. """
        mu = -100
        s = 1
        X = Lognormal (mu, s)
        assert (X.rvs () > 0)

    
    def test_convergence_to_mean (self):
        """ Tests if random values has mean as Lognormal should have. 
            WARNING: this test can fail with some small probability. """
        normal_mu = 3
        normal_s2 = .001
        # convergence to mean is very slow if the normal distribution
        # has large variance (which means even larger variance of the
        # lognormal)
        X = Lognormal (normal_mu, normal_s2)
        N = 100000
        mean = 0
        for i in range (N):
            mean += X.rvs ()
        mean /= N
        analytical_mean = np.exp (normal_mu + normal_s2 / 2)
        assert (abs (analytical_mean - mean) < 1e-1)


    def test_get_pdf (self):
        """ Tests if we can get a point of the probability density
            function of this random variable. """
        mu = 1
        s = 1
        X = Lognormal (mu, s)
        x = 2
        # If the underlying Normal distribution is N (mu, s²) then,
        # the pdf of X is f(x):
        # 1/sqrt (2pi * s²) exp (-(ln(x) - mu)²/(2 * s²)) * 1/x
        analytic = (1 / np.sqrt (2 * np.pi * s * s)) \
                    * np.exp (-(np.log (x) - mu) ** 2 / (2 * s * s)) \
                    * (1 / x)
        assert (abs (X.pdf (x) - analytic) < 1e-4)
