import sys
sys.path.insert(0, '..')

import unittest
import numpy as np
from distributions.MultivariateNormal import MultivariateNormal

class TestMultivariateLognormal (unittest.TestCase):

    def test_get_mean (self):
        """ Tests if one can get the mean of the random variable. """
        mu = [1, 2, 3, 4]
        S = np.array ([[2, 0, 0, 0], 
             [0, 1, 0, 0], 
             [0, 0, 4, 1],
             [0, 0, 1, 2]])
        X = MultivariateNormal (mu, S)
        assert all (abs (X.mean () - np.array (mu)) < 1e-4)


    def test_convergence_to_mean (self):
        """ Tests if the randomly generated values has a sample mean
            that converges to the correct mean. """
        mu = [2, 1]
        s = [[.01, 0], [0, .01]]
        X = MultivariateNormal (mu, s)
        N = 1000
        mean = np.array ([.0, .0])
        for i in range (N):
            mean += X.rvs ()
        mean /= N
        assert all (abs (X.mean () - mean) < 1e-1)
        
        mu = [2, 1]
        s = [[.01, -.005], [-.005, .01]]
        X = MultivariateNormal (mu, s)
        N = 1000
        mean = np.array ([.0, .0])
        for i in range (N):
            mean += X.rvs ()
        mean /= N
        assert all (abs (X.mean () - mean) < 1e-1)



    def test_get_pdf (self):
        """ Tests if we can get a point of the probability density
            function of this random variable. """
        mu = [0, 0, 0]
        S = [[1, 0, 0], 
             [0, 1, 0], 
             [0, 0, 1]]
        X = MultivariateNormal (mu, S)
        x = [1, 1, 1]
        analytic = 1 / np.sqrt ((2 * np.pi) ** 3) * np.exp (-.5 * 3)
        assert (abs (X.pdf (x) - analytic) < 1e-4)


    def test_copy (self):
        """ Tests if we can copy an object. """
        mu = [1, 2, 3, 4]
        S = np.array ([[2, 0, 0, 0], 
             [0, 1, 0, 0], 
             [0, 0, 4, 1],
             [0, 0, 1, 2]])
        X = MultivariateNormal (mu, S)
        Y = X.copy ()
        x = [1, 1, 1, 1]
        assert all (abs (Y.mean () - X.mean ()) < 1e-4)
        assert (abs (X.pdf (x) - Y.pdf (x)) < 1e-4)

