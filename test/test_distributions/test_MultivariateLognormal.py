import sys
sys.path.insert(0, '..')

import unittest
import numpy as np
from distributions.MultivariateLognormal import MultivariateLognormal
from parallel_map import parallel_map

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


    def test_convergence_with_map (self):
        """ Tests if we can still achieve convergence when using 
            parallel_map. """
        mu = [2, 1]
        s = [[.01, 0], [0, .01]]
        X = MultivariateLognormal (mu, s)
        N = 2000
        rvs_caller = lambda _: X.rvs ()
        rand_vals = parallel_map (rvs_caller, range (N), 5)
        mean = np.array ([.0, .0])
        for val in rand_vals:
            mean += val
        mean /= N
        assert all (abs (X.mean () - mean) / X.mean () < 1e-1)
        
        mu = [2, 1]
        s = [[.01, 0], [0, .01]]
        X = MultivariateLognormal (mu, s)
        N = 2000

        rvs_caller = lambda _: X.rvs ()
        rand_vals = parallel_map (rvs_caller, range (N), 5)
        mean = np.array ([.0, .0])
        for val in rand_vals:
            mean += val
        mean /= N
        assert all (abs (X.mean () - mean) / X.mean () < 1e-1)
        
        # mu = [1e-12, 1]
        # s = [[10, 0], [0, 3]]
        # X = MultivariateLognormal.create_lognormal_with_shape (mu, s)
        # N = 80000
        # rvs_caller = lambda _: X.rvs ()
        # rand_vals = parallel_map (rvs_caller, range (N), 5)
        # mean = np.array ([.0, .0])
        # for val in rand_vals:
            # mean += val
        # mean /= N
        # print (mean)
        # assert all (abs (X.mean () - mean) / X.mean () < 1e-1)

    
    def test_convergence_to_mean (self):
        """ Tests if the randomly generated values has a sample mean
            that converges to the correct mean. """
        mu = [2, 1]
        s = [[.01, 0], [0, .01]]
        X = MultivariateLognormal (mu, s)
        N = 2000
        mean = np.array ([.0, .0])
        for i in range (N):
            mean += X.rvs ()
        mean /= N
        assert all (abs (X.mean () - mean) / X.mean () < 1e-1)
        
        mu = [2, 1]
        s = [[.01, 0], [0, .01]]
        X = MultivariateLognormal (mu, s)
        N = 2000
        mean = np.array ([.0, .0])
        for i in range (N):
            mean += X.rvs ()
        mean /= N
        assert all (abs (X.mean () - mean) / X.mean () < 1e-1)

        # mu = [1e-12, 1]
        # s = [[10, 0], [0, 3]]
        # X = MultivariateLognormal.create_lognormal_with_shape (mu, s)
        # N = 80000
        # mean = np.array ([.0, .0])
        # for i in range (N):
            # mean += X.rvs ()
        # mean /= N
        # assert all (abs (X.mean () - mean) / X.mean () < 1e-1)


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


    def test_get_pdf_of_zero_prob_point (self):
        """ Tests if we can get the pdf of a point with pdf equal to
            zero. """
        mu = [0, 0, 0]
        S = [[1, 0, 0], 
             [0, 1, 0], 
             [0, 0, 1]]
        X = MultivariateLognormal (mu, S)
        x = [1, 0, 1]
        analytic = 0
        assert (abs (X.pdf (x) - analytic) < 1e-4)


    def test_get_pdf_of_small_prob_point (self):
        """ Tests if we can get the pdf of a point with small pdf value.  
        """
        mu = [0, 0, 0]
        S = [[1, 0, 0], 
             [0, 1, 0], 
             [0, 0, 1]]
        X = MultivariateLognormal (mu, S)
        x = [1e-230, 1e-120, 1e-130]
        analytic = 0
        assert (abs (X.pdf (x) - analytic) < 1e-4)


    def test_create_shaped_distribution (self):
        """ Tests if we can create a Multivariate Lognormal 
            distribution with a specified mean and variance. """
        mu = [1, 2, .1]
        S = np.array ([[.1,  0,  0], 
                       [ 0, .2,  0], 
                       [ 0,  0, .2]])
        X = MultivariateLognormal.create_lognormal_with_shape (mu, S)
        N = 2000
        mean = np.array ([.0, .0, .0])
        for i in range (N):
            x = X.rvs ()
            mean += x
        mean /= N
        assert all (abs (mu - mean) / mean < 2e-1)
        
