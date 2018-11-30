import sys
sys.path.insert(0, '../src/distributions/')

import unittest
import numpy as np
from DiscreteLaplacian import DiscreteLaplacian

class TestDiscreteLaplacian (unittest.TestCase):

    def test_get_random_value (self):
        """ Tests if the distribution returns valid values. """
        rep = 100
        n = 11
        i = 4
        discrete_laplacian = DiscreteLaplacian (n, i)
        values = discrete_laplacian.rvs (rep) 
        assert all (np.array (values) > 0)
        assert all (np.array (values) <= n)


    def test_zero_probability_point (self):
        """ The pdf should be zero for the point i, and points smaller
            than 1 and points greater than n. """ 
        n = 13
        i = 4
        discrete_laplacian = DiscreteLaplacian (n, i)
        self.assertEqual (discrete_laplacian.pdf (0), 0)
        self.assertEqual (discrete_laplacian.pdf (i), 0)
        self.assertEqual (discrete_laplacian.pdf (n + 1), 0)
        self.assertEqual (discrete_laplacian.pdf (-1), 0)
        self.assertEqual (discrete_laplacian.pdf (n + 2), 0)


    def test_probability_density (self):
        """ A probability density function should sum 1 for all point
            desity. """
        n = 39
        i = 2
        discrete_laplacian = DiscreteLaplacian (n, i)
        p_sum = 0
        for j in range (1, n + 1):
            p_sum += discrete_laplacian.pdf (j)
        assert (abs (p_sum - 1) < 1e-4)
