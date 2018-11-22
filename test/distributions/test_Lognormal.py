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
        self.assertEqual (X.mean (), np.exp (mu + s * s / 2))

    def test_fail (self):
        self.assertEqual (1, 2)
