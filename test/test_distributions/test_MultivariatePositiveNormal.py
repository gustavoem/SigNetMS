import sys
sys.path.insert(0, '..')

import unittest
import numpy as np
from distributions.MultivariatePositiveNormal import MultivariatePositiveNormal

class TestMultivariateLognormal (unittest.TestCase):

    def test_get_pdf (self):
        """ Tests if we can get a point of the probability density
            function of this random variable. """
        mu = [1, 1, 1]
        S = [[  1, .1, -.3], 
             [ .1,  1,   0], 
             [-.3,  0,   1]]
        X = MultivariatePositiveNormal (mu, S)
        x = [1.5, 1.5, 1.5]
        y = [0.5, 0.5, 0.5]

        assert (abs (X.pdf_ratio (x, y) - 1) < 1e-4)

    def test_is_simmetric (self):
        """ Tests if the distribution is symmetric. """
        mu1 = [1, 1, 1]
        mu2 = [3, 3, 3]
        S = [[  2, .1, -.3], 
             [ .1,  2,   0], 
             [-.3,  0,   1]]
        X1 = MultivariatePositiveNormal (mu1, S)
        X2 = MultivariatePositiveNormal (mu2, S)
        p1 = X1.log_pdf (mu2)
        p2 = X2.log_pdf (mu1)
        assert (abs (p1 - p2) < 1e-4)

