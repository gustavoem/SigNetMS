import sys
sys.path.insert (0, '..')

import os
import unittest
import numpy as np
from marginal_likelihood.PriorsReader import read_priors_file
from marginal_likelihood.RandomParameter import RandomParameter
from distributions.Gamma import Gamma


class TestPriors (unittest.TestCase):

    def is_nan (self, x):
        return x != x

    def test_get_priors_log_prob (self):
        """ Given theta, we should be able to get
            log {p(theta)}. """
        priors_file = 'input/kolch_model.priors'  
        theta = read_priors_file (priors_file)
        # theta has been randomly generated by the joint distribution,
        # and there's a high (but of course not 1) probability of the
        # log of the pdf being representable in floating point (it 
        # could not be if the value is really small).
        log_prior = theta.log_p ()
        assert (not self.is_nan (log_prior))
        assert (log_prior > float ('-inf'))


    def test_get_param_log_p (self):
        """ Tests if the log p is what we expect when analytically
            calculating it. """
        param = RandomParameter ("", Gamma (2, 1))
        param.value = 5
        log_p = param.get_log_p ()
        analytic = np.log (param.value) - param.value
        assert (abs (log_p - analytic) < 1e-8)
