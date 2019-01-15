import sys
sys.path.insert (0, '../src/')
sys.path.insert (0, '../src/samplers/')

import unittest
import numpy as np
from MetropolisHastings import MetropolisHastings
from RandomParameterList import RandomParameterList
from RandomParameter import RandomParameter
from Gamma import Gamma
from MultivariateLognormal import MultivariateLognormal


class MockMH (MetropolisHastings):

    def _create_jump_dist (self, theta_t):
        n = theta_t.get_size ()
        s2 = np.eye (n) / 100
        variances = s2.diagonal ()
        theta_values = np.array (theta_t.get_values ())
        mu = np.log (theta_values) - variances / 2
        return MultivariateLognormal (mu, s2)

class TestMetropolisHastings (unittest.TestCase):

    def test_jumps_centered_on_current_theta (self):
        n = 10
        N = 1000
        theta_values = np.ones (n) / 2
        theta = RandomParameterList ()
        for p_val in theta_values:
            gamma = Gamma (2, 1)
            rand_par = RandomParameter ('p', gamma)
            rand_par.value = p_val
            theta.append (rand_par)

        mocked_mh = MockMH (theta)
        mean_jump = np.zeros (n)
        for i in range (N):
            jump = mocked_mh.propose_jump (theta)
            mean_jump += jump.get_values ()
        mean_jump /= N
        print (mean_jump)
        assert (False)
