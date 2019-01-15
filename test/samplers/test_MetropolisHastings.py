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



class TestMetropolisHastings (unittest.TestCase):

    def test_jumps_centered_on_current_theta (self):
        """ Tests if, when using suitable jump distribution, the 
            proposed values are centered on the current theta. """
        n = 10
        N = 1000
        theta_values = np.ones (n) / 3.2
        theta = RandomParameterList ()
        for p_val in theta_values:
            gamma = Gamma (2, 1)
            rand_par = RandomParameter ('p', gamma)
            rand_par.value = p_val
            theta.append (rand_par)
        
        class MockMH (MetropolisHastings):

            def _create_jump_dist (self, theta_t):
                n = theta_t.get_size ()
                s2 = np.eye (n) / 100
                variances = s2.diagonal ()
                theta_values = np.array (theta_t.get_values ())
                mu = np.log (theta_values) - variances / 2
                return MultivariateLognormal (mu, s2)

        mocked_mh = MockMH (theta)
        mean_jump = np.zeros (n)
        for i in range (N):
            jump = mocked_mh.propose_jump (theta)
            mean_jump += jump.get_values ()
        mean_jump /= N
        assert all (abs (mean_jump - theta_values) < 1e-1)


    def test_sample_start (self):
        """ Tests if the start sample is centered on the mean of the 
            prior distribution. """
        n = 10
        N = 1000
        theta = RandomParameterList ()
        for i in range (n):
            gamma = Gamma (2, 0.15)
            rand_par = RandomParameter ('p', gamma)
            theta.append (rand_par)

        class MockMH (MetropolisHastings):
            def _calc_log_likelihood (self, t):
                return 1
        
        start_mean = np.zeros (n)
        analytical_mean = np.ones (n) * 0.3
        for i in range (N):
            mocked_mh = MockMH (theta)
            mocked_mh.start_sample_from_prior ()
            sample = mocked_mh.get_last_sampled (1)[0]
            t = sample[0]
            start_mean += t.get_values ()
        start_mean /= N
        assert all (abs (start_mean - analytical_mean) < 1e-1)


    def test_acceptance_ratio (self):
        """ Tests if the acceptance ratio converges to 0.5 when the 
            mh ratio in a jump proposal is always 0.5. """
        n = 10
        N = 1000
        theta = RandomParameterList ()
        for i in range (n):
            gamma = Gamma (2, 2)
            rand_par = RandomParameter ('p', gamma)
            theta.append (rand_par)

        class MockMH (MetropolisHastings):
            def _calc_log_likelihood (self, t):
                return 1

            def _calc_mh_ratio (self, new_t, new_l, old_t, old_l):
                return 0.5

        mocked_mh = MockMH (theta)
        mocked_mh.get_sample (N)
        acceptance_ratio = mocked_mh.get_acceptance_ratio ()
        assert (abs (acceptance_ratio - .5) < 1e-4)
