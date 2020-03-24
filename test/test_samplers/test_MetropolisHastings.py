import sys
sys.path.insert (0, '..')

import unittest
import numpy as np
from marginal_likelihood.samplers.MetropolisHastings import \
        MetropolisHastings
from model.RandomParameterList import RandomParameterList
from model.RandomParameter import RandomParameter
from distributions.Gamma import Gamma
from distributions.MultivariateLognormal import MultivariateLognormal


class MHJumpMock (MetropolisHastings):
    def _create_jump_dist (self, theta_t):
        n = theta_t.get_size ()
        s2 = np.eye (n) / 100
        variances = s2.diagonal ()
        theta_values = np.array (theta_t.get_values ())
        mu = np.log (theta_values) - variances / 2
        return MultivariateLognormal (mu, s2)


class MHLikelihoodMock (MetropolisHastings):
    def _calc_log_likelihood (self, t):
        return 1


class MHFullMock (MetropolisHastings):
    def _calc_log_likelihood (self, t):
        return 1

    def _calc_mh_ratio (self, new_t, new_l, old_t, old_l):
        return 0.5
            
    def _create_jump_dist (self, theta_t):
        n = theta_t.get_size ()
        s2 = np.eye (n) / 100
        variances = s2.diagonal ()
        theta_values = np.array (theta_t.get_values ())
        mu = np.log (theta_values) - variances / 2
        return MultivariateLognormal (mu, s2)


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
        
        mocked_mh = MHJumpMock (theta)
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

                
        start_mean = np.zeros (n)
        analytical_mean = np.ones (n) * 0.3
        for i in range (N):
            mocked_mh = MHLikelihoodMock (theta)
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
        N = 3000
        theta = RandomParameterList ()
        for i in range (n):
            gamma = Gamma (2, 2)
            rand_par = RandomParameter ('p', gamma)
            theta.append (rand_par)

        mocked_mh = MHFullMock (theta)
        mocked_mh.start_sample_from_prior ()
        mocked_mh.get_sample (N)
        acceptance_ratio = mocked_mh.get_acceptance_ratio ()
        assert (abs (acceptance_ratio - .5) < 1e-1)

    def test_manual_jump (self):
        """ Tests if one can perform a manual jump. """
        n = 10
        theta = RandomParameterList ()
        for i in range (n):
            gamma = Gamma (2, 2)
            rand_par = RandomParameter ('p', gamma)
            theta.append (rand_par)
        
        mocked_mh = MHLikelihoodMock (theta)
        mocked_mh.start_sample_from_prior ()
        mocked_mh.manual_jump (theta, 1)

        sample = mocked_mh.get_last_sampled (2)[0]
        last_theta = sample[-1]
        theta_vals = np.array (theta.get_values ())
        last_theta_vals = np.array (last_theta.get_values ())
        assert all (theta_vals == last_theta_vals)

    
    def test_get_last_sampled (self):
        """ Tests if one can get the last N sampled parameters. """
        n = 10
        N = 20
        theta = RandomParameterList ()
        for i in range (n):
            gamma = Gamma (2, 2)
            rand_par = RandomParameter ('p', gamma)
            theta.append (rand_par)

        mocked_mh = MHFullMock (theta)
        mocked_mh.start_sample_from_prior ()
        mocked_mh.get_sample (2 * N)
        sample = mocked_mh.get_last_sampled (N)[0]
        assert (len (sample) > N * .25)


    def test_samples_target (self):
        """ With a controlled target distribution, tests if the MH
            successfully generates a sample from the target. """
        n = 10
        N = 5000
        theta = RandomParameterList ()
        for i in range (n):
            gamma = Gamma (1, 1)
            p = RandomParameter ('p', gamma)
            theta.append (p)

        class MHMultivariateLognormalTargetMock (MetropolisHastings):
     
            def __init__ (self, theta):
                super ().__init__ (theta, verbose=False)
                n = theta.get_size ()
                mu = np.ones (n)
                Sigma = np.eye (n) / 10
                self.target_distribution = MultivariateLognormal (mu, 
                        Sigma)

            def _calc_log_likelihood (self, t):
                t_values = t.get_values ()
                l = self.target_distribution.pdf (t_values)
                return np.log (l)

            def _create_jump_dist (self, theta_t):
                n = theta_t.get_size ()
                S = np.eye (n) / 100
                variances = S.diagonal ()
                theta_values = np.array (theta_t.get_values ())
                mu = np.log (theta_values) - variances / 2
                return MultivariateLognormal (mu, S)

            def _calc_mh_ratio (self, new_t, new_l, old_t, old_l):
                J_gv_new = self._create_jump_dist (new_t)
                J_gv_old = self._create_jump_dist (old_t)
                p_old_gv_new = J_gv_new.pdf (old_t.get_values ())
                p_new_gv_old = J_gv_old.pdf (new_t.get_values ())
                l_ratio = np.exp (new_l - old_l)
                r = (l_ratio) * (p_old_gv_new / p_new_gv_old)
                return r


        mocked_mh = MHMultivariateLognormalTargetMock (theta)
        mocked_mh.start_sample_from_prior ()
        mocked_mh.get_sample (N)[0]
        
        # let's throw away the first 10% of samples
        sample = mocked_mh.get_last_sampled (int (N * .9))[0]

        mu = np.ones (n)
        S  = np.eye (n) / 10
        analytic_mean = np.exp (mu + S.diagonal () / 2)
        sample_mean = np.zeros (n) 
        for t in sample:
            sample_mean += t.get_values ()
        sample_mean /= len (sample)
        
        diff = analytic_mean - sample_mean
        diff_norm2 = np.sqrt (diff.dot (diff))
        assert (diff_norm2 < 1)
