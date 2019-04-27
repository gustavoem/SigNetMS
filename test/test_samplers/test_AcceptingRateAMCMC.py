import sys
sys.path.insert (0, '..')

import unittest
import numpy as np
from marginal_likelihood.SBML import SBML
from marginal_likelihood.SBMLtoODES import sbml_to_odes
from marginal_likelihood.PriorsReader import define_sbml_params_priors
from marginal_likelihood.samplers.AcceptingRateAMCMC import AcceptingRateAMCMC
from marginal_likelihood.RandomParameter import RandomParameter
from marginal_likelihood.RandomParameterList import RandomParameterList
from experiment.ExperimentSet import ExperimentSet
from distributions.Gamma import Gamma


class AlwaysAcceptMock (AcceptingRateAMCMC):

    def _calc_mh_ratio (self, new_t, new_l, old_t, old_l):
        return 1

    def get_jump_S (self):
        return list (self._jump_S)

class AlwaysRejectMock (AcceptingRateAMCMC):

    def _calc_mh_ratio (self, new_t, new_l, old_t, old_l):
        return 0

    def get_jump_S (self):
        return list (self._jump_S)


class TestAcceptingRateAMCMC (unittest.TestCase):
    
    def setUp (self):
        sbml = SBML ()
        sbml.load_file ('input/simple_enzymatic.xml')
        self.__model = sbml_to_odes (sbml)
        self.__experiments = ExperimentSet ('input/simple_enzymatic.data')
        self.__theta_priors = define_sbml_params_priors (sbml, 
                'input/simple_enzymatic.priors')


    def test_decreasing_jump (self):
        """ Tests if jump_S actually increases when acceptance ratio
            is too big. """
        model = self.__model
        experiments = self.__experiments
        theta = self.__theta_priors

        mock_mh = AlwaysAcceptMock (theta, model, experiments, 10)
        mock_mh.start_sample_from_prior ()
        
        before = mock_mh.get_jump_S ()
        t, l = mock_mh.get_sample (1)
        for i in range (8):
            mock_mh.manual_jump (t[-1], l[-1])
        after = mock_mh.get_jump_S ()
        assert all (np.array (before) == np.array (after))
        mock_mh.get_sample (1)
        after = mock_mh.get_jump_S ()
        assert all (np.array (before) < np.array (after))


    def test_increasing_jump (self):
        """ Test if jump_S decreases when acceptance ratio is to 
            small. """
        model = self.__model
        experiments = self.__experiments
        theta = self.__theta_priors
            
        mock_mh = AlwaysRejectMock (theta, model, experiments, 10)
        mock_mh.start_sample_from_prior ()

        before = mock_mh.get_jump_S ()
        mock_mh.get_sample (10)
        after = mock_mh.get_jump_S ()
        assert all (np.array (before) > np.array (after))


    def test_get_sample (self):
        """ Tests if we can get a sample. """
        model = self.__model
        experiments = self.__experiments
        theta = self.__theta_priors

        acc_rate_amh = AcceptingRateAMCMC (theta, model, experiments, 
                20, verbose=False)
        acc_rate_amh.start_sample_from_prior ()
        sample = acc_rate_amh.get_sample (60)[0]
        self.assertEqual (len (sample), 60)


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
        
        class JumpMock (AcceptingRateAMCMC):
            
            def __init__ (self, theta):
                n = theta.get_size ()
                self._jump_S = [0.001 for x in range (n)]

        mocked_mh = JumpMock (theta)
        mean_jump = np.zeros (n)
        for i in range (N):
            jump = mocked_mh.propose_jump (theta)
            mean_jump += jump.get_values ()
        mean_jump /= N
        assert all (abs (mean_jump - theta_values) < 1e-1)

