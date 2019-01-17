import sys
sys.path.insert (0, '../src/')
sys.path.insert (0, '../src/samplers/')
sys.path.insert (0, '../src/distributions/')

import unittest
import numpy as np
from SBML import SBML
from AcceptingRateAMCMC import AcceptingRateAMCMC
from SBMLtoODES import sbml_to_odes
from ExperimentSet import ExperimentSet
from PriorsReader import define_sbml_params_priors
from Gamma import Gamma


class AlwaysAcceptMock (AcceptingRateAMCMC):

    def _calc_mh_ratio (self, new_t, new_l, old_t, old_l):
        return 1

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
        """ Tests if jump_S actually decreases when acceptance ratio
            is too small. """
        model = self.__model
        experiments = self.__experiments
        theta = self.__theta_priors

        mock_mh = AlwaysAcceptMock (theta, model, experiments, 100)
        mock_mh.start_sample_from_prior ()
        
        before = mock_mh.get_jump_S ()
        mock_mh.get_sample (99)
        after = mock_mh.get_jump_S ()
        assert all (np.array (before) == np.array (after))
        mock_mh.get_sample (1)
        after = mock_mh.get_jump_S ()
        assert all (np.array (before) != np.array (after))