import sys
sys.path.insert (0, '../src/')
sys.path.insert (0, '../src/samplers/')
sys.path.insert (0, '../src/distributions/')

import unittest
import numpy as np
from SBML import SBML
from AdaptingCovarianceMCMC import AdaptingCovarianceMCMC
from SBMLtoODES import sbml_to_odes
from ExperimentSet import ExperimentSet
from PriorsReader import define_sbml_params_priors
from Gamma import Gamma

class TestAdaptingCovarianceMCMC (unittest.TestCase):
    
    def setUp (self):
        sbml = SBML ()
        sbml.load_file ('input/simple_enzymatic.xml')
        self.__model = sbml_to_odes (sbml)
        self.__experiments = ExperimentSet ('input/simple_enzymatic.data')
        self.__theta_priors = define_sbml_params_priors (sbml, 
                'input/simple_enzymatic.priors')

    def test_except_ifnot_start_sample (self):
        model = self.__model
        experiments = self.__experiments
        theta = self.__theta_priors

        mh = AdaptingCovarianceMCMC (theta, model, experiments)
        self.assertRaises (ValueError, mh.get_sample (1))
