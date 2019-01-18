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
from MultivariateLognormal import MultivariateLognormal
from Gamma import Gamma
from RandomParameterList import RandomParameterList
from RandomParameter import RandomParameter

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
        self.assertRaises (ValueError, mh.get_sample, 1)
        mh.start_sample_from_prior ()
        self.assertRaises (ValueError, mh.get_sample, 1)


    def test_get_sample (self):
        model = self.__model
        experiments = self.__experiments
        theta = self.__theta_priors
        start_sample, start_likels = self.create_starting_sample ()

        mh = AdaptingCovarianceMCMC (theta, model, experiments)
        mh.start_sample_from_prior ()
        mh.define_start_sample (start_sample, start_likels)
        sample = mh.get_sample (20)[0]
        self.assertEqual (len (sample), 20)


    def create_starting_sample (self):
        my_artificial_sample = []
        log_likelihoods = []
        sample_mean = [0.05, 0.1, 0.2, .3]
        S = np.eye (4) / 5
        mu = np.log (np.array (sample_mean)) - S.diagonal () / 2
        my_sample_dist = MultivariateLognormal (mu, S)
        for i in range (50):
            theta = RandomParameterList ()
            values = my_sample_dist.rvs ()
            for v in values:
                p = RandomParameter ('p', Gamma (2, 2))
                p.value = v
                theta.append (p)
            log_likelihoods.append (1)
            my_artificial_sample.append (theta)
        return (my_artificial_sample, log_likelihoods)
