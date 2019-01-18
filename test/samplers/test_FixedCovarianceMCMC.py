import sys
sys.path.insert (0, '../src/')
sys.path.insert (0, '../src/samplers/')
sys.path.insert (0, '../src/distributions/')

import unittest
import numpy as np
from SBML import SBML
from FixedCovarianceMCMC import FixedCovarianceMCMC
from SBMLtoODES import sbml_to_odes
from ExperimentSet import ExperimentSet
from PriorsReader import define_sbml_params_priors
from MultivariateLognormal import MultivariateLognormal
from Gamma import Gamma
from RandomParameterList import RandomParameterList
from RandomParameter import RandomParameter
from CovarianceMatrix import calc_covariance

class TestFixedCovarianceMCMC (unittest.TestCase):
    
    def setUp (self):
        sbml = SBML ()
        sbml.load_file ('input/simple_enzymatic.xml')
        self.__model = sbml_to_odes (sbml)
        self.__experiments = ExperimentSet ('input/simple_enzymatic.data')
        self.__theta_priors = define_sbml_params_priors (sbml, 
                'input/simple_enzymatic.priors')

    def test_get_sample (self):
        model = self.__model
        experiments = self.__experiments
        theta = self.__theta_priors
        covar = self.create_covar_matrix ()

        mh = FixedCovarianceMCMC (theta, model, experiments, covar)
        mh.start_sample_from_prior ()
        sample = mh.get_sample (20)[0]
        self.assertEqual (len (sample), 20)


    def create_covar_matrix (self):
        my_artificial_sample = []
        sample_mean = [0.05, 0.1, 0.2, .3]
        S = np.eye (4) / 5
        mu = np.log (np.array (sample_mean)) - S.diagonal () / 2
        my_sample_dist = MultivariateLognormal (mu, S)
        for i in range (50):
            values = my_sample_dist.rvs ()
            my_artificial_sample.append (values)
        covar = calc_covariance (my_artificial_sample)
        return covar
