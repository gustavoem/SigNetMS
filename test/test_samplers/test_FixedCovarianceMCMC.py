import sys
sys.path.insert (0, '..')

import unittest
import numpy as np
from model.SBML import SBML
from model.SBMLtoODES import sbml_to_odes
from model.PriorsReader import define_sbml_params_priors
from experiment.ExperimentSet import ExperimentSet
from distributions.Gamma import Gamma
from distributions.MultivariateLognormal import MultivariateLognormal
from model.RandomParameterList import RandomParameterList
from model.RandomParameter import RandomParameter
from covariance_estimate import calc_covariance_diagonal
from marginal_likelihood.samplers.FixedCovarianceMCMC import \
        FixedCovarianceMCMC

class TestFixedCovarianceMCMC (unittest.TestCase):
    
    def setUp (self):
        sbml = SBML ()
        sbml.load_file ('input/simple_enzymatic.xml')
        self.__model = sbml_to_odes (sbml)
        self.__experiments = \
                ExperimentSet ('input/simple_enzymatic.data')
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
        assert len (sample) > 1

    def test_covariance_is_fixed (self):
        model = self.__model
        experiments = self.__experiments
        theta = self.__theta_priors
        covar = self.create_covar_matrix ()
        mh = FixedCovarianceMCMC (theta, model, experiments, covar)
        mh.start_sample_from_prior ()
        mh.get_sample (2)
        new_covar = mh.get_jump_covariance ()
        assert np.array_equal (new_covar, covar)
        covar = new_covar
        mh.get_sample (1)
        new_covar = mh.get_jump_covariance ()
        assert np.array_equal (new_covar, covar)


    def create_covar_matrix (self):
        my_artificial_sample = []
        sample_mean = [0.05, 0.1, 0.2, .3]
        S = np.eye (4) / 5
        mu = np.log (np.array (sample_mean)) - S.diagonal () / 2
        my_sample_dist = MultivariateLognormal (mu, S)
        for i in range (50):
            values = my_sample_dist.rvs ()
            my_artificial_sample.append (values)
        covar = calc_covariance_diagonal (my_artificial_sample)
        return covar
