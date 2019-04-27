import sys
sys.path.insert (0, '..')

import unittest
import numpy as np
from marginal_likelihood.SBML import SBML
from marginal_likelihood.samplers.FixedCovarianceMCMC import \
        FixedCovarianceMCMC
from marginal_likelihood.samplers.PopulationalMCMC import \
        PopulationalMCMC
from marginal_likelihood.SBMLtoODES import sbml_to_odes
from marginal_likelihood.RandomParameterList import RandomParameterList
from marginal_likelihood.RandomParameter import RandomParameter
from marginal_likelihood.CovarianceMatrix import calc_covariance
from marginal_likelihood.PriorsReader import define_sbml_params_priors
from distributions.MultivariateLognormal import MultivariateLognormal
from distributions.Gamma import Gamma
from experiment.ExperimentSet import ExperimentSet

class TestFixedCovarianceMCMC (unittest.TestCase):
    
    def setUp (self):
        sbml = SBML ()
        sbml.load_file ('input/simple_enzymatic.xml')
        self.__model = sbml_to_odes (sbml)
        self.__experiments = ExperimentSet ('input/simple_enzymatic.data')
        self.__theta_priors = define_sbml_params_priors (sbml, 
                'input/simple_enzymatic.priors')


    def test_get_sample (self):
        n_strata = 10
        strata_size = 2
        fcmcmcs = self.create_list_of_fcmcmc (n_strata * \
                strata_size + 2)
        pop_mcmc = PopulationalMCMC (n_strata, strata_size, fcmcmcs,
                verbose=False)
        sample = pop_mcmc.get_sample (10)[0]
        self.assertEqual (len (sample), 22)


    def create_list_of_fcmcmc (self, n):
        """ Creates a list of FixedCovarianceMCMC. """
        model = self.__model
        experiments = self.__experiments
        theta = self.__theta_priors
        fcmcmcs = []
        for i in range (n):
            S = self.create_covar_matrix ()
            sampler = FixedCovarianceMCMC (theta, model, experiments, S,
                    verbose=False)
            sampler.start_sample_from_prior ()
            fcmcmcs.append (sampler)
        return fcmcmcs


    def create_covar_matrix (self):
        my_artificial_sample = []
        sample_mean = [0.05, 0.1, 0.2, .3]
        S = np.eye (4) / 5
        my_sample_dist = MultivariateLognormal.\
                create_lognormal_with_shape (sample_mean, S)
        for i in range (50):
            values = my_sample_dist.rvs ()
            my_artificial_sample.append (values)
        covar = calc_covariance (my_artificial_sample)
        return covar
