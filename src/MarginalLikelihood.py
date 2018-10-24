# In this module we calculate the marginal likelihood of models given
# observed data using the methodology presented in "Inferring Signaling 
# Pathway Topologies from Multiple Perturbation Measurements of Specific 
# Biochemical Species", Tian-Rui Xu, et. al.

from RandomParameter import RandomParameter
from RandomParameterList import RandomParameterList
from LikelihoodFunction import LikelihoodFunction
from MCMCInitialization import MCMCInitialization
from AdaptiveMCMC import AdaptiveMCMC
from ODES import ODES
import numpy as np
import random
import re

# NUMBER_OF_STRATA = 100
# NUMBER_OF_STRATA = 5
# SCHEDULE_POWER = 5
# POPULATION_ITERATIONS = 100000
# POPULATION_ITERATIONS = 1000


class MarginalLikelihood:
    """ This class is able to perform an adaptive MCMC sampling to 
        estimate the likelihood of a model given experimental data. """
    
    def __init__ (self, init_iterations, adaptive_iterations, 
            fixed_iterations, n_strata, strata_size):
        """ Default constructor. init_iterations is the number of 
            iterations performed by the MCMCInitialize sampler, which is
            an adaptive sampler that performs independent MCMC on each
            system variable. adaptive_iterations is the number of 
            iterations performed in the adaptive phase of AdaptiveMCMC
            object and fixed_iterations is the number of iterations in 
            the fixed phase of the same object. n_strata is the number 
            of strata used in the populational phase of AdaptiveMCMC 
            (fixed phase), and strata_size is the number of individuals
            per strata."""
        self.__init_iterations = init_iterations
        self.__adaptive_iterations = adaptive_iterations
        self.__fixed_iterations = fixed_iterations
        self.__n_strata = n_strata
        self.__strata_size = strata_size


    def estimate_marginal_likelihood (self, experiments, sbml, model):
        """ This function estimates the marginal likelihood of a  model.
        """
        theta = self.__get_theta (sbml, model)
        mcmc_init = MCMCInitialization (theta, model, experiments)
        start_sample = mcmc_init.get_sample (100)
        amcmc = AdaptiveMCMC (model, experiments, start_sample)
        betas, thetas = amcmc.get_sample (100, 200)
        ml = calculate_marignal_likelihood (betas, thetas)


    def __get_theta (self, sbml, model):
        """ Given a model, construct a list containing all parameters of 
        the model as random variables. """
        theta = RandomParameterList ()
        params = model.get_all_parameters ()
        for param in params:
            param_original_name = sbml.get_original_param_name (param)
            if re.search ("Km", param_original_name):
                rand_param = RandomParameter (param, 2.0, 3333.0)
            else:
                rand_param = RandomParameter (param, 1.1, 9.0)

            if param_original_name == "k1":
                rand_param = RandomParameter (param, 2.0, 0.01)
            if param_original_name == "d1" or param_original_name == "kcat":
                rand_param = RandomParameter (param, 2.0, 0.1)
            theta.append (rand_param)
        sigma = RandomParameter ("experimental_sigma", 2.0, 2.6)
        theta.set_experimental_error_parameter (sigma)
        return theta


    def __calculate_marginal_likelihood (self, betas, thetas):
        """ Given a list with samples, calculates the marginal 
            likelihood. """
        ml = 0


