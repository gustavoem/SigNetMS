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


def estimate_marginal_likelihood (experiments, sbml, model):
    """ This function estimates the marginal likelihood of a  model """
    M_n = 20
    theta = get_theta (sbml, model)
    mcmc_init = MCMCInitialization (theta, model, experiments)
    start_sample = mcmc_init.get_sample (100)
    amcmc = AdaptiveMCMC (model, experiments, start_sample)
    amcmc.get_sample (100, 200)


def get_theta (sbml, model):
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
