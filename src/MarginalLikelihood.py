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
    start_sample = mcmc_init.get_sample (1000)
    print (start_sample)
    amcmc = AdaptiveMCMC (model, experiments, start_sample)


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
    return theta


# def sample_betas (M_n):
    # """ Samples M_n beta from each of the NUMBER_OF_STRATA strata."""
    # betas = []
    # for i in range (NUMBER_OF_STRATA):
        # strata_start = (i / NUMBER_OF_STRATA) ** SCHEDULE_POWER
        # strata_end = ((i + 1) / NUMBER_OF_STRATA) ** SCHEDULE_POWER
        # for j in range (M_n):
            # x = np.random.uniform (strata_start, strata_end)
            # betas.append (x)
        # betas.sort ()
    # return betas


# def iterate_thetas (model, thetas, betas, experiments):
    # sigma = np.random.gamma (2.0, 3333.0)
    # for i in range (POPULATION_ITERATIONS):
        # Local move
        # j = random.choice (range (len (thetas)))
        # new_theta = propose_theta_jump (thetas[j])
        # print ("Calculating likelihood with old theta.")
        # old_l = set_of_experiments_likelihood (experiments, thetas[j], \
                # model, sigma)
        # print ("Calculating likelihood with new theta.")
        # new_l = set_of_experiments_likelihood (experiments, new_theta, \
                # model, sigma)
        # print ("old theta: ")
        # for r in thetas[j]:
            # print (r.value, end=' ')
        # print ("\nnew theta: ")
        # for r in new_theta:
            # print (r.value, end=' ')
        # print ("\nOld theta likelihood " + str (old_l))
        # print ("New theta likelihood " + str (new_l), end='\n\n')
        # r = (new_l / old_l) ** betas[j]
        # print ("Local move prob" + str (r) + "\n\n")
        # if (np.random.uniform () <= r):
            # thetas[j] = new_theta

        # Global move
        # j = random.choice (range (len (thetas) - 1))
        # theta1 = thetas[j]
        # theta2 = thetas[j + 1]
        # theta1_l = set_of_experiments_likelihood (experiments, theta1, \
                # model, sigma)
        # theta2_l = set_of_experiments_likelihood (experiments, theta2, \
                # model, sigma)
        # t1ot2 = theta1_l / theta2_l
        # t2ot1 = theta2_l / theta1_l
        # r = (t2ot1) ** betas[j] * (t1ot2) ** betas[j + 1]
        # print ("Global move prob " + str (r))
        # if (np.random.uniform () <= r):
            # aux = thetas[j]
            # thetas[j] = thetas[j + 1]
            # thetas[j + 1] = aux
