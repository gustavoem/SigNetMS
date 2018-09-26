# In this module we calculate the marginal likelihood of models given
# observed data using the methodology presented in "Inferring Signaling 
# Pathway Topologies from Multiple Perturbation Measurements of Specific 
# Biochemical Species", Tian-Rui Xu, et. al.

from RandomParameter import RandomParameter
from LikelihoodFunction import LikelihoodFunction
from ODES import ODES
import numpy as np
import random
import re

# NUMBER_OF_STRATA = 100
NUMBER_OF_STRATA = 5
SCHEDULE_POWER = 5
# POPULATION_ITERATIONS = 100000
POPULATION_ITERATIONS = 1000
PROPOSAL_DISTR_STD = 0.3


def estimate_marginal_likelihood (experiments, sbml, model):
    """ This function estimates the marginal likelihood of a  model """
    M_n = 20
    betas = sample_betas (M_n)
    thetas = get_theta_chains (sbml, model, betas)
    iterate_thetas (model, thetas, betas, experiments)


def sample_betas (M_n):
    """ Samples M_n beta from each of the NUMBER_OF_STRATA strata."""
    betas = []
    for i in range (NUMBER_OF_STRATA):
        strata_start = (i / NUMBER_OF_STRATA) ** SCHEDULE_POWER
        strata_end = ((i + 1) / NUMBER_OF_STRATA) ** SCHEDULE_POWER
        for j in range (M_n):
            x = np.random.uniform (strata_start, strata_end)
            betas.append (x)
        betas.sort ()
    return betas


def get_theta_chains (sbml, model, betas):
    """ Given a model, construct a list containing all parameters of 
        the model as random variables. """
    theta = []
    params = model.get_all_parameters ()
    for param in params:
        param_original_name = sbml.get_original_param_name (param)
        if re.search ("Km", param_original_name):
            rand_param = RandomParameter (param, 2.0, 3333.0)
        else:
            rand_param = RandomParameter (param, 1.1, 9.0)
        theta.append (rand_param)
    
    thetas = []
    for b in betas:
        thetas.append (theta.copy ())
    return thetas


# TODO: should cache somehow p (y | theta)
def iterate_thetas (model, thetas, betas, experiments):
    sigma = np.random.gamma (2.0, 3333.0)
    for i in range (POPULATION_ITERATIONS):
        # Local move
        j = random.choice (range (len (thetas)))
        new_theta = propose_theta_jump (thetas[j])
        print ("Calculating likelihood with old theta.")
        old_l = set_of_experiments_likelihood (experiments, thetas[j], \
                model, sigma)
        print ("Calculating likelihood with new theta.")
        new_l = set_of_experiments_likelihood (experiments, new_theta, \
                model, sigma)
        print ("old theta: ")
        for r in thetas[j]:
            print (r.value, end=' ')
        print ("\nnew theta: ")
        for r in new_theta:
            print (r.value, end=' ')
        print ("\nOld theta likelihood " + str (old_l))
        print ("New theta likelihood " + str (new_l), end='\n\n')
        r = (new_l / old_l) ** betas[j]
        # print ("Local move prob" + str (r) + "\n\n")
        if (np.random.uniform () <= r):
            thetas[j] = new_theta

        # Global move
        j = random.choice (range (len (thetas) - 1))
        theta1 = thetas[j]
        theta2 = thetas[j + 1]
        theta1_l = set_of_experiments_likelihood (experiments, theta1, \
                model, sigma)
        theta2_l = set_of_experiments_likelihood (experiments, theta2, \
                model, sigma)
        t1ot2 = theta1_l / theta2_l
        t2ot1 = theta2_l / theta1_l
        r = (t2ot1) ** betas[j] * (t1ot2) ** betas[j + 1]
        # print ("Global move prob " + str (r))
        if (np.random.uniform () <= r):
            aux = thetas[j]
            thetas[j] = thetas[j + 1]
            thetas[j + 1] = aux


def set_of_experiments_likelihood (experiments, theta, model, sigma):
    likelihood_f = LikelihoodFunction (model, sigma)
    l = likelihood_f.get_experiments_likelihood (experiments, theta)
    return l


def propose_theta_jump (theta):
    new_theta = list (theta)
    for i in range (len (theta)):
        ti = theta[i]
        rand_param = RandomParameter (ti.name, ti.get_a (), ti.get_b ())
        rand_param.value = ti.value
        std = rand_param.value
        move = np.random.normal (0, std)
        if rand_param.value + move > 0:
            rand_param.value += move
        new_theta[i] = rand_param
    return new_theta