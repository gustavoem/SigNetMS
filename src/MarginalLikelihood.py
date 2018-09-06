# In this module we calculate the marginal likelihood of models given
# observed data using the methodology presented in "Inferring Signaling 
# Pathway Topologies from Multiple Perturbation Measurements of Specific 
# Biochemical Species", Tian-Rui Xu, et. al.

from RandomParameter import RandomParameter
from LikelihoodFunction import LikelihoodFunction
from ODES import ODES
import numpy as np
import random

# NUMBER_OF_STRATA = 100
NUMBER_OF_STRATA = 5
SCHEDULE_POWER = 5
# POPULATION_ITERATIONS = 100000
POPULATION_ITERATIONS = 1000
PROPOSAL_DISTR_STD = 0.3

def estimate_marginal_likelihood (experiments, model):
    """ This function estimates the marginal likelihood of an ODE model,
    """
    M_n = 20
    betas = sample_betas (M_n)
    thetas = get_theta_chains (model, betas)
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


def get_theta_chains (model, betas):
    """ Given a model, construct a list containing all parameters of 
        the model as random variables. """
    theta = []
    params = model.get_all_parameters ()
    for param in params:
        rand_param = RandomParameter (param, 2, 3333.0)
        theta.append (rand_param)
    
    thetas = []
    for b in betas:
        thetas.append (theta.copy ())
    return thetas


# TODO: should cache somehow p (y | theta)
def iterate_thetas (model, thetas, betas, experiments):
    for i in range (POPULATION_ITERATIONS):
        # Local move
        j = random.choice (range (len (thetas)))
        new_theta = propose_theta_jump (thetas[j])
        old_l = set_of_experiments_likelihood (experiments, thetas[j], \
                model)
        new_l = set_of_experiments_likelihood (experiments, new_theta, \
                model)
        r = (new_l / old_l) ** betas[j]
        print ("Local move prob" + str (r))
        if (np.random.uniform () <= r):
            thetas[j] = new_theta

        # Global move
        j = random.choice (range (len (thetas) - 1))
        theta1 = thetas[j]
        theta2 = thetas[j + 1]
        theta1_l = set_of_experiments_likelihood (experiments, theta1, \
                model)
        theta2_l = set_of_experiments_likelihood (experiments, theta2, \
                model)
        t1ot2 = theta1_l / theta2_l
        t2ot1 = theta2_l / theta1_l
        r = (t2ot1) ** betas[j] * (t1ot2) ** betas[j + 1]
        print ("Global move prob" + str (r))
        if (np.random.uniform () <= r):
            aux = thetas[j]
            thetas[j] = thetas[j + 1]
            thetas[j + 1] = aux


def set_of_experiments_likelihood (experiments, theta, model):
    likelihood_f = LikelihoodFunction (model, 0.3)
    params_hash = random_params_to_hash (theta)
    l = 1
    for experiment in experiments:
        X = experiment.values
        var = experiment.var
        t = experiment.times
        l *= likelihood_f.get_experiment_likelihood (X, var, t, \
                params_hash)
    return l


def propose_theta_jump (theta):
    for i in range (len (theta)):
        rand_param = theta[i]
        move = np.random.normal (0, PROPOSAL_DISTR_STD)
        if rand_param.value + move > 0:
            rand_param.value += move
        theta[i] = rand_param
    return theta


def random_params_to_hash (rand_params):
    param_hash = {}
    for param in rand_params:
        param_hash[param.name] = param.value
    return param_hash
