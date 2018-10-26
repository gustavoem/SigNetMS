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
        start_sample = mcmc_init.get_sample (self.__init_iterations)
        N1, N2 = self.__adaptive_iterations, self.__fixed_iterations
        amcmc = AdaptiveMCMC (model, experiments, start_sample, 
                self.__n_strata, self.__strata_size)
        betas, thetas = amcmc.get_sample (N1, N2)
        ml = self.__calculate_marginal_likelihood (model, experiments,
                betas, thetas)
        return ml


    def __get_theta (self, sbml, model):
        """ Given a model, construct a list containing all parameters of 
        the model as random variables. """
        theta = RandomParameterList ()
        params = model.get_all_parameters ()
        for param in params:
            param_original_name = sbml.get_original_param_name (param)
            # if re.search ("Km", param_original_name):
                # rand_param = RandomParameter (param, 2.0, 3333.0)
            # else:
                # rand_param = RandomParameter (param, 1.1, 9.0)

            if param_original_name == "k1":
                rand_param = RandomParameter (param, 2.0, 0.01)
            if param_original_name == "d1" or param_original_name == "kcat":
                rand_param = RandomParameter (param, 2.0, 0.1)

            # rand_param = RandomParameter (param, 4, .5)
            theta.append (rand_param)
        # sigma = RandomParameter ("experimental_sigma", 2.0, 2.6)
        sigma = RandomParameter ("experimental_sigma", 1, 1)
        theta.set_experimental_error_parameter (sigma)
        return theta


    def __calculate_marginal_likelihood (self, model, experiments, 
            betas, thetas):
        """ Given a list with samples, calculates the marginal 
            likelihood. """
        ml = 0
        j = 0 
        n_strata = self.__n_strata
        strata_size = self.__strata_size
        sched_power = AdaptiveMCMC.get_sched_power ()
        likeli_f = LikelihoodFunction (model)
    
        print ("Estimating marginal likelihood")
        
        for i in range (n_strata):
            strata_start = (n_strata) ** sched_power
            strata_end = ((i + 1) / n_strata) ** sched_power
            del_strata = strata_end - strata_start
                
            strat_sum = 0
            while j < len (betas) and betas[j] < strata_end:
                p_y_given_theta = likeli_f.get_experiments_likelihood (
                    experiments, thetas[j])
                # print ("\tTheta: ", end='')
                # for r in thetas[j]:
                    # print (r.value, end=' ')
                # print ("\n\tLikelihood: " + str (p_y_given_theta) + "\n")
                if p_y_given_theta > 0:
                    strat_sum += np.log (p_y_given_theta)
                j += 1

            strat_sum *= (del_strata / strata_size)
            ml += strat_sum
        return ml
