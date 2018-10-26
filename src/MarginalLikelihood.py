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


    def estimate_marginal_likelihood (self, experiments, model, 
            theta_prior):
        """ This function estimates the marginal likelihood of a  model.
        """
        n_init = self.__init_iterations
        n_adap = self.__adaptive_iterations
        n_fixed = self.__fixed_iterations
        n_strata = self.__n_strata
        strata_size = self.__strata_size

        mcmc_init = MCMCInitialization (experiments, model, theta_prior)
        start_sample = mcmc_init.get_sample (n_init)

        amcmc = AdaptiveMCMC (model, experiments, start_sample, 
                n_strata, strata_size)
        betas, thetas = amcmc.get_sample (n_adap, n_fixed)

        ml = self.__calculate_marginal_likelihood (model, experiments,
                betas, thetas)
        return ml


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
