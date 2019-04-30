# In this module we calculate the marginal likelihood of models given
# observed data using the methodology presented in "Inferring Signaling 
# Pathway Topologies from Multiple Perturbation Measurements of Specific 
# Biochemical Species", Tian-Rui Xu, et. al.

from marginal_likelihood.RandomParameter import RandomParameter
from marginal_likelihood.RandomParameterList import RandomParameterList
from marginal_likelihood.LikelihoodFunction import LikelihoodFunction
from marginal_likelihood.ODES import ODES
from marginal_likelihood.samplers.AcceptingRateAMCMC import \
        AcceptingRateAMCMC
from marginal_likelihood.samplers.AdaptingCovarianceMCMC import \
        AdaptingCovarianceMCMC
from marginal_likelihood.samplers.FixedCovarianceMCMC import \
        FixedCovarianceMCMC
from marginal_likelihood.samplers.PopulationalMCMC import \
        PopulationalMCMC
from utils import plot_theta_var_sample
import numpy as np
import random
import re


class MarginalLikelihood:
    """ This class is able to perform an adaptive MCMC sampling to 
        estimate the likelihood of a model given experimental data. """
    
    def __init__ (self, phase1_iterations, sigma_update_n, 
            phase2_iterations, phase3_iterations, n_strata, 
            strata_size, verbose=False):
        """ Default constructor. phase1_iterations is the number of 
            iterations performed by the AcceptingRateAMCMC, which is
            an adaptive sampler that performs independent MCMC on each
            system variable and adapts the jump covariance with the goal
            of having acceptance rate between 0.25 and 0.75. 
            sigma_update_n is the number of iterations before updating 
            sigma in intial phase. phase2_iterations is the number of 
            iterations performed in the second phase of sampling with 
            an adaptive sampling algorithm that has jumps distributed 
            with a probability function that has as covariance matrix
            a matrix of covariance estimated from the current sample.
            phase3_iterations is the number of iterations performed by
            the PopulationalMCMC algorithm. n_strata is the number of 
            strata used in the populational algorithm, and strata_size 
            is the number of individuals per strata."""
        self.__phase1_iterations = phase1_iterations
        self.__phase2_iterations = phase2_iterations
        self.__phase3_iterations = phase3_iterations
        self.__sigma_update_n = sigma_update_n
        self.__n_strata = n_strata
        self.__strata_size = strata_size
        self.__verbose = verbose


    def estimate_marginal_likelihood (self, experiments, model, 
            theta_prior):
        """ This function estimates the marginal likelihood of a  model.
        """
        n_acc = self.__phase1_iterations
        n_adap_cov = self.__phase2_iterations
        n_pop = self.__phase3_iterations
        n_strata = self.__n_strata
        strata_size = self.__strata_size
        verbose = self.__verbose
    
        # Phase 1
        acc_mcmc = AcceptingRateAMCMC (theta_prior, model, experiments,
                self.__sigma_update_n, verbose=verbose)
        acc_mcmc.start_sample_from_prior ()
        sample, likelis = acc_mcmc.get_sample (n_acc)
        
        if self.__verbose:
            print ("Phase 2 starts.")
        # Phase 2
        adap_cov_mcmc = AdaptingCovarianceMCMC (theta_prior, model, 
                experiments, verbose=verbose)
        adap_cov_mcmc.define_start_sample (sample, likelis)
        sample, likelis = adap_cov_mcmc.get_sample (n_adap_cov)
        
        if self.__verbose:
            print ("Phase 3 starts.")
        # Phase 3
        fc_mcmcs = self.__create_fcmcmc_samplers (adap_cov_mcmc, 
                n_strata * strata_size, experiments, model, theta_prior)
        pop_mcmc = PopulationalMCMC (n_strata, strata_size, fc_mcmcs,
                verbose=verbose)
        pop_mcmc.get_sample (n_pop)
        betas, thetas, log_ls = pop_mcmc.get_last_sampled (n_pop // 4)
        
        if self.__verbose:
            print ("Sampling ended. Here are the sampled parameters" + \
                    " separated by temperature")
            for i in range (len (betas)):
                temp = betas[i]
                print ("Sample for t = " + str (temp))
                for j in range (len (thetas[i])):
                    theta = thetas[i][j]
                    likeli = log_ls[i][j]
                    print ("parameter = " + str (theta.get_values ()) +\
                            " likelihood = " + str (likeli))

        ml = self.__calculate_marginal_likelihood (betas, thetas, 
                log_ls)
        return ml


    def __create_fcmcmc_samplers (self, adap_cov_mcmc, m, experiments,
            model, theta_prior):
        """ Creates a list of FixedCovarianceMCMC objects that has the
            same jump covariance matrix as AdaptiveCovarianceMatrix. """
        S = adap_cov_mcmc.get_jump_covariance ()
        start_t, start_l = adap_cov_mcmc.get_last_sampled (1)
        fc_mcmcs = []
        for i in range (m):
            sampler = FixedCovarianceMCMC (theta_prior, model, 
                    experiments, S)
            sampler.define_start_sample (start_t, start_l)
            fc_mcmcs.append (sampler)
        return fc_mcmcs


    def __calculate_marginal_likelihood (self, betas, thetas, \
            likelihoods):
        """ Given a list with samples, calculates the marginal 
            likelihood. """
        ml = 0
        
        if self.__verbose:
            print ("Estimating marginal likelihood")

        exp_gv_tip1 = 0
        power_sample_ip1 = thetas[0]
        for j in range (len (power_sample_ip1)):
            log_l = likelihoods[0][j]
            exp_gv_tip1 += log_l
        exp_gv_tip1 /= len (power_sample_ip1)

        for i in range (len (betas) - 1):
            tip1 = betas[i + 1]
            ti = betas[i]
                        
            exp_gv_ti = exp_gv_tip1
            exp_gv_tip1 = 0
            power_sample_ip1 = thetas[i + 1]
            for j in range (len (power_sample_ip1)):
                log_l = likelihoods[i + 1][j]
                exp_gv_tip1 += log_l
            exp_gv_tip1 /= len (power_sample_ip1)

            ml += (tip1 - ti) * (exp_gv_tip1 + exp_gv_ti)
            
            if self.__verbose:
                print ("first temperature = " + str (ti))
                print ("second temperature = " + str (tip1))
                print ("log ml given first temperature = " + \
                        str (exp_gv_ti))
                print ("log ml given second temperature = " + \
                        str (exp_gv_tip1))

        ml /= 2
        print ("Calculated log marginal likelihood: " + str (ml))
        return ml
