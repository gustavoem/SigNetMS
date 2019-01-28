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

        betas = PopulationalMCMC.sample_betas (n_strata, strata_size)
        fc_mcmcs = []

        print ("Phase 1 and 2 starts.")
        for b in betas:
            if verbose:
                print ("Temperature t = " + str (b))

            # Phase 1
            acc_mcmc = AcceptingRateAMCMC (theta_prior, model, 
                    experiments, self.__sigma_update_n, verbose=verbose)
            acc_mcmc.set_temperature (b)
            acc_mcmc.start_sample_from_prior ()
            sample, likelis = acc_mcmc.get_sample (n_acc)

            # Phase 2
            adap_cov_mcmc = AdaptingCovarianceMCMC (theta_prior, model, 
                    experiments, verbose=verbose)
            adap_cov_mcmc.set_temperature (b)
            adap_cov_mcmc.define_start_sample (sample, likelis)
            sample, likelis = adap_cov_mcmc.get_sample (n_adap_cov)

            # Construct phase 3 local temperature sampler
            S = adap_cov_mcmc.get_jump_covariance ()
            fc_mcmc = FixedCovarianceMCMC (theta_prior, model, 
                    experiments, S, t=b, verbose=verbose)
            theta = sample[-1]
            log_likeli = likelis[-1]
            fc_mcmc.define_start_sample ([theta], [log_likeli])
            fc_mcmcs.append (fc_mcmc)
           
        if verbose:
            print ("Phase 3 starts.")
        # Phase 3
        # fc_mcmcs = self.__create_fcmcmc_samplers (adap_cov_mcmc, 
                # n_strata * strata_size, experiments, model, theta_prior)
        pop_mcmc = PopulationalMCMC (n_strata, strata_size, fc_mcmcs,
                betas=betas, verbose=verbose)
        betas, thetas, log_ls = pop_mcmc.get_sample (n_pop)
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
        j = 0 
        n_strata = self.__n_strata
        strata_size = self.__strata_size
        sched_power = PopulationalMCMC.get_sched_power ()
        
        if self.__verbose:
            print ("Estimating marginal likelihood")
        for i in range (n_strata):
            strata_start = (i / n_strata) ** sched_power
            strata_end = ((i + 1) / n_strata) ** sched_power
            del_strata = strata_end - strata_start
            
            if self.__verbose:
                print ("strata_start = " + str (strata_start))
                print ("strata_end = " + str (strata_end))
                print ("del_strata = " + str (del_strata))
                
            strat_sum = 0
            while j < len (betas) and betas[j] <= strata_end:
                log_p_y_given_theta = likelihoods[j]
                strat_sum += log_p_y_given_theta
                if self.__verbose:
                    print ("\tTheta: ", end='')
                    for r in thetas[j]:
                        print (r.value, end=' ')
                    print ("\n\tLikelihood: " + str (log_p_y_given_theta) \
                            + "\n")
                    print ("Strat_sum = " + str (strat_sum))
                j += 1

            strat_sum *= (del_strata / strata_size)
            ml += strat_sum
        print ("Calculated marginal likelihood: " + str (ml))
        return ml
