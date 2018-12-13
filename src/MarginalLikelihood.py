# In this module we calculate the marginal likelihood of models given
# observed data using the methodology presented in "Inferring Signaling 
# Pathway Topologies from Multiple Perturbation Measurements of Specific 
# Biochemical Species", Tian-Rui Xu, et. al.

from RandomParameter import RandomParameter
from RandomParameterList import RandomParameterList
from LikelihoodFunction import LikelihoodFunction
from MCMCInitialization import MCMCInitialization
from DiscreteLaplacian import DiscreteLaplacian
from AdaptiveMCMC import AdaptiveMCMC
from ODES import ODES
import numpy as np
import random
import re


class MarginalLikelihood:
    """ This class is able to perform an adaptive MCMC sampling to 
        estimate the likelihood of a model given experimental data. """
    
    def __init__ (self, init_iterations, sigma_update_n, 
            adaptive_iterations, fixed_iterations, n_strata, 
            strata_size):
        """ Default constructor. init_iterations is the number of 
            iterations performed by the MCMCInitialize sampler, which is
            an adaptive sampler that performs independent MCMC on each
            system variable. sigma_update_n is the number of iterations
            before updating sigma in intial phase. adaptive_iterations 
            is the number of iterations performed in the adaptive phase 
            of AdaptiveMCMC object and fixed_iterations is the number of 
            iterations in the fixed phase of the same object. n_strata 
            is the number of strata used in the populational phase of 
            AdaptiveMCMC (fixed phase), and strata_size is the number of 
            individuals per strata."""
        self.__init_iterations = init_iterations
        self.__sigma_update_n = sigma_update_n
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
    
        # Choose betas:
        betas = AdaptiveMCMC.__sample_betas (n_strata, strata_size)
        thetas = []
        cov_matrices = []

        
        for beta in betas:
            # First sampling step
            mcmc_init = MCMCInitialization (theta_prior, model, experiments,
                    self.__sigma_update_n, beta)
            start_sample = mcmc_init.get_sample (n_init)

            # Second sampling step
            amcmc = AdaptiveMCMC (model, experiments, start_sample, 
                    n_strata, strata_size, beta)
            theta, cov_matrix = amcmc.get_sample (n_adap, n_fixed)
            
            thetas.append (theta)
            cov_matrices.append (cov_matrix)
        
        # Third step 
        betas, thetas, likelihoods = self.__fixed_phase (self.__fixed_iterations, 
                betas, thetas, model, cov_matrices)

        ml = self.__calculate_marginal_likelihood (betas, thetas, \
                likelihoods)
        return ml


    def __get_proposal_dist (self, theta, var_covar):
        """ Returns the proposal distribution given a current point. 
            This distribution is MultivariateLognormal and its 
            parameters are set in such manner that the mean of this 
            distribution is exactly theta. """
        values = theta.get_values ()
        variances = var_covar.diagonal ()
        # The mean of the normal distribution should be 
        #      log (theta) - 1/2 diagonal (S)
        # to guarantee that the jump has theta as its mean value
        mean = np.log (values) - variances / 2
        jump_dist = MultivariateLognormal (mean, var_covar)
        return jump_dist

    def __propose_jump (self, theta, jump_dist):
        """ Proposes a jump from theta. This jump is log multivariate
            normal with covariance matrix self.__covar_matrix. """
        n = theta.get_size ()
        new_values = jump_dist.rvs ()
        new_t = theta.get_copy ()

        print ("\n\n\n------------------------------")
        for i in range (n):
            new_t[i].value = new_values[i]
        return new_t


    def __perform_jump (self, old_t, new_t, old_l, new_l, \
            old_proposal, var_covar, power=None):
        """ Given a current theta, decides to jump or not to new_t 
            according to Metropolis-Hastings. If the jump is accepted
            this function returns a tuple (new_t likelihood, 
            new_t jumping distribution); otherwise the method returns 
            False. """
        # Debugging #
        print ("\nCurrent theta: ", end='')
        for r in old_t:
            print (r.value, end=' ')
        print ("\nNew theta:     ", end='')
        for r in new_t:
            print (r.value, end=' ')
        print ("")
        print ("Current likelihood: " + str (old_l))
        print ("New likelihood: " + str (new_l), end='\n\n\n')
        # Debugging #

        if not new_l > 0:
            return False
        
        # Jumping probabilities
        new_proposal = self.__get_proposal_dist (new_t, var_covar)
        new_values = new_t.get_values ()
        old_values = old_t.get_values ()
        new_gv_old = old_proposal.pdf (new_values)
        old_gv_new = new_proposal.pdf (old_values)

        # p (old_t) and p (new_t)
        old_prior = old_t.get_p ()
        new_prior = new_t.get_p ()
        
        # ratio calculation
        l_ratio = new_l / old_l
        if power is not None:
            l_ratio = safe_power (l_ratio, power)
        theta_prior_ratio = new_prior / old_prior
        jump_ratio = old_gv_new / new_gv_old
        r = l_ratio * theta_prior_ratio  * jump_ratio
        
        if np.random.uniform () <= r:
            return (new_l, new_proposal)
        else:
            return False


    def __fixed_phase (self, N2, betas, theta_chains, model, var_covars):
        """ Performs the fixed phase of the algorithm. """
        experiments = self.__experiments
        n_strata = self.__n_strata
        strata_size = self.__strata_size
        likeli_f = LikelihoodFunction (model)

        pop_likelihoods = []
        for i in range (len (theta_chains)):
            theta = theta_chains[i]
            l = likeli_f.get_experiments_likelihood (experiments, theta)
            pop_likelihoods.append (l)
        
        for i in range (N2):
            # Local move
            for j in range (len (theta_chains)):
                var_covar_j = var_covars[j]
                old_t = theta_chains[j]
                old_proposal = self.__get_proposal_dist (old_t, var_covar_j)
                new_t = self.__propose_jump (old_t, old_proposal)
                old_l = likeli_f.get_experiments_likelihood \
                        (experiments, theta_chains[j])
                new_l = likeli_f.get_experiments_likelihood \
                        (experiments, new_t)
                result = self.__perform_jump (old_t, new_t, old_l, \
                        new_l, var_covar_j, betas[i])
                if result:
                    new_l = result[0]
                    theta_chains[j] = new_t
                    pop_likelihoods[j] = new_l
            
            # Global move
            j = np.random.choice (range (len (theta_chains) - 1))
            jump_distr = DiscreteLaplacian (len (betas), j + 1)
            k = jump_distr.rvs () - 1 
            thetaj = theta_chains[j]
            thetak = theta_chains[k]
            thetaj_l = likeli_f.get_experiments_likelihood (experiments, 
                    thetaj)
            thetak_l = likeli_f.get_experiments_likelihood (experiments, 
                    thetak)
            
            if thetaj_l == 0 or thetak_l == 0:
                continue

            inv_jump_dist = DiscreteLaplacian (len (betas), k + 1)
            j_gv_k = inv_jump_dist.pdf (j + 1)
            k_gv_j = jump_distr.pdf (k + 1)
            tjotk = thetaj_l / thetak_l
            tkotj = thetak_l / thetaj_l
            r = safe_power (tkotj, betas[j]) * \
                    safe_power (tjotk, betas[k]) * \
                    (j_gv_k / k_gv_j)
            
            if np.random.uniform () <= r:
                aux = theta_chains[j]
                theta_chains[j] = theta_chains[k]
                theta_chains[k] = aux
                aux = pop_likelihoods[j]
                pop_likelihoods[j] = pop_likelihoods[k]
                pop_likelihoods[k] = aux
        return (betas, theta_chains, pop_likelihoods)

    def __calculate_marginal_likelihood (self, betas, thetas, \
            likelihoods):
        """ Given a list with samples, calculates the marginal 
            likelihood. """
        ml = 0
        j = 0 
        n_strata = self.__n_strata
        strata_size = self.__strata_size
        sched_power = AdaptiveMCMC.get_sched_power ()
    
        print ("Estimating marginal likelihood")
        for i in range (n_strata):
            strata_start = (i / n_strata) ** sched_power
            strata_end = ((i + 1) / n_strata) ** sched_power
            del_strata = strata_end - strata_start

            print ("strata_start = " + str (strata_start))
            print ("strata_end = " + str (strata_end))
            print ("del_strata = " + str (del_strata))
                
            strat_sum = 0
            while j < len (betas) and betas[j] < strata_end:
                p_y_given_theta = likelihoods[j]
                print ("\tTheta: ", end='')
                for r in thetas[j]:
                    print (r.value, end=' ')
                print ("\n\tLikelihood: " + str (p_y_given_theta) + "\n")
                if p_y_given_theta > 0:
                    strat_sum += np.log (p_y_given_theta)
                print ("Strat_sum = " + str (strat_sum))
                j += 1

            strat_sum *= (del_strata / strata_size)
            ml += strat_sum
        print ("Calculated marginal likelihood: " + str (ml))
        return ml
