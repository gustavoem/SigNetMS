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
from MultivariateLognormal import MultivariateLognormal
from DiscreteLaplacian import DiscreteLaplacian
from utils import safe_power
import numpy as np
import random
import re


class MarginalLikelihood:
    """ This class is able to perform an adaptive MCMC sampling to 
        estimate the likelihood of a model given experimental data. """
    
    __SCHEDULE_POWER = 4

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
    
    
    @staticmethod
    def __sample_betas (n_strata, strata_size):
        """ Samples strata_size beta from each of the n_strata 
            strata."""
        betas = []
        sched_power = MarginalLikelihood.__SCHEDULE_POWER
        for i in range (n_strata + 1):
            betas.append ((i / n_strata) ** sched_power)
        return betas

       
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
        betas = self.__sample_betas (n_strata, strata_size)
        thetas = []
        cov_matrices = []

        for beta in betas:
            # First sampling step
            print ("\nFIRST STEP | T = " + str (beta) + "\n")
            mcmc_init = MCMCInitialization (theta_prior, model, 
                    experiments, self.__sigma_update_n, beta)
            start_sample = mcmc_init.get_sample (n_init)

            # Second sampling step
            amcmc = AdaptiveMCMC (model, experiments, start_sample, 
                    n_strata, strata_size, beta)
            theta, cov_matrix = amcmc.get_sample (n_adap, n_fixed)
            
            thetas.append (theta)
            cov_matrices.append (cov_matrix)
        
        # Third step 
        betas, thetas, log_ls = self.__fixed_phase (self.__fixed_iterations, 
                betas, thetas, model, cov_matrices, experiments)

        ml = self.__calculate_marginal_likelihood (betas, thetas, 
                log_ls)
        return ml


    def __get_proposal_dist (self, theta, var_covar):
        """ Returns the proposal distribution given a current point. 
            This distribution is MultivariateLognormal and its 
            parameters are set in such manner that the mean of this 
            distribution is exactly theta. """
        values = theta.get_values ()
        var_covar = np.array (var_covar)
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

        for i in range (n):
            new_t[i].value = new_values[i]
        return new_t


    def __perform_jump (self, old_t, new_t, old_log_l, new_log_l, \
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
        print ("Current likelihood: " + str (old_log_l))
        print ("New likelihood: " + str (new_log_l), end='\n\n\n')
        # Debugging #

        if not new_log_l > float ("-inf"):
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
        l_ratio = np.exp (new_log_l - old_log_l)
        if power is not None:
            l_ratio = safe_power (l_ratio, power)
        theta_prior_ratio = new_prior / old_prior
        jump_ratio = old_gv_new / new_gv_old
        r = l_ratio * theta_prior_ratio  * jump_ratio
        
        if np.random.uniform () <= r:
            return (new_log_l, new_proposal)
        else:
            return False


    def __fixed_phase (self, N2, betas, theta_chains, model, var_covars,
            experiments):
        """ Performs the fixed phase of the algorithm. """
        n_strata = self.__n_strata
        strata_size = self.__strata_size
        l_f = LikelihoodFunction (model)

        pop_log_likelds = []
        sample_log_likelds = [[] for i in range (len (betas))]
        thetas_samples = [[] for i in range (len (betas))]
        for i in range (len (theta_chains)):
            theta = theta_chains[i]
            l = l_f.get_log_likelihood (experiments, theta)
            pop_log_likelds.append (l)
            sample_log_likelds[i].append(l)
            thetas_samples[i].append (theta)

        
        for i in range (N2):
            # Local move
            print ("\n\n\n------------------------------")
            print ("Iteration: " + str (i))
            for j in range (len (theta_chains)):
                var_covar_j = var_covars[j]
                old_t = thetas_samples[j][-1]
                old_log_l = sample_log_likelds[j][-1]
                old_proposal = self.__get_proposal_dist (old_t, var_covar_j)
                new_t = self.__propose_jump (old_t, old_proposal)
                new_log_l = l_f.get_log_likelihood (experiments, \
                        new_t)
                result = self.__perform_jump (old_t, new_t, old_log_l, \
                        new_log_l, old_proposal, var_covar_j, betas[j])
                if result:
                    new_log_l = result[0]
                    theta_chains[j] = new_t
                    pop_log_likelds[j] = new_log_l
                    sample_log_likelds[j].append (new_log_l)
                    thetas_samples[j].append (new_t)
                else:
                    sample_log_likelds[j].append (old_log_l)
                    thetas_samples[j].append (old_t)

            
            # Global move
            j = np.random.choice (range (len (theta_chains) - 1))
            jump_distr = DiscreteLaplacian (len (betas), j + 1)
            k = jump_distr.rvs () - 1 
            thetaj = thetas_samples[j][-1]
            thetak = thetas_samples[k][-1]
            thetaj_log_l = sample_log_likelds[j][-1]
            thetak_log_l = sample_log_likelds[k][-1]
            
            if not thetaj_log_l > float ("-inf") \
                    or not thetak_log_l > float ("-inf"):
                continue

            inv_jump_dist = DiscreteLaplacian (len (betas), k + 1)
            j_gv_k = inv_jump_dist.pdf (j + 1)
            k_gv_j = jump_distr.pdf (k + 1)
            tjotk = np.exp (thetaj_log_l - thetak_log_l)
            tkotj = np.exp (thetak_log_l - thetaj_log_l)

            r = safe_power (tkotj, betas[j]) * \
                    safe_power (tjotk, betas[k]) * \
                    (j_gv_k / k_gv_j)
            
            if np.random.uniform () <= r:
                aux = theta_chains[j]
                theta_chains[j] = theta_chains[k]
                theta_chains[k] = aux
                aux = pop_log_likelds[j]
                pop_log_likelds[j] = pop_log_likelds[k]
                pop_log_likelds[k] = aux
                sample_log_likelds[j].append (pop_log_likelds[j])
                sample_log_likelds[k].append (pop_log_likelds[k])
                thetas_samples[j].append (theta_chains[j])
                thetas_samples[k].append (theta_chains[k])

        return (betas, thetas_samples, sample_log_likelds)


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
        for i in range (len (thetas) - 1):
            t_i = betas[i] 
            t_ip1 = betas[i + 1] 
                
            exp_t_i = 0
            for log_l in likelihoods[i]:
                exp_t_i += log_l
            exp_t_i /= len (likelihoods[i])
            
            exp_t_ip1 = 0
            for log_l in likelihoods[i + 1]:
                exp_t_ip1 += log_l
            exp_t_ip1 /= len (likelihoods[i + 1])
            
            ml += (t_ip1 - t_i) * (exp_t_ip1 + exp_t_i) / 2
        print ("Calculated marginal likelihood: " + str (ml))
        return ml
