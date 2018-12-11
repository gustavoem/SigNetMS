import numpy as np
from LikelihoodFunction import LikelihoodFunction
from CovarianceMatrix import calc_covariance
from MultivariateLognormal import MultivariateLognormal
from DiscreteLaplacian import DiscreteLaplacian
from utils import safe_power
from utils import plot_theta_var_sample

import numpy as np

class AdaptiveMCMC:
    """ This class receives a current sample from a target distribution
        on creation of an object and, usign an adaptive MCMC generates
        a new sample from the target distribution. The adaptive MCMC is
        described on "Supplementary Materials for Inferring Signaling 
        Pathway Topologies from Multiple Perturbation Measurements of 
        Specific Biochemical Species", Tian-Rui Xu et. al. Altough we
        run an adaptive MCMC, the returned sample is sampled with a 
        fixed proposal distribution, as recommended on "Bayesian Data
        Analysis", Gelman. """


    # This is a good scheduling power according to "Estimating Bayes 
    # Factors via Thermodynamic Integration and Population MCMC"
    __SCHEDULE_POWER = 5


    def __init__ (self, model, experiments, start_sample, n_strata, 
            strata_size):
        self.__model = model
        self.__experiments = experiments
        self.__sampled_params = start_sample
        self.__covar_matrix = self.__estimate_cov_matrix ()
        self.__n_strata = n_strata
        self.__strata_size = strata_size

        self.__sampled_params = start_sample[len (start_sample) // 2:]
        print (self.__covar_matrix)
    

    def __estimate_cov_matrix (self):
        """ Given the current sample, estimates the covariance matrix
            of the model variables. """
        sample_values = []
        for theta in self.__sampled_params:
            theta_values = theta.get_values ()
            sample_values.append (theta_values)
        return calc_covariance (sample_values)
    
    
    def __get_proposal_dist (self, theta):
        """ Returns the proposal distribution given a current point. 
            This distribution is MultivariateLognormal and its 
            parameters are set in such manner that the mean of this 
            distribution is exactly theta. """
        values = theta.get_values ()
        var_covar = self.__covar_matrix
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
            old_proposal, power=None):
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
        new_proposal = self.__get_proposal_dist (new_t)
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

    
    def __adapting_phase (self, N1):
        """ Performs the adaptive phase of the algorithm. """
        experiments = self.__experiments
        likeli_f = LikelihoodFunction (self.__model)

        current_t = self.__sampled_params[-1].get_copy ()
        current_l = likeli_f.get_experiments_likelihood (experiments, 
                current_t)
        current_proposal = self.__get_proposal_dist (current_t)

        for i in range (N1):
            new_t = self.__propose_jump (current_t, current_proposal)
            new_l = likeli_f.get_experiments_likelihood (experiments, \
                    new_t)
            result = self.__perform_jump (current_t, new_t, current_l, \
                    new_l, current_proposal)
            if result:
                (new_likelihood, new_proposal) = result
                self.__sampled_params.append (new_t)
                self.__covar_matrix = self.__estimate_cov_matrix ()
                current_t = new_t
                current_l = new_l
                current_proposal = new_proposal


    @staticmethod
    def get_sched_power ():
        return AdaptiveMCMC.__SCHEDULE_POWER;


    @staticmethod
    def __sample_betas (n_strata, strata_size):
        """ Samples strata_size beta from each of the n_strata 
            strata."""
        betas = []
        sched_power = AdaptiveMCMC.__SCHEDULE_POWER
        for i in range (n_strata):
            strata_start = (i /  n_strata) ** sched_power
            strata_end = ((i + 1) / n_strata) ** sched_power
            for j in range (strata_size):
                x = np.random.uniform (strata_start, strata_end)
                betas.append (x)
            betas.sort ()
        return betas


    def __init_population (self, betas):
        """ Sample from the identified posterior. """
        n = self.__sampled_params[-1].get_size ()
        theta_pop = []
        for b in betas:
            theta = self.__sampled_params[-1].get_copy ()
            theta_pop.append (theta)
        return theta_pop

    
    def __fixed_phase (self, N2):
        """ Performs the fixed phase of the algorithm. """
        experiments = self.__experiments
        n_strata = self.__n_strata
        strata_size = self.__strata_size
        betas = AdaptiveMCMC.__sample_betas (n_strata, strata_size)
        theta_chains = self.__init_population (betas)
        likeli_f = LikelihoodFunction (self.__model)

        pop_likelihoods = []
        for i in range (len (theta_chains)):
            theta = theta_chains[i]
            l = likeli_f.get_experiments_likelihood (experiments, theta)
            pop_likelihoods.append (l)
        
        for i in range (N2):
            # Local move
            for j in range (len (theta_chains)):
                old_t = theta_chains[j]
                old_proposal = self.__get_proposal_dist (old_t)
                new_t = self.__propose_jump (old_t, old_proposal)
                old_l = likeli_f.get_experiments_likelihood \
                        (experiments, theta_chains[j])
                new_l = likeli_f.get_experiments_likelihood \
                        (experiments, new_t)
                result = self.__perform_jump (old_t, new_t, old_l, \
                        new_l, old_proposal)
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


    def get_sample (self, N1, N2):
        """ Returns a sample of theta. The sampling occurs in two 
            phases. The first is an adaptive MCMC, which has N1 
            iterations. The second, with a fixed proposal distribution 
            generated on the first phase, is a simple MCMC with N2 
            iterations. """
        print ("ADAPTING PHASE")
        self.__adapting_phase (N1) 
        
        print ("FIXED PHASE")
        return self.__fixed_phase (N2)
