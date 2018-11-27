import numpy as np
from LikelihoodFunction import LikelihoodFunction
from CovarianceMatrix import calc_covariance
from MultivariateLognormal import MultivariateLognormal
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


    def __propose_jump (self, current_t):
        """ Proposes a jump on current_t. This jump is log multivariate
            normal with covariance matrix self.__covar_matrix. """
        n = current_t.get_size ()
        cov = self.__covar_matrix
        # The mean of the normal distribution should be 
        #      log (current_t) + 1/2 diagonal (S)
        # to guarantee that the jump has current_t as its mean value
        current_t_values = np.array (current_t.get_values ())
        variances = cov.diagonal ()
        mean = np.log (current_t_values) - variances / 2
        jump_dist = MultivariateLognormal (mean, cov)
        new_values = jump_dist.rvs ()
        new_t = current_t.get_copy ()
        for i in range (n):
            new_t[i].value = new_values[i]
        return new_t

    
    def __adapting_phase (self, N1):
        """ Performs the adaptive phase of the algorithm. """
        experiments = self.__experiments
        likeli_f = LikelihoodFunction (self.__model)

        current_t = self.__sampled_params[-1].get_copy ()
        current_l = likeli_f.get_experiments_likelihood (experiments, 
                current_t)
        mean = current_t.get_values ()
        covar = self.__covar_matrix
        current_jump_dist = MultivariateLognormal (mean, covar)

        for i in range (N1):
            new_t = self.__propose_jump (current_t)
            new_l = likeli_f.get_experiments_likelihood (experiments, 
                    new_t)
            
            # Debugging #
            print ("\nCurrent theta: ", end='')
            for r in current_t:
                print (r.value, end=' ')
            print ("\nNew theta:     ", end='')
            for r in new_t:
                print (r.value, end=' ')
            print ("")
            print ("Current likelihood: " + str (current_l))
            print ("New likelihood: " + str (new_l), end='\n\n\n')
            # Debugging #

            if new_l > 0:
                mean = new_t.get_values ()
                covar = self.__covar_matrix
                new_jump_dist = MultivariateLognormal (mean, covar)
                
                new_values = new_t.get_values ()
                current_values = current_t.get_values ()
                new_gv_old = current_jump_dist.pdf (new_values)
                old_gv_new = new_jump_dist.pdf (current_values)
                
                old_prior = current_t.get_p ()
                new_prior = new_t.get_p ()

                r = (new_l / current_l) *  (new_prior / old_prior) * \
                        (old_gv_new / new_gv_old)
                        
                if np.random.uniform () <= r:
                    current_t = new_t
                    current_l = new_l
                    self.__sampled_params.append (current_t)
                    self.__covar_matrix = self.__estimate_cov_matrix ()


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
        theta_mean = self.__sampled_params[-1].get_copy ()
        for p in theta_mean:
            p.value = 0
        for t in self.__sampled_params:
            for i in range (t.get_size ()):
                theta_mean[i].value += t[i].value
        for p in theta_mean:
            p.value /= len (self.__sampled_params)
        
        theta_pop = []
        for b in betas:
            theta = self.__propose_jump (theta_mean)
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
        
        for i in range (N2):
            # Local move
            j = np.random.choice (range (len (theta_chains)))
            new_t = self.__propose_jump (theta_chains[j])
            old_l = likeli_f.get_experiments_likelihood (experiments, 
                    theta_chains[j])
            new_l = likeli_f.get_experiments_likelihood (experiments, 
                    new_t)
            
            print ("\nCurrent theta: ", end='')
            for r in theta_chains[j]:
                print (r.value, end=' ')
            print ("\nNew theta:     ", end='')
            for r in new_t:
                print (r.value, end=' ')
            print ("")
            print ("Current likelihood: " + str (old_l))
            print ("New likelihood: " + str (new_l), end='\n\n\n')
            
            if old_l != 0:
                r = (new_l / old_l) ** betas[j]

            if old_l == 0 or np.random.uniform () <= r:
                theta_chains[j] = new_t
            
            # Global move
            j = np.random.choice (range (len (theta_chains) - 1))
            theta1 = theta_chains[j]
            theta2 = theta_chains[j + 1]
            theta1_l = likeli_f.get_experiments_likelihood (experiments, 
                    theta1)
            theta2_l = likeli_f.get_experiments_likelihood (experiments, 
                    theta2)
            
            if theta1_l == 0 or theta2_l == 0:
                continue

            t1ot2 = theta1_l / theta2_l
            t2ot1 = theta2_l / theta1_l
            r = safe_power (t2ot1, betas[j]) * \
                    safe_power (t1ot2, betas[j + 1])
            if (np.random.uniform () <= r):
                aux = theta_chains[j]
                theta_chains[j] = theta_chains[j + 1]
                theta_chains[j + 1] = aux
        
        # print posterior
        figname = "posterior_after_third_phase_" + \
                theta_chains[0][0].name + ".png"
        plot_theta_var_sample (theta_chains, 0, figname)
        
        return (betas, theta_chains)


    def get_sample (self, N1, N2):
        """ Returns a sample of theta. The sampling occurs in two 
            phases. The first is an adaptive MCMC, which has N1 
            iterations. The second, with a fixed proposal distribution 
            generated on the first phase, is a simple MCMC with N2 
            iterations. """
        print ("ADAPTING PHASE")
        self.__adapting_phase (N1)
        
        # print posterior
        figname = "posterior_after_second_phase_" + \
                self.__sampled_params[0][0].name + ".png"
        plot_theta_var_sample (self.__sampled_params, 0, figname)

        
        print ("FIXED PHASE")
        return self.__fixed_phase (N2)
