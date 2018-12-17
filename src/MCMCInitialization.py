import sys
sys.path.insert (0, './distributions/')

import numpy as np
from LikelihoodFunction import LikelihoodFunction
from Lognormal import Lognormal
from utils import safe_power


class MCMCInitialization:
    """ This class is able to return a sample of theta using an adaptive
        algorithm that samples each parameter independently and with
        a log-normal proposal distribution with adapting variance. This
        method is inspired in the one described on "Supplementary 
        Materials for Inferring Signaling Pathway Topologies from 
        Multiple Perturbation Measurements of Specific Biochemical 
        Species", Tian-Rui Xu et. al. """

    def __init__ (self, param_list, model, experiments, sigma_update_n, beta):
        """ Default constructor. """
        self.__params = param_list
        self.__model = model
        self.__experiments = experiments
        self.__start_params ()
        self.__sampled_params = []
        self.__sigma_update_n = sigma_update_n
        self.__beta = beta


    def __start_params (self):
        """ Starts the parameter with a random value sampled from its
            prior. """
        for param in self.__params:
            param.set_rand_value ()


    def __propose_independent_t (self, old_t, jump_sigma):
        """ Given theta, propose a new theta, jumping each parameter
            indepenntly accorind to a lognormal distribution with
            mean zero and sqrt (variance) = sigma. Returns a tuple 
            containing a new list of parameter and its likelihood."""
        new_t = old_t.get_copy ()
        for i in range (new_t.get_size ()):
            p = new_t[i]
            jump_s = jump_sigma[i]
            normal_mean = np.log (p.value) - jump_s * jump_s / 2
            # If X is lognormal (mu, sigma), then E[X] is 
            #     exp (mu + sigma * sigma / 2)
            # If we set X as lognormal (log (mu) - sigma*sigma/2, sigma)
            # then E[X] is 
            #     exp (log (mu)) = mu
            jump_dist = Lognormal (normal_mean, jump_s)
            new_value = jump_dist.rvs ()
            p.value = new_value
        return new_t


    def __calc_jump_prob (self, t1, t2, jump_sigma):
        """ Calculates the probability of jumping from one set of 
            parameters to another. That is, returns 
            P (theta^t = t2 | theta^{t-1} = t1). """
        p = 1
        for i in range (t1.get_size ()):
            t1_param = t1[i].value
            t2_param = t2[i].value
            s = jump_sigma[i]
            mean = np.log (t1_param) - s * s / 2
            jump_dist = Lognormal (mean, s)
            p *= jump_dist.pdf (t2_param)
        return p


    def __init_jump_sigma (self):
        """ Gives an initial value for the proposal distribution sigma
            on phase one of inferece. """
        params = self.__params
        jump_sigma = []
        for p in params:
            param_dist = p.get_distribution ()
            prior_variance = param_dist.variance ()
            sigma2 = np.log (np.sqrt (prior_variance) + 1)
            jump_sigma.append (np.sqrt (sigma2))
        return jump_sigma


    def __update_sigma (self, jump_sigma, iterations, accepted_jumps):
        """ Updates sigma according to current acceptance rate. """
        acceptance_rate = accepted_jumps / iterations
        for i in range (len (jump_sigma)):
            if acceptance_rate > .4 and jump_sigma[i] < 1e8:
                jump_sigma[i] += jump_sigma[i] * .5
            if acceptance_rate < .25 and jump_sigma[i] > 1e-4:
                jump_sigma[i] -= jump_sigma[i] * .5
        return jump_sigma

    
    def __sample_theta (self, N):
        jump_sigma = self.__init_jump_sigma ()
        experiments = self.__experiments
        l_f = LikelihoodFunction (self.__model)
        accepted_jumps = 0

        old_t = self.__params
        old_log_l = l_f.get_log_likelihood (experiments, old_t)

        for i in range (N):
            new_t = self.__propose_independent_t (old_t, jump_sigma)
            new_log_l = l_f.get_log_likelihood (experiments, new_t)
            
            # Debugging #
            if i + 1 % 1000 == 0:
                print ("Iteration " + str (i + 1) + \
                        " on MCMCInitialization.")
            print ("\nCurrent theta: ", end='')
            for r in old_t[:10]:
                print (r.value, end=' ')
            print ("\nNew theta:     ", end='')
            for r in new_t[:10]:
                print (r.value, end=' ')
            print ("\nCurrent sigma: " + str (jump_sigma[:10]))
            print ("Current likelihood: " + str (old_log_l))
            print ("New likelihood: " + str (new_log_l), end='\n')
            # Debugging #

            minus_inf = float ("-inf")
            if new_log_l > minus_inf:
                new_gv_old = self.__calc_jump_prob (old_t, new_t, 
                        jump_sigma)
                old_gv_new = self.__calc_jump_prob (new_t, old_t, 
                        jump_sigma)
                new_prior = new_t.get_p ()
                old_prior = old_t.get_p ()
                l_ratio = safe_power (np.exp (new_log_l - old_log_l), self.__beta )
                print ("Likelihood ratio: " + str (l_ratio))
                prior_ratio = new_prior / old_prior
                jump_ratio = old_gv_new / new_gv_old
                r = l_ratio * prior_ratio * jump_ratio
                
                if not old_log_l > minus_inf \
                        or np.random.uniform () <= r:
                    accepted_jumps += 1
                    old_t = new_t
                    old_log_l = new_log_l
                    self.__sampled_params.append (old_t)

            if (i + 1) % self.__sigma_update_n == 0 : 
                jump_sigma = self.__update_sigma (jump_sigma, 
                        self.__sigma_update_n + 1, accepted_jumps)
                accepted_jumps = 0


    def get_sample (self, N):
        """ Samples theta with an Adaptive MCMC. The proposal 
            distribution on this method is log normal and has an 
            adaptive sigma value. """
        self.__sample_theta (N)
        sample = []
        for t in self.__sampled_params:
            t_cpy = t.get_copy ()
            sample.append (t_cpy)
        return sample

