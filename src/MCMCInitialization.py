import numpy as np
from LikelihoodFunction import LikelihoodFunction

class MCMCInitialization:
    """ This class is able to return a sample of theta using an adaptive
        algorithm that samples each parameter independently and with
        a log-normal proposal distribution with adapting variance. This
        method is inspired in the one described on "Supplementary 
        Materials for Inferring Signaling Pathway Topologies from 
        Multiple Perturbation Measurements of Specific Biochemical 
        Species", Tian-Rui Xu et. al. """

    def __init__ (self, param_list, model, experiments):
        """ Default constructor. """
        self.__params = param_list
        self.__model = model
        self.__experiments = experiments
        self.__start_params ()
        self.__sampled_params = []
        self.__sigma_update_n = 100


    def __start_params (self):
        """ Starts the parameter with a random value sampled from its
            prior. """
        for param in self.__params:
            param.set_rand_value ()


    def __propose_independent_t (self, current_t, jump_sigma):
        """ Given theta, propose a new theta, jumping each parameter
            indepenntly accorind to a lognormal distribution with
            mean zero and sqrt (variance) = sigma. """
        new_t = current_t.get_copy ()
        for i in range (new_t.get_size ()):
            p = new_t[i]
            jump_s = jump_sigma[i]
            pjump = np.random.lognormal (0, jump_s) - 1
            if p.value + pjump >= 0:
                p.value += pjump
        return new_t


    def __init_jump_sigma (self):
        """ Gives an initial value for the proposal distribution sigma
            on phase one of inferece. """
        params = self.__params
        jump_sigma = []
        for p in params:
            sigma = np.random.gamma (2.0, p.value)
            jump_sigma.append (sigma)
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
        likeli_f = LikelihoodFunction (self.__model)
        accepted_jumps = 0

        current_t = self.__params
        old_l = likeli_f.get_experiments_likelihood (experiments, 
                current_t)

        for i in range (N):
            print ("\nCurrent theta: ", end='')
            for r in current_t:
                print (r.value, end=' ')
            print ("\nNew theta:     ", end='')
            new_t = self.__propose_independent_t (current_t, jump_sigma)
            for r in new_t:
                print (r.value, end=' ')
            
            new_l = likeli_f.get_experiments_likelihood (experiments, 
                    new_t)
            print ("\nCurrent sigma: " + str (jump_sigma))
            print ("Current likelihood: " + str (old_l))
            print ("New likelihood: " + str (new_l), end='\n\n\n')

            if new_l > 0:
                r = new_l / old_l
                if np.random.uniform () <= r:
                    accepted_jumps += 1
                    current_t = new_t
                    old_l = new_l
                    self.__sampled_params.append (current_t)
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

