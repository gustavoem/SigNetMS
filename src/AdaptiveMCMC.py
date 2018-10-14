import numpy as np
from LikelihoodFunction import LikelihoodFunction

class AdaptiveMCMC:
    """ This class performs MCMC sampling of theta with an adaptive
        proposal function. """

    def __init__ (self, param_list, model, experiments):
        """ Default constructor. """
        self.__params = param_list
        self.__model = model
        self.__experiments = experiments
        self.__start_params ()
        self.__sampled_params = []
        self.__sigma_update_n = 20


    def __start_params (self):
        """ Starts the parameter with a random value sampled from its
            prior. """
        for param in self.__params:
            param.set_rand_value ()


    def __propose_idependent_t (self, current_t, sigma):
        """ Given theta, propose a new theta, jumping each parameter
            indepenntly accorind to a lognormal distribution with
            mean zero and sqrt (variance) = sigma. """
        new_t = []
        for p in current_t:
            new_p = p.copy ()
            pjump = np.random.lognormal (0, sigma) - 1
            # print ("\nB = " + str (new_p.value))
            # print ("J = " + str (pjump))
            if new_p.value + pjump >= 0:
                new_p.value = pjump + new_p.value
            # print ("A = " + str (new_p.value))
            new_t.append (new_p)
        return new_t


    def __init_jump_sigma (self):
        """ Gives an initial value for the proposal distribution sigma
            on phase one of inferece. """
        params = self.__params
        param_mean = 0.0
        for p in params:
            param_mean += p.value
        param_mean /= len (params)
        sigma = np.random.gamma (2.0, param_mean)
        return sigma


    def __update_sigma (self, sigma, iterations, accepted_jumps):
        """ Updates sigma according to current acceptance rate. """
        acceptance_rate = accepted_jumps / iterations
        print ("Acceptance rate: " + str (acceptance_rate))
        if acceptance_rate > .4:
            sigma += sigma * .5
        if acceptance_rate < .25:
            if (sigma > 1e-7):
                sigma -= sigma * .5
        print ("New sigma: " + str (sigma))
        return sigma

    
    def __first_phase (self):
        jump_sigma = self.__init_jump_sigma ()
        error_sigma = .5
        experiments = self.__experiments
        likeli_f = LikelihoodFunction (self.__model, error_sigma)
        accepted_jumps = 0

        current_t = self.__params
        old_l = likeli_f.get_experiments_likelihood (experiments, 
                current_t)

        for i in range (100000):
            print ("Current theta: ", end='')
            for r in current_t:
                print (r.value, end=' ')
            print ("\nNew theta:     ", end='')
            new_t = self.__propose_idependent_t (current_t, jump_sigma)
            for r in new_t:
                print (r.value, end=' ')
            
            new_l = likeli_f.get_experiments_likelihood (experiments, 
                    new_t)

            print ("\nCurrent sigma: " + str (jump_sigma))
            print ("\nCurrent likelihood: " + str (old_l))
            print ("New likelihood: " + str (new_l), end='\n\n\n')

            if new_l > 0:
                r = new_l / old_l
                if np.random.uniform () <= r:
                    accepted_jumps += 1
                    current_t = new_t
                    old_l = new_l
                    self.__sampled_params.append (current_t)
            if (i + 1) % 20 == 0 : 
                jump_sigma = self.__update_sigma (jump_sigma, 
                        self.__sigma_update_n, accepted_jumps)
                accepted_jumps = 0

            
        
    def __second_phase (self):
        pass


    def get_sample (self, betas):
        self.__first_phase ()
        self.__second_phase ()
        return None
