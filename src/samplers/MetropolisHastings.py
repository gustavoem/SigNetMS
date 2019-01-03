import numpy as np

class MetropolisHastings:
    """ This class is an interface that should be used as base for 
        other classes that implement MetropolisHastings. We use the
        name theta for the object that is being sampled. We assume
        that there is some distribution that involves theta that is
        the target distribution. """
    
    def __init__ (self, theta):
        """ Default constructor. """
        self.__theta = theta
        self.__sample = []
        self.__sample_log_likelds = []
        self.__n_accepted = 0
        self.__n_jumps = 0
        self.__is_verbose = True
        
    
    def __create_jump_dist (self, theta_t):
        """ This method should return a distribution object (such as
            MultivariateLognormal) that is the distribution of the 
            random jump theta* that will be proposed given that the 
            current parameter is theta_t. """
        raise NotImplementedError


    def propose_jump (self, c_theta):
        """ Propose a new parameter theta* given that the current theta
            is c_theta. """
        proposal_distribution = self.__create_jump_dist (c_theta)
        new_theta = proposal_distribution.rvs ()
        return new_theta


    def get_acceptance_ratio (self):
        """ Returns the ratio  # accepted jumps / # jumps. """
        return self.__n_accepted / self.__n_jumps 


    def define_start_sample (self, sample, log_likelds):
        """ Inserts sampled parameters and its log-likelihoods at the 
            beginning of the self.__sample array. """
        if len (sample) != len (log_likelds):
            raise ValueError ("sample and log_likelds should have " \
                    + "same dimensions.")

        self.__sample = sample + self.__sample
        self.__sample_log_likelds = log_likelds + \
                self.__sample_log_likelds
    
    
    def start_sample_from_prior (self):
        """ Create a first sample based on the prior distribution of
            the parameter. """
        new_t = self.__theta.get_copy ()
        for p in new_t:
            val = p.set_rand_value ()
        new_l = self.__calc_log_likelihood (new_t)
        self.__sample.append (new_t)
        self.__sample_log_likelds.append (new_l)


    def manual_jump (self, theta, log_likeli):
        """ Manually jump from current theta to theta. If there's no
            current parameter, then theta becomes the first sample. """
        self.__sample.append (theta)
        self.__sample_log_likelds.append (log_likeli)
        self.__n_jumps += 1
        self.__n_accepted += 1


    def get_sample (self, N):
        """ Get a sample of size N. """
        if len (self.__sample[-1]) == 0:
            raise ValueError ("The sample can't be empty. Try using " \
                    + "the start_sample_from_prior () method.")

        for i in range (N):
            old_t = self.__sample[-1]
            old_l = self.__sample_log_likelds[-1]
            new_t = self.__propose_jump (old_t)
            new_l = self.__calc_log_likelihood (new_t)

            r = self.__calc_mh_ratio (old_t, old_l, new_t, new_l)
            if np.random.uniform () <= r:
                old_t = new_t
                old_l = new_l
            self.__sample.append (old_t)
            self.__sample_log_likelds.append (old_l)
        
        return get_last_sampled (N)
    

    def get_last_sampled (self, N):
        """ Returns the N last sampled parameters and a list of its
            log-likelihoods. """
        sample = []
        log_likelds = []
        i = -1
        while i < len (self.__)
        for t in self.__sample[-1:-(N + 1):-1]:
            sample.append (t.get_copy ())
            log_likelds.append (self.__sample_log_likelds[i])
            i -= 1
        return (sample, log_likelds)


    def __calc_mh_ratio (self, new_t, new_l, old_t, old_l):
        """ Returns the ratio that is going to be used on the 
            Metropolis-Hastings algorithm as the probability of 
            accepting the proposed jump new_t, given that the current
            parameter is old_t. """
        raise NotImplementedError


    def __calc_log_likelihood (self, theta):
        """ Should calculate the log-likelihood of a parameter theta. 
        """
        raise NotImplementedError


    def __iteration_update (self):
        """ Method called at the end of each iteration on get_sample.
        """
        continue
