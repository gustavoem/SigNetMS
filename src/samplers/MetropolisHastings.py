import numpy as np

class MetropolisHastings:
    """ This class is an interface that should be used as base for 
        other classes that implement MetropolisHastings. We use the
        name theta for the object that is being sampled. We assume
        that there is some distribution that involves theta that is
        the target distribution. """
    
    def __init__ (self, theta, verbose=False):
        """ Default constructor. """
        self._theta = theta.get_copy ()
        self.__sample = []
        self.__sample_log_likelds = []
        self._n_accepted = 0
        self._n_jumps = 0
        self._is_verbose = verbose
        
    
    def _create_jump_dist (self, theta_t):
        """ This method should return a distribution object (such as
            MultivariateLognormal) that is the distribution of the 
            random jump theta* that will be proposed given that the 
            current parameter is theta_t. """
        raise NotImplementedError


    def propose_jump (self, c_theta):
        """ Propose a new parameter theta* given that the current theta
            is c_theta. """
        proposal_distribution = self._create_jump_dist (c_theta)
        new_theta_values = proposal_distribution.rvs ()
        new_theta = c_theta.get_copy ()
        for i in range (len (new_theta_values)):
            p = new_theta[i] 
            p.value = new_theta_values[i]
        return new_theta


    def get_acceptance_ratio (self):
        """ Returns the ratio  # accepted jumps / # jumps. """
        return self._n_accepted / self._n_jumps 


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
        new_t = self._theta.get_copy ()
        for p in new_t:
            val = p.set_rand_value ()
        new_l = self._calc_log_likelihood (new_t)
        self.__sample.append (new_t)
        self.__sample_log_likelds.append (new_l)


    def manual_jump (self, theta, log_likeli):
        """ Manually jump from current theta to theta. If there's no
            current parameter, then theta becomes the first sample. """
        self.__sample.append (theta)
        self.__sample_log_likelds.append (log_likeli)
        self._n_jumps += 1
        self._n_accepted += 1


    def get_sample (self, N):
        """ Get a sample of size N. """
        if len (self.__sample) == 0:
            raise ValueError ("The current sample can't be empty. " \
                    + "Try using the start_sample_from_prior () " \
                    + "method.")

        for i in range (N):
            old_t = self.__sample[-1]
            old_l = self.__sample_log_likelds[-1]
            new_t = self.propose_jump (old_t)
            new_l = self._calc_log_likelihood (new_t)

            if self._is_verbose:
                print ("old_t: ")
                print ("[", end='')
                for p in old_t:
                    print (p.value, end=' ')
                print ("]")
                print ("new_t")
                print ("[", end='')
                for p in new_t:
                    print (p.value, end=' ')
                print ("]")
                print ("old_l: " + str (old_l))
                print ("new_l: " + str (new_l))


            r = self._calc_mh_ratio (new_t, new_l, old_t, old_l)
            if np.random.uniform () <= r:
                old_t = new_t
                old_l = new_l
                self._n_accepted += 1
            self.__sample.append (old_t)
            self.__sample_log_likelds.append (old_l)
            self._n_jumps += 1
            self._iteration_update ()
        
        return self.get_last_sampled (N)
    

    def get_last_sampled (self, N):
        """ Returns the N last sampled parameters and a list of its
            log-likelihoods. """
        sample = []
        log_likelds = []
        i = -1
        for t in self.__sample[-N:]:
            sample.append (t.get_copy ())
            log_likelds.append (self.__sample_log_likelds[i])
            i -= 1
        return (sample, log_likelds)


    def _calc_mh_ratio (self, new_t, new_l, old_t, old_l):
        """ Returns the ratio that is going to be used on the 
            Metropolis-Hastings algorithm as the probability of 
            accepting the proposed jump new_t, given that the current
            parameter is old_t. new_l and old_l should be the log
            likelihood of new_t and old_t, respectively. """
        raise NotImplementedError


    def _calc_log_likelihood (self, theta):
        """ Should calculate the log-likelihood of a parameter theta. 
        """
        raise NotImplementedError


    def _iteration_update (self):
        """ Method called at the end of each iteration on get_sample.
        """
        pass
