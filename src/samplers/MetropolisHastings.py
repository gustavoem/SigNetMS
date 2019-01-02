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
        self.__sample_likelds = []
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


    def define_start_sample (self, sample, likelds):
        """ Inserts sampled parameters and its likelihoods at the 
            beginning of the self.__sample array. """
        if len (sample) != len (likelds):
            raise ValueError ("sample and likelds should have same " \
                    + "dimensions.")

        self.__sample = sample + self.__sample
        self.__sample_likelds = likelds + self.__sample_likelds


    def get_sample (self, N):
        """ Get a sample of size N. """
        # TODO


    def plot_sampled_distribution (self, N):
        """ Plots the last N sampled parameters. """
        # TODO 


    def manual_jump (self, theta, likeli):
        """ Manually jump from current theta to theta. """
        self.__sample.append (theta)
        self.__sample_likelds.append (likeli)
        self.__n_jumps += 1
        self.__n_accepted += 1


    def __calc_mh_ratio (self):
        """ Returns the ratio that is going to be used on the 
            Metropolis-Hastings algorithm as the probability of 
            accepting the proposed jump. """
        raise NotImplementedError
