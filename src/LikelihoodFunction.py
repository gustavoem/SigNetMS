#The likelihood function
#follows the same definition as seen in Tian-Rui Xu, et. al, 
#Supplementary Materials for Inferring Signaling Pathway 
#Topologies from Multiple Perturbation Measurements of Specific 
#Biochemical Species

import numpy as np

class LikelihoodFunction:
    """ This class defines a likelihood function for experimental data
        observed from dynamical systems that can be modeled by a sytem
        of ODEs. """ 

        
    def __init__ (self, ode, sigma):
        """ Default constructor. ode is the system that rules the 
            observed system. """
        self.__ode = ode
        self.__sigma = sigma
        
    def __point_likelihood (self, mu, x):
        exp = np.exp (-0.5 * ((x - mu) / self.__sigma) ** 2)
        return exp * (1 / (self.__sigma * np.sqrt (2 * np.pi)))

    
    def get_experiment_likelihood (self, experiment, theta):
        """ Given the observed X of values of variable var on time t, 
            what is the probability that X was observed given that the 
            system parameters are theta. Initial variable values are
            stored in the ode object. """
        if (theta is not None):
            for param in theta:
                self.__ode.define_parameter (param, theta[param])

        t = experiment.times
        var = experiment.var
        X = experiment.values
        system_states = self.__ode.evaluate_on (t)
        noiseless_X = system_states[var]
        likeli = 1
        for i in range (len (X)):
            observed_x = X[i]
            system_x = noiseless_X[i]
            likeli *= self.__point_likelihood (system_x, observed_x)
        return likeli
