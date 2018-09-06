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


    def __calculate_likelihood (self, X_sys, X_obs):
        """ Calculates the likelihood of X_obs given that the real 
            solution i s X_sys. """
        likeli = 1
        for i in range (len (X_obs)):
            observed_x = X_obs[i]
            system_x = X_sys[i]
            likeli *= self.__point_likelihood (system_x, observed_x)
        return likeli


    def __get_system_state (self, var, t, theta):
        """ Calculates the state of the variable var on times t and with 
            theta parameters. """
        if (theta is not None):
            for param in theta:
                self.__ode.define_parameter (param, theta[param])

        system_states = self.__ode.evaluate_on (t)
        X_sys = system_states[var]
        return X_sys
    

    def get_experiment_likelihood (self, experiment, theta):
        """ Given the observed X of values of variable var on time t, 
            what is the probability that X was observed given that the 
            system parameters are theta. Initial variable values are
            stored in the ode object. """
        t = experiment.times
        var = experiment.var
        X_sys = self.__get_system_state (var, t, theta)
        var = experiment.var
        X_obs = experiment.values
        return self.__calculate_likelihood (X_sys, X_obs)
        

    def get_experiments_likelihood (self, experiments, theta):
        """ Given a list of experiments that happens all with the same
            time intervals and with respect to the same variable, 
            calculates the likelihood of all expeirments. """
        t = experiments[0].times
        var = experiments[0].var
        X_sys = self.__get_system_state (var, t, theta)
        l = 1
        for exp in experiments:
            X_obs = exp.values
            l *= self.__calculate_likelihood (X_sys, X_obs)
        return l
