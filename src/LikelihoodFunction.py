#The likelihood function
#follows the same definition as seen in Tian-Rui Xu, et. al, 
#Supplementary Materials for Inferring Signaling Pathway 
#Topologies from Multiple Perturbation Measurements of Specific 
#Biochemical Species

import numpy as np
from scipy.stats import lognorm
import RandomParameter

class LikelihoodFunction:
    """ This class defines a likelihood function for experimental data
        observed from dynamical systems that can be modeled by a sytem
        of ODEs. """ 

        
    def __init__ (self, ode):
        """ Default constructor. ode is the system that rules the 
            observed system. """
        self.__ode = ode
        
    def __point_likelihood (self, mu, x, sigma):
        exp = np.exp (-0.5 * ((x - mu) / sigma) ** 2)
        return exp * (1 / (sigma * np.sqrt (2 * np.pi)))


    def __calculate_likelihood (self, X_sys, X_obs, sigma):
        """ Calculates the likelihood of X_obs given that the real 
            solution i s X_sys. """
        likeli = 1
        for i in range (len (X_obs)):
            observed_x = X_obs[i]
            system_x = X_sys[i]
            likeli *= self.__point_likelihood (system_x, observed_x, 
                    sigma)

        return likeli


    def __get_system_state (self, var, t, theta):
        """ Calculates the state of the variable var on times t and with 
            theta parameters. theta should be a RandomParameter object. """
        if (theta is not None):
            for param in theta.get_model_parameters ():
                self.__ode.define_parameter (param.name, param.value)

        # TODO: the likelihood should be on a measurement of the system
        # and not directly on the state.
        # We could create this measure as a variable on the sbml model!
        # That would make the result of this function just the value of
        # this variable
        system_states = self.__ode.evaluate_on (t)

        # returns the state of a variable
        if var == "ERKPP":
            return np.asarray (system_states[var]) / 100
        else:
            return np.asarray (system_states[var])
    

    def get_experiment_likelihood (self, experiment, theta):
        """ Given an experiment, what is the probability that the values
            of this experiment were observed given that the system 
            random parameters are theta. Initial variable values are
            stored in the ode object. """
        t = experiment.times
        var = experiment.var
        X_sys = self.__get_system_state (var, t, theta)
        var = experiment.var
        X_obs = experiment.values
        sigma = theta.get_experimental_error ()
        return self.__calculate_likelihood (X_sys, X_obs, sigma)
        

    def get_experiments_likelihood (self, experiments, theta):
        """ Given a list of independent experiments that happens all 
            with the same time intervals and with respect to the same 
            variable, calculates the likelihood of all expeirments. """
        t = experiments[0].times
        var = experiments[0].var
        X_sys = self.__get_system_state (var, t, theta)
        sigma = theta.get_experimental_error ()
        # print ("X_sys: " + str (X_sys))
        l = 1
        for exp in experiments:
            X_obs = exp.values
            # print ("\tX_obs: " + str (X_obs))
            l *= self.__calculate_likelihood (X_sys, X_obs, sigma)
            # print ("\tpartial likelihood: " + str (l))
        return l
