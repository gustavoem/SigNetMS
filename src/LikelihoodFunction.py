#The likelihood function
#follows the same definition as seen in Tian-Rui Xu, et. al, 
#Supplementary Materials for Inferring Signaling Pathway 
#Topologies from Multiple Perturbation Measurements of Specific 
#Biochemical Species

import numpy as np
import RandomParameter

# TODO: theta should be RandomParameter

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
            theta parameters. theta should be a RandomParameter object. """
        if (theta is not None):
            for param in theta:
                self.__ode.define_parameter (param.name, param.value)

        # TODO: the likelihood should be on a measurement of the system
        # and not directly on the state.
        # We could create this measure as a variable on the sbml model!
        # That would make the result of this function just the value of
        # this variable
        system_states = self.__ode.evaluate_on (t)

        # returns the state of a variable
        # return system_states[var]

        # This is temporary.
        X_sys = []
        for i in range (len (system_states["ERK"])):
            ratio = system_states["ERKPP"][i] / (system_states["ERK"][i] + system_states["ERKPP"][i])
            X_sys.append (ratio * 100)
        return X_sys
    

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
        return self.__calculate_likelihood (X_sys, X_obs)
        

    def get_experiments_likelihood (self, experiments, theta):
        """ Given a list of independent experiments that happens all 
            with the same time intervals and with respect to the same 
            variable, calculates the likelihood of all expeirments. """
        t = experiments[0].times
        var = experiments[0].var
        X_sys = self.__get_system_state (var, t, theta)
        print ("X_sys: " + str (X_sys))
        l = 1
        for exp in experiments:
            X_obs = exp.values
            print ("\tX_obs: " + str (X_obs))
            l *= self.__calculate_likelihood (X_sys, X_obs)
            print ("\tpartial likelihood: " + str (l))
        return l
