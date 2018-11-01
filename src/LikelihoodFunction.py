#The likelihood function
#follows the same definition as seen in Tian-Rui Xu, et. al, 
#Supplementary Materials for Inferring Signaling Pathway 
#Topologies from Multiple Perturbation Measurements of Specific 
#Biochemical Species

import numpy as np
from scipy.stats import lognorm
import RandomParameter
from asteval import Interpreter

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


    def __get_sys_measure (self, measure_expression, t, theta):
        """ Calculates the values of the measure on times t and with 
            theta parameters. theta should be a RandomParameter object. """
        if (theta is not None):
            for param in theta.get_model_parameters ():
                self.__ode.define_parameter (param.name, param.value)

        system_states = self.__ode.evaluate_on (t)
        aeval = Interpreter ()

        sys_measures = []
        for i in range (len (t)):
            for var in system_states:
                aeval.symtable[var] = system_states[var][i]
            measure_value = aeval (measure_expression)
            sys_measures.append (measure_value)
            
        return sys_measures 

    def get_experiment_likelihood (self, experiment, theta):
        """ Given an experiment, what is the probability that the values
            of this experiment were observed given that the system 
            random parameters are theta. Initial variable values are
            stored in the ode object. """
        t = experiment.times
        measure_expression = experiment.measure_expression
        X_sys = self.__get_sys_measure (measure_expression, t, theta)
        X_obs = experiment.values
        sigma = theta.get_experimental_error ()
        return self.__calculate_likelihood (X_sys, X_obs, sigma)
        

    def get_experiments_likelihood (self, experiments, theta):
        """ Given a list of independent experiments that happens all 
            with the same time intervals and with respect to the same 
            measure, calculates the likelihood of all expeirments. """
        t = experiments[0].times
        measure_expression = experiments[0].measure_expression
        X_sys = self.__get_sys_measure (measure_expression, t, theta)
        sigma = theta.get_experimental_error ()
        # print ("\nX_sys: " + str (X_sys))
        l = 1
        for exp in experiments:
            X_obs = exp.values
            # print ("\tX_obs: " + str (X_obs))
            l *= self.__calculate_likelihood (X_sys, X_obs, sigma)
            # print ("\tpartial likelihood: " + str (l))
        return l
