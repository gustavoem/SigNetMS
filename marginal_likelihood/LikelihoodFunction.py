# The likelihood function
# follows the same definition as seen in Tian-Rui Xu, et. al, 
# Supplementary Materials for Inferring Signaling Pathway 
# Topologies from Multiple Perturbation Measurements of Specific 
# Biochemical Species

import numpy as np
import math

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


    def __point_log_likelihood (self, mu, x, sigma):
        term1 = -0.5 * ((x - mu) / sigma) ** 2
        term2 = np.log (1 / (sigma * np.sqrt (2 * np.pi)))
        return term1 + term2


    def __calculate_likelihood (self, X_sys, X_obs, sigma):
        """ Calculates the log-likelihood of observing X_obs given that 
            the real value is X_sys. """
        log_l = 0
        for i in range (len (X_obs)):
            obs_xi = X_obs[i]
            sys_xi = X_sys[i]
            x_ll = self.__point_log_likelihood (sys_xi, obs_xi, sigma)
            log_l += x_ll    
        return log_l


    def __get_sys_measure (self, measure_expression, t, theta):
        """ Calculates the values of the measure on times t and with 
            theta parameters. theta should be a RandomParameter object. """
        if (theta is not None):
            for param in theta.get_model_parameters ():
                self.__ode.define_parameter (param.name, param.value)
    
        return self.__ode.evaluate_exp_on (measure_expression, t)


    def get_log_likelihood (self, experiments, theta):
        """ Given a list of independent experiments that happens all 
            with the same time intervals and with respect to the same 
            measure, calculates the likelihood of all expeirments. """
        t = experiments[0].times
        measure_expression = experiments[0].measure_expression
        X_sys = self.__get_sys_measure (measure_expression, t, theta)

        for sys_val in X_sys:
            if math.isnan (sys_val) or sys_val == float ("inf") or \
                    sys_val == float ("-inf"):
                return float ("-inf")

        sigma = theta.get_experimental_error ()
        #print ("\nX_sys: " + str (X_sys))
        log_l = 0
        for exp in experiments:
            X_obs = exp.values
            #print ("\tX_obs: " + str (X_obs))
            log_l += self.__calculate_likelihood (X_sys, X_obs, sigma)
            #print ("\tpartial log-likelihood: " + str (log_l))
        return log_l
