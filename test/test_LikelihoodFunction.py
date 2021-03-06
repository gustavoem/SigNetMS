import sys
sys.path.insert(0, '..')

import unittest
import math
import numpy as np
from model.ODES import ODES
from marginal_likelihood.LikelihoodFunction import LikelihoodFunction
from distributions.Gamma import Gamma
from experiment.Experiment import Experiment
from model.RandomParameterList import RandomParameterList
from model.RandomParameter import RandomParameter

class TestLikelihoodFunction (unittest.TestCase):
    
    def setUp (self):
        # dx1 (t)/dt = x1 (t), x1 (0) = 1
        # x1 (t) = e ^ t
        self.odes = ODES ()
        self.odes.add_equation ("x1", "x1")
        self.odes.define_initial_value ("x1", 1.0)
        self.theta = RandomParameterList ()
        sigma_dist = Gamma (1, 1)
        sigma = RandomParameter ("sigma", sigma_dist)
        sigma.value = 1.0
        self.theta.set_experimental_error (sigma)
        


    def __gaussian (self, mu, sigma, x):
        """ Calculates a point of the gaussian function. """
        l = np.exp (-.5 * ((x - mu) / sigma) ** 2) * \
                (1 / (sigma * np.sqrt (2 * np.pi)))
        return l

        
    def test_get_likelihood_point (self):
        """ Tests if the likelihood can return the correct value when
            there's only one observation. """
        # x1(0) = 1.0
        # D ~ Gaussian (x1(0), 1) ~ Gaussian (1, 1)
        # f_D (1) = e ^ -{[(0) ^ 2] / [2 * 1]} * {1 * sqrt (2pi)} ^ -1
        f_D = self.__gaussian (1, 1, 1)
        analytic = np.log (f_D)
        
        t = [0]
        values = [1.0]
        var = "x1"
        experiments = [Experiment (t, values, var)]
        likelihood_f = LikelihoodFunction (self.odes)
        l = likelihood_f.get_log_likelihood (experiments, self.theta)
        assert (abs (analytic - l) < 1e-8)


    def test_get_likelihood_over_time (self):
        """ Tests if the likelihood can return the correct value when 
            the observation occurs in multiple time points. """
        t = [0, .25, .5, .75, 1]
        D = [np.exp (x) for x in t]
        experiments = [Experiment (t, D, "x1")]
        f_D = 1
        for y in D:
            f_D *= self.__gaussian (y, 1, y)
        analytic = np.log (f_D)
        
        likelihood_f = LikelihoodFunction (self.odes) 
        l = likelihood_f.get_log_likelihood (experiments, self.theta)
        assert (abs (analytic - l) < 1e-2)

    
    def test_get_likelihood_experiment_set (self):
        """ Tests if can return the likelihood of data in respect to
            multiple experimnets. """
        t = [0, .25, .5, .75, 1]
        D = [np.exp (x) for x in t]
        experiment = Experiment (t, D, "x1")
        experiments = [experiment, experiment]
        f_D = 1
        for y in D:
            f_D *= self.__gaussian (y, 1, y)
        f_D **= 2
        analytic = np.log (f_D)
        
        likelihood_f = LikelihoodFunction (self.odes)
        l = likelihood_f.get_log_likelihood (experiments, self.theta)
        assert (abs (analytic - l) < 1e-2)




