import numpy as np
from LikelihoodFunction import LikelihoodFunction
from CovarianceMatrix import calc_covariance

class AdaptiveMCMC:
    """ This class receives a current sample from a target distribution
        on creation of an object and, usign an adaptive MCMC generates
        a new sample from the target distribution. The adaptive MCMC is
        described on "Supplementary Materials for Inferring Signaling 
        Pathway Topologies from Multiple Perturbation Measurements of 
        Specific Biochemical Species", Tian-Rui Xu et. al. Altough we
        run an adaptive MCMC, the returned sample is sampled with a 
        fixed proposal distribution, as recommended on "Bayesian Data
        Analysis", Gelman. """

    def __init__ (self, model, experiments, start_sample):
        self.__model = model
        self.__experiments = experiments
        self.__sampled_params = start_sample
        self.__covar_matrix = self.__estimate_cov_matrix ()
        print (self.__covar_matrix)
    
    
    def __estimate_cov_matrix (self):
        """ Given the current sample, estimates the covariance matrix
            of the model variables. """
        sample_values = []
        for theta in self.__sampled_params:
            theta_values = []
            for p in theta:
                theta_values.append (p.value)
            sample_values.append (theta_values)
        return calc_covariance (sample_values)
