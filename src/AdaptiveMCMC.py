import numpy as np
from LikelihoodFunction import LikelihoodFunction

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
        self.__sampled_params = []
        self.__covar_matrix = []
    
    

        



