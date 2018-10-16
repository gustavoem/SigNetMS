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
            theta_values = theta.get_values ()
            sample_values.append (theta_values)
        return calc_covariance (sample_values)


    def __propose_jump (current_t):
        """ Proposes a jump on current_t. This jump is log multivariate
            normal with covariance matrix self.__covar_matrix. """
        n = current_t.get_size ()
        mean = np.zeros (n)
        cov = self.__covar_matrix
        jump = np.random.multivariate_normal (mean, cov)

    
    def __adapting_phase (self, N1):
        """ Performs the adaptive phase of the algorithm. """
        experiments = self.__experiments
        error_sigma = .5
        likeli_f = LikelihoodFunction (self.__model, error_sigma)

        current_t = self.__sampled_params[-1].get_copy ()
        current_l = likeli_f.get_experiments_likelihood (experiments, 
                current_t)

        for i in range (N1):
            print ("\nCurrent theta: ", end='')
            for r in current_t:
                print (r.value, end=' ')
            print ("\nNew theta:     ", end='')
            new_t = self.__propose_jump (current_t)
            for r in new_t:
                print (r.value, end=' ')
            
            new_l = likeli_f.get_experiments_likelihood (experiments, 
                    new_t)
            print ("Current likelihood: " + str (old_l))
            print ("New likelihood: " + str (new_l), end='\n\n\n')

            if new_l > 0:
                r = new_l / old_l
                if np.random.uniform () <= r:
                    accepted_jumps += 1
                    current_t = new_t
                    old_l = new_l
                    self.__sampled_params.append (current_t)
                    self.__covar_matrix = self.__estimate_cov_matrix ()

    
    def __fixed_phase (self, N2):
        """ Performs the fixed phase of the algorithm. """


    def get_sample (self, N1, N2):
        """ Returns a sample of theta. The sampling occurs in two 
            phases. The first is an adaptive MCMC, which has N1 
            iterations. The second, with a fixed proposal distribution 
            generated on the first phase, is a simple MCMC with N2 
            iterations. """
        self.__adapting_proposal (N1)
