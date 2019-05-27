import sys
sys.path.insert (0, '..')

import numpy as np
from numpy.random import multivariate_normal as np_multivar_normal
from utils import safe_log

class MultivariateNormal:
    """ This class implements a multivariate multivariate normal random 
        variable. """

    def __init__ (self, mu, Sigma):
        """ Default constructor. mu and sigma are the mean and variance
            covariance matrix, respectively, of the underlying 
            multivariate normal distribution. """
        self.__mu = np.array (mu)
        self.__S = np.array (Sigma)
        self.__inv_S = None
        self.__det_S = None

    
    def copy (self):
        """ Returns a copy of this object. """
        cpy = MultivariateNormal (self.__mu, self.__S)
        cpy.__inv_S = self.__inv_S
        cpy.__det_S = self.__det_S
        return cpy


    def __get_S_inverse (self):
        """ Returns the inverse of the S matrix. """
        if self.__inv_S is None:
            self.__inv_S = np.linalg.inv (self.__S)
        return self.__inv_S


    def __get_S_determinant (self):
        """ Returns the determinant of the S matrix. """
        if self.__det_S is None:
            self.__det_S = np.linalg.det (self.__S)
        return self.__det_S


    def mean (self):
        """ Returns the mean of this random variable. """
        return self.__mu


    def rvs (self, n=None):
        """ Returns independent observations of this random variable. 
        """
        mu = self.__mu
        S = self.__S
        if n is None:
            normal_values = np_multivar_normal (mu, S)
        else:        
            normal_values = np_multivar_normal (mu, S, n)
        return np.array (normal_values)


    def pdf (self, x):
        """ Returns the value of the probability density function of 
            this random variable on point x. """
        if any (xi <= 0 for xi in x):        
            return 0
            
        mu = self.__mu
        S = self.__S
        n = len (mu)
        inv_S = self.__get_S_inverse ()
        det_S = self.__get_S_determinant ()
        x_minus_mu = x - mu
        x_minus_mu.shape = (n, 1)
        x_minus_mu_t = x_minus_mu.transpose ()
        term1 = 1 / np.sqrt (np.power (2 * np.pi, n) * det_S)
        term2 = float (np.exp (-.5 * np.dot (np.dot (x_minus_mu_t, 
            inv_S), x_minus_mu)))
        return term1 * term2

    
    def log_pdf (self, x):
        """ Returns the value of the probability density function of 
            this random variable on point x. """
        # TODO: simplify calculations
        return safe_log (self.pdf (x))
