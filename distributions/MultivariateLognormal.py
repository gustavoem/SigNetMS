import sys
sys.path.insert (0, '..')

import numpy as np
from numpy.random import multivariate_normal as multivar_normal
from utils import safe_log

class MultivariateLognormal:
    """ This class implements a multivariate lognormal random variable. 
    """

    def __init__ (self, mu, Sigma):
        """ Default constructor. mu and sigma are the mean and variance
            covariance matrix, respectively, of the underlying 
            multivariate normal distribution. """
        self.__mu = np.array (mu)
        self.__S = np.array (Sigma)
        self._inv_S = None
        self._det_S = None


    def copy (self):
        """ Returns a copy of this object. """
        cpy = MultivariateLognormal (self.__mu, self.__S)
        cpy.set_S_inverse (self._inv_S)
        cpy.set_S_determinant (self._det_S)
        return cpy


    def set_S_inverse (self, inv_S):
        """ Sets the inverse of S, if you already have it. """
        self._inv_S = inv_S


    def set_S_determinant (self, det_S):
        """ Sets the determinant of S, if you already have it. """
        self._det_S = det_S


    def __get_S_inverse (self):
        """ Returns the inverse of the S matrix. """
        if self._inv_S is None:
            self._inv_S = np.linalg.inv (self.__S)
        return self._inv_S


    def __get_S_determinant (self):
        """ Returns the determinant of the S matrix. """
        if self._det_S is None:
            self._det_S = np.linalg.det (self.__S)
        return self._det_S


    def mean (self):
        """ Returns the mean of this random variable. """
        mu = self.__mu
        S = self.__S
        return np.exp (mu + S.diagonal () / 2)


    def rvs (self, n=None):
        """ Returns independent observations of this random variable. 
        """
        mu = self.__mu
        S = self.__S
        if n is None:
            normal_values = multivar_normal (mu, S)
        else:        
            normal_values = multivar_normal (mu, S, n)
        lognormal_values = np.exp (normal_values)
        return lognormal_values


    def pdf (self, x):
        """ Returns the value of the probability density function of 
            this random variable on point x. """
        # This is not very clean... we are assuming that the 
        # distribution is not concentrated around very small numbers.
        # What we should do instead is to continue calculations until we
        # can determine that the probability of the point is actually 
        # very small...
        if any (xi <= 1e-20 for xi in x):        
            return 0
            
        mu = self.__mu
        n = len (mu)
        inv_S = self.__get_S_inverse ()
        det_S = self.__get_S_determinant ()
        logx_minus_mu = np.log (x) - mu
        logx_minus_mu.shape = (n, 1)
        logx_minus_mu_t = logx_minus_mu.transpose ()
        
        term1 = 1 / np.sqrt (np.power (2 * np.pi, n) * abs (det_S))
        term2 = 1 
        for xi in x:
            term2 *= xi
        term2 = 1 / term2
        term3 = float (np.exp (-.5 * np.dot (np.dot (logx_minus_mu_t, 
            inv_S), logx_minus_mu)))
        
        return term1 * term2 * term3

    
    def log_pdf (self, x):
        """ Returns the value of the probability density function of 
            this random variable on point x. """
        # TODO: simplify calculations
        p = self.pdf (x)
        logp = safe_log (p)
        return logp


    @staticmethod
    def create_lognormal_with_shape (mu, S):
        """ Creates a Lognormal distribution with mean mu and 
            covariance S (diagnoal matrix). """
        mu = np.array (mu)
        S = np.array (S)
        S_diagonal = S.diagonal ()
        mu2 = mu * mu
        n = len (mu)
        n_ones = np.ones (n)
        
        # normal_mu_i = ln (mu_i / sqrt[S_ii / mu_i^2 + 1])
        normal_mu = np.log (mu / np.sqrt (S_diagonal / mu2 + n_ones))

        # normal_S_diag_i = ln (S_ii / mu_i^2 + 1)
        normal_S_diagonal = np.log (S_diagonal / (mu * mu) + n_ones)               
        normal_S = normal_S_diagonal * np.eye (n)
        return MultivariateLognormal (normal_mu, normal_S)
