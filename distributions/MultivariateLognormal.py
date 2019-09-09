import sys
sys.path.insert (0, '..')

import numpy as np
# from numpy.random import multivariate_normal as multivar_normal
from scipy.stats import multivariate_normal as multivar_normal
import scipy
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
            lognormal_values = np.array ([scipy.stats.lognorm.rvs ( \
                    scale=np.exp (mu[i]), s=np.sqrt(S[i, i]), \
                    random_state=np.random.RandomState ()) \
                    for i in range (len (mu))])
        else:        
            lognormal_values = []
            for _ in range (n):
                lognormal_values.append ([scipy.stats.lognorm.rvs ( \
                        scale=np.exp (mu[i]), s=np.sqrt(S[i, i]), \
                        size=n, random_state=np.random.RandomState ()) \
                        for i in range (len (mu))])
        return lognormal_values


    def pdf (self, x):
        """ Returns the value of the probability density function of 
            this random variable on point x. """
        # This is not very clean... we are assuming that the 
        # distribution is not concentrated around very small numbers.
        # What we should do instead is to continue calculations until we
        # can determine that the probability of the point is actually 
        # very small...
        # Now that I decided to tackle this problem I realized that 
        # there are some distributions which can have enourmeous pdf for
        # small values of x (when mean of lognormal is negative). 
        if any (xi <= 0 for xi in x):        
            return 0

        log_p = self.log_pdf (x)
        if log_p > 707:
            return 0
            return float ("inf")
        return float (np.exp (log_p))
    

    def log_pdf (self, x):
        """ Returns the value of the probability density function of 
            this random variable on point x. """
        if any (xi <= 0 for xi in x):        
            return 0

        mu = self.__mu
        n = len (mu)
        inv_S = self.__get_S_inverse ()
        det_S = self.__get_S_determinant ()
        logx_minus_mu = np.log (x) - mu
        logx_minus_mu.shape = (n, 1)
        logx_minus_mu_t = logx_minus_mu.transpose ()

        term1 = -np.log (np.sqrt (np.power (2 * np.pi, n) * \
                abs (det_S)))
        term2 = 0
        for xi in x:
            term2 -= np.log (xi)
        term3 = -.5 * np.dot (np.dot (logx_minus_mu_t, inv_S), \
                logx_minus_mu)
        return float (term1 + term2 + term3)


    @staticmethod
    def create_lognormal_with_shape (mu, S):
        """ Creates a Lognormal distribution with mean mu and 
            covariance S (diagnoal matrix). 
        
        Parameters
            mu: a array with size n. Every component of mu must be 
                greater than 1e-150.
            S: a matrix of size n x n.

        """
        # None of the components in mu can be smaller than 1e-150
        mu = np.array (mu)
        S = np.array (S)
        S_diagonal = S.diagonal ()
        mu2 = mu * mu
        n = len (mu)

        normal_mu = np.log (mu2 / np.sqrt (S_diagonal + mu2))
        normal_S_diagonal = np.log ((S_diagonal + mu2) / mu2)               
        normal_S = normal_S_diagonal * np.eye (n)
        return MultivariateLognormal (normal_mu, normal_S)
