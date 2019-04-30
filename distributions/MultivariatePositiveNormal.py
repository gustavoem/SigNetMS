import sys
sys.path.insert (0, '..')

import numpy as np
from numpy.random import multivariate_normal as np_multivar_normal
from utils import safe_log
from utils import safe_exp_ratio
from utils import safe_exp

class MultivariatePositiveNormal:
    """ This class implements a multivariate Normal random variable that
        is positive. This is a special case of a Multivariate Truncated
        Normal distribution.

        Warning! This class has an improper pdf function! Since we don't
        know how to compute the actual pdf, we will just return a value 
        of a function that is proportional to the pdf. """

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
        raise NotImplementedError ("This method is not implemented.")


    def rvs (self):
        """ Returns independent observations of this random variable. 
        """
        mu = self.__mu
        S = self.__S
        positive = False
        counter = 0
        while not positive:
            normal_values = np_multivar_normal (mu, S)
            if all (normal_values > 0):
                positive = True
            counter += 1
        if counter > 100:
            print ("maybe you are wasting your time :(")
        return normal_values


    def pdf (self, x):
        """ Returns the value of the probability density function of 
            this random variable on point x. """
        return safe_exp (self.log_pdf (x))
    
    def log_pdf (self, x):
        """ Returns the value of the probability density function of 
            this random variable on point x. """
        mu = self.__mu
        S = self.__S
        inv_S = self.__get_S_inverse ()
        n = len (mu)
        x_minus_mu = x - mu
        x_minus_mu.shape = (n, 1)
        x_minus_mu_t = x_minus_mu.transpose ()
        log_p = float (-.5 * np.dot (np.dot (x_minus_mu_t, inv_S), 
                x_minus_mu))
        return log_p
        

    def __normal_log_pdf (self, x):
        """ Returns the pdf value of x as if this instance were Normal. 
        """
        mu = self.__mu
        S = self.__S
        n = len (mu)
        inv_S = self.__get_S_inverse ()
        det_S = self.__get_S_determinant ()
        x_minus_mu = x - mu
        x_minus_mu.shape = (n, 1)
        x_minus_mu_t = x_minus_mu.transpose ()
        term1 = 1 / np.sqrt (np.power (2 * np.pi, n) * det_S)
        log_term2 = float (-.5 * np.dot (np.dot (x_minus_mu_t, inv_S), 
                x_minus_mu))
        return safe_log (term1) + log_term2


    def pdf_ratio (self, x, y):
        """ Calculates the ratio of the pdf. """
        logpx = self.__normal_log_pdf (x)
        logpy = self.__normal_log_pdf (y)
        ratio = safe_exp_ratio (logpx, logpy)
        return ratio
