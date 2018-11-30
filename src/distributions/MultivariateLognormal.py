import numpy as np
from numpy.random import multivariate_normal as multivar_normal

class MultivariateLognormal:
    """ This class implements a multivariate lognormal random variable. 
    """

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
        cpy = MultivariateLognormal (self.__mu, self.__S)
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
        if any (xi <= 0 for xi in x):        
            return 0
            
        print ("\t\tLognorm of mean: " + str (self.__mu))
        print ("\t\tAnd variance of: " + str (self.__S.diagonal ()))

        mu = self.__mu
        S = self.__S
        n = len (mu)
        inv_S = self.__get_S_inverse ()
        det_S = self.__get_S_determinant ()
        logx_minus_mu = np.log (x) - mu
        logx_minus_mu.shape = (n, 1)
        logx_minus_mu_t = logx_minus_mu.transpose ()
    
        term1 = 1 / np.sqrt (np.power (2 * np.pi, n) * det_S)
        term2 = 1 
        for xi in x:
            term2 *= xi
        term2 = 1 / term2


        term3 = float (np.exp (-.5 * np.dot (np.dot (logx_minus_mu_t, 
            inv_S), logx_minus_mu)))

        #print ("\t\tlogx - mu = " + str(logx_minus_mu))
        #print ("\t\tinv_S = " + str(inv_S))
        #print ("\t\tS = " + str(S))
        #print ("\t\tterm3: " + str(term3))
        
        return term1 * term2 * term3

