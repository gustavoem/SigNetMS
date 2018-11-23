import numpy as np

class MultivariateLognormal:
    """ This class implements a multivariate lognormal random variable. 
    """

    def __init__ (self, mu, Sigma):
        """ Default constructor. mu and sigma are the mean and variance
            covariance matrix, respectively, of the underlying 
            multivariate normal distribution. """
        self.__mu = mu
        self.__S = sigma


    def mean (self):
        """ Returns the mean of this random variable. """
        mu = self.__mu
        S = self.__S
        return np.exp (mu + s * s / 2)


    def rvs (self, n=None):
        """ Returns independent observations of this random variable. 
        """
        return []


    def pdf (self, x):
        """ Returns the value of the probability density function of 
            this random variable on point x. """
        mu = self.__mu
        S = self.__S
