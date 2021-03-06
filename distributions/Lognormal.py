import numpy as np
from utils import safe_log

class Lognormal:
    """ This class implements a lognormal random variable. """

    def __init__ (self, mu, sigma):
        """ Default constructor. mu and sigma are the mean and standard 
            deviation, respectively, of the underlying normal 
            distribution. """
        self.__mu = mu
        self.__s = sigma

    
    def copy (self):
        """ Returns a copy of this object. """
        cpy = Lognormal (self.__mu, self.__s)
        return cpy
        


    def mean (self):
        """ Returns the mean of this random variable. """
        mu = self.__mu
        s = self.__s
        return np.exp (mu + s * s / 2)

    def variance (self):
        """ Returns the variance of the Lognormal distribution. """
        s2 = self.__s * self.__s
        var = (np.exp (s2) - 1) * np.exp (2 * self.__mu + s2)
        return var

    def rvs (self, n=None):
        """ Returns independent observations of this random variable. 
        """
        mu = self.__mu
        s = self.__s
        if n is None:
            return np.random.lognormal (mu, s)
        else:
            return [np.random.lognormal (mu, s) for i in range (n)]


    def pdf (self, x):
        """ Returns the value of the probability density function of 
            this random variable on point x. """
        if x <= 0:
            return 0        

        mu = self.__mu
        s = self.__s
        return (1 / x) * (1 / np.sqrt (2 * np.pi * s * s)) \
                * np.exp (- (np.log (x) - mu) ** 2 / (2 * s * s))
    

    def log_pdf (self, x):
        """ Returns the log value of the probability density function of 
            this random variable on point x. """
        # TODO: simplify calculations
        return safe_log (self.pdf (x))
