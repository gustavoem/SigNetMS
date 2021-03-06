import numpy as np
from scipy.stats import gamma as gamma
from utils import safe_log


class Gamma:
    """ This class implements a gamma random variable. """
    
    def __init__ (self, a, b):
        """ Default constructor. The parameter a is the alpha parameter
            and b is the scale parameter. The parameter a should be 
            an integer. """
        if int (a) - a != 0:
            print ("WARNING! Wrong usage of Gamma function. The " +
                    " parameter a should be an integer!")
        a = int (a)
        self.__a = a
        self.__b = b
        gamma_a = self.__gamma_function (a)
        self.__pdf_const = 1 / (gamma_a * (b ** a))


    def copy (self):
        """ Return a copy of this object. """
        cpy = Gamma (self.__a, self.__b)
        return cpy


    @staticmethod   
    def __gamma_function (a):
        """ Calculates the gamma function of an integer a. """
        if a == 1:
            return 1
        elif a == 2:
            return 1
        else:
            return (a - 1) * Gamma.__gamma_function (a - 1)
        

    def mean (self):
        """ Returns the mean of this random variable. """
        return self.__a * self.__b


    def variance (self):
        """ Returns the variance of this random variable. """
        return self.__a * self.__b * self.__b


    def rvs (self, n=None):
        """ Returns independent observations of this random variable.
        """
        a = self.__a
        b = self.__b
        if n is None:
            return gamma.rvs (a, scale=b)
        else:
            return [gamma.rvs (a, scale=b) for i in range (n)]


    def get_a (self):
        """ Retuns the parameter a. """
        return self.__a
    
    
    def get_b (self):
        """ Retuns the parameter b. """
        return self.__b


    def pdf (self, x):
        """ Returns the value of the probability density function of 
            this random variable on point x. """
        if x <= 0:
            return 0
    
        a = self.__a
        b = self.__b
        term1 = pow (x, a - 1)
        term2 = np.exp (- x / b)
        return self.__pdf_const * term1 * term2


    def log_pdf (self, x):
        """ Returns the log value of the probability density function of 
            this random variable on point x. """
        # TODO: simplify calculations
        return safe_log (self.pdf (x))
