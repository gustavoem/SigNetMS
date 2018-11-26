import numpy as np

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
        cpy.__pdf_const = self.__pdf_const
        return cpy


    @staticmethod   
    def __gamma_function (a):
        """ Calculates the gamma function of an integer a. """
        if a == 1:
            return 1
        elif a == 2:
            return 1
        else:
            return (a - 1) * self.__gamma_function (a - 1)
        

    def mean (self):
        """ Returns the mean of this random variable. """
        return self.__a * self.__b


    def rvs (self, n=None):
        """ Returns independent observations of this random variable.
        """
        a = self.__a
        b = self.__b
        if n is None:
            return np.random.gamma (a, b)
        else:
            return [np.random.gamma (a, b) for i in range (n)]


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
