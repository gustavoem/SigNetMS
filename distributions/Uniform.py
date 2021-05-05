import numpy as np
from utils import safe_log

class Uniform:
    """ This class implements a uniform random variable. """

    def __init__ (self, a, b):
        """ Default constructor. a and b are intervals used to construc
            the non-zero probability space. This space is the interval
            [a, b). """
        if b <= a:
            raise ValueError ("invalid uniform range, b should be" 
                +  "greater than a")

        self.__a = a
        self.__b = b

    
    def copy (self):
        """ Returns a copy of this object """
        cpy = Uniform (self.__a, self.__b)
        return cpy


    def mean (self):
        a, b = self.__a, self.__b
        return (a + b) / 2


    def variance (self):
        a, b = self.__a, self.__b
        return (b - a) * (b - a) / 12


    def rvs (self, n=None):
        if n is None:
            return np.random.uniform (self.__a, self.__b)
        else:
            return [np.random.uniform (self.__a, self.__b) for i in \
                    range (n)]


    def get_a (self):
        return self.__a


    def get_b (self):
        return self.__b


    def pdf (self, x):
        a, b = self.__a, self.__b
        if x >= a and x < b:
            return 1 / (b - a)
        else:
            return 0


    def log_pdf (self, x):
        return safe_log (self.pdf (x))
