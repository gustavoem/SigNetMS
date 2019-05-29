import sys
sys.path.insert (0, '..')

import numpy as np
from distributions.Gamma import Gamma
from utils import safe_log


class DiscreteLaplacian:
    """ This class implements a Laplacian discrete random variable.
        A variable X is distributed as a discrete Laplacian if 
        p_i (X = j) is proportinal to e^{beta | i - j |}. X variable
        takes values from 1 to n. This implementation is specific to
        beta = 0.5. """

    def __init__ (self, N, i):
        """ Default constructor. """
        self.__N = N
        self.__i = i
        self.__constant = (np.exp (.5) - 1) \
                / (2 - (np.exp (-(N - i) / 2)  \
                + (np.exp (-(i - 1) / 2))))

    def copy (self):
        """ Return a copy of this object. """
        cpy = Gamma (self.__N, self.__i)
        return cpy


    def mean (self):
        """ Returns the mean of this random variable. """
        raise NotImplementedError


    def variance (self):
        """ Returns the variance of this random variable. """
        raise NotImplementedError
        
    
    def __randomly_variate (self):
        """ Returns one observation of the random variable. """
        N = self.__N
        i = self.__i
        if i == 1:
            signal = 1
        elif i == N:
            signal = -1
        else:
            signal = np.random.choice ([1, -1])
        j = np.random.geometric (np.exp (- 1 / 2))
        if signal > 0:
            if i + j > N:
                return N
            else:
                return i + j
        else:
            if i - j < 1:
                return 1
            else:
                return i - j


    def rvs (self, n=None):
        """ Returns independend observations of this random variable.
        """
        if n is None:
            return self.__randomly_variate ()
        else:
            return [self.__randomly_variate () for x in range (n)]

    def pdf (self, j):
        """ Returns p_i (j). """
        i = self.__i
        N = self.__N
        if i == j:
            return 0
        elif j > N:
            return 0
        elif j <= 0:
            return 0 

        abs_diff = abs (i - j)
        return self.__constant * np.exp (-.5 * abs_diff)


    def log_pdf (self, j):
        """ Returns log p_j (j). """
        # TODO: simplify calculations
        return safe_log (self.pdf (j))
